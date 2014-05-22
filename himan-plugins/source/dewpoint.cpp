/**
 * @file dewpoint.cpp
 *
 * @date Jan 21, 2012
 * @author partio
 */

#include "dewpoint.h"
#include "plugin_factory.h"
#include "logger_factory.h"
#include <boost/lexical_cast.hpp>
#include "NFmiGrid.h"

#define HIMAN_AUXILIARY_INCLUDE

#include "fetcher.h"

#undef HIMAN_AUXILIARY_INCLUDE

using namespace std;
using namespace himan::plugin;

const double Rw_div_L = himan::constants::kRw / himan::constants::kL;

dewpoint::dewpoint()
{
	itsClearTextFormula = "Td = T / (1 - (T * ln(RH)*(Rw/L)))";
	itsCudaEnabledCalculation = true;

	itsLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("dewpoint"));

}

void dewpoint::Process(shared_ptr<const plugin_configuration> conf)
{

	Init(conf);
	
	/*
	 * Set target parameter to dewpoint.
	 *
	 */

	vector<param> params;

	param requestedParam ("TD-C", 10);
	requestedParam.Unit(kK);
	
	// GRIB 2

	requestedParam.GribDiscipline(0);
	requestedParam.GribCategory(0);
	requestedParam.GribParameter(6);

	params.push_back(requestedParam);

	SetParams(params);

	Start();
}

/*
 * Calculate()
 *
 * This function does the actual calculation.
 */

void dewpoint::Calculate(shared_ptr<info> myTargetInfo, unsigned short threadIndex)
{

	
	shared_ptr<fetcher> f = dynamic_pointer_cast <fetcher> (plugin_factory::Instance()->Plugin("fetcher"));

	// Required source parameters

	param TParam("T-K");
	param RHParam("RH-PRCNT");

	unique_ptr<logger> myThreadedLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("dewpointThread #" + boost::lexical_cast<string> (threadIndex)));
	
	ResetNonLeadingDimension(myTargetInfo);

	myTargetInfo->FirstParam();

	bool useCudaInThisThread = compiled_plugin_base::GetAndSetCuda(threadIndex);
	
	while (AdjustNonLeadingDimension(myTargetInfo))
	{

		myThreadedLogger->Debug("Calculating time " + myTargetInfo->Time().ValidDateTime()->String("%Y%m%d%H%M") +
								" level " + boost::lexical_cast<string> (myTargetInfo->Level().Value()));

		double TBase = 0;
		double RHScale = 1;
		
		shared_ptr<info> TInfo;
		shared_ptr<info> RHInfo;

		try
		{
			TInfo = f->Fetch(itsConfiguration,
								myTargetInfo->Time(),
								myTargetInfo->Level(),
								TParam,
								itsConfiguration->UseCudaForPacking() && useCudaInThisThread);

			RHInfo = f->Fetch(itsConfiguration,
								myTargetInfo->Time(),
								myTargetInfo->Level(),
								RHParam,
								itsConfiguration->UseCudaForPacking() && useCudaInThisThread);


		}
		catch (HPExceptionType& e)
		{
			switch (e)
			{
				case kFileDataNotFound:
					itsLogger->Info("Skipping step " + boost::lexical_cast<string> (myTargetInfo->Time().Step()) + ", level " + boost::lexical_cast<string> (myTargetInfo->Level().Value()));
					myTargetInfo->Data()->Fill(kFloatMissing); // Fill data with missing value

					if (itsConfiguration->StatisticsEnabled())
					{
						itsConfiguration->Statistics()->AddToMissingCount(myTargetInfo->Grid()->Size());
						itsConfiguration->Statistics()->AddToValueCount(myTargetInfo->Grid()->Size());
					}

					continue;
					break;

				default:
					throw runtime_error(ClassName() + ": Unable to proceed");
					break;
			}
		}

		assert(TInfo->Grid()->AB() == RHInfo->Grid()->AB());
		
		SetAB(myTargetInfo, TInfo);

		if (RHInfo->Param().Unit() != kPrcnt)
		{
			itsLogger->Warning("Unable to determine RH unit, assuming percent");
		}

		// Formula assumes T == Celsius

		if (TInfo->Param().Unit() == kK)
		{
			TBase = -himan::constants::kKelvin;
		}

		size_t missingCount = 0;
		size_t count = 0;

		bool equalGrids = (*myTargetInfo->Grid() == *TInfo->Grid() && *myTargetInfo->Grid() == *RHInfo->Grid());

		string deviceType;

#ifdef HAVE_CUDA

		// If we read packed data but grids are not equal we cannot use cuda
		// for calculations (our cuda routines do not know how to interpolate)

		if (!equalGrids && (TInfo->Grid()->IsPackedData() || RHInfo->Grid()->IsPackedData()))
		{
			myThreadedLogger->Debug("Unpacking for CPU calculation");

			Unpack({TInfo, RHInfo});
		}

		if (useCudaInThisThread && equalGrids)
		{

			deviceType = "GPU";

			auto opts = CudaPrepare(myTargetInfo, TInfo, RHInfo); 

			dewpoint_cuda::Process(*opts);

			missingCount = opts->missing;
			count = opts->N;

			CudaFinish(move(opts), myTargetInfo, TInfo, RHInfo);

		}
		else
#endif
		{
			deviceType = "CPU";

			shared_ptr<NFmiGrid> targetGrid(myTargetInfo->Grid()->ToNewbaseGrid());
			shared_ptr<NFmiGrid> TGrid(TInfo->Grid()->ToNewbaseGrid());
			shared_ptr<NFmiGrid> RHGrid(RHInfo->Grid()->ToNewbaseGrid());

			assert(targetGrid->Size() == myTargetInfo->Data()->Size());

			myTargetInfo->ResetLocation();

			targetGrid->Reset();
			
			while (myTargetInfo->NextLocation() && targetGrid->Next())
			{
				count++;

				double T = kFloatMissing;
				double RH = kFloatMissing;

				InterpolateToPoint(targetGrid, TGrid, equalGrids, T);
				InterpolateToPoint(targetGrid, RHGrid, equalGrids, RH);

				if (T == kFloatMissing || RH == kFloatMissing)
				{
					missingCount++;

					myTargetInfo->Value(kFloatMissing);
					continue;
				}

				T += TBase;
				RH *= RHScale;

				double TD = kFloatMissing;

				if (RH > 50)
				{
					TD = T - ((100 - RH) * 0.2) + constants::kKelvin;
				}
				else
				{
					TD = T / (1 - (T * log(RH) * (Rw_div_L))) + constants::kKelvin;
				}

				if (!myTargetInfo->Value(TD))
				{
					throw runtime_error(ClassName() + ": Failed to set value to matrix");
				}

			}

			/*
			 * Newbase normalizes scanning mode to bottom left -- if that's not what
			 * the target scanning mode is, we have to swap the data back.
			 */

			SwapTo(myTargetInfo, kBottomLeft);
		}

		if (itsConfiguration->StatisticsEnabled())
		{
			itsConfiguration->Statistics()->AddToMissingCount(missingCount);
			itsConfiguration->Statistics()->AddToValueCount(count);
		}
		
		/*
		 * Now we are done for this level
		 *
		 * Clone info-instance to writer since it might change our descriptor places
		 */

		myThreadedLogger->Info("[" + deviceType + "] Missing values: " + boost::lexical_cast<string> (missingCount) + "/" + boost::lexical_cast<string> (count));

		if (itsConfiguration->FileWriteOption() != kSingleFile)
		{
			WriteToFile(myTargetInfo);
		}
	}
}

#ifdef HAVE_CUDA

unique_ptr<dewpoint_cuda::options> dewpoint::CudaPrepare(shared_ptr<info> myTargetInfo, shared_ptr<info> TInfo, shared_ptr<info> RHInfo)
{
	unique_ptr<dewpoint_cuda::options> opts(new dewpoint_cuda::options);

	opts->t = TInfo->ToSimple();
	opts->rh = RHInfo->ToSimple();
	opts->td = myTargetInfo->ToSimple();

	opts->N = opts->td->size_x * opts->td->size_y;

	if (TInfo->Param().Unit() == kK)
	{
		opts->t_base = -himan::constants::kKelvin;
	}

	if (RHInfo->Param().Unit() != kPrcnt)
	{
		// If unit cannot be detected, assume the values are from 0 .. 1
		opts->rh_scale = 100;
	}

	return opts;
}

void dewpoint::CudaFinish(unique_ptr<dewpoint_cuda::options> opts, shared_ptr<info> myTargetInfo, shared_ptr<info> TInfo, shared_ptr<info> RHInfo)
{
	// Copy data back to infos

	CopyDataFromSimpleInfo(myTargetInfo, opts->td, false);

	assert(TInfo->Grid()->ScanningMode() == RHInfo->Grid()->ScanningMode());

	if (TInfo->Grid()->IsPackedData())
	{
		CopyDataFromSimpleInfo(TInfo, opts->t, true);
	}

	if (RHInfo->Grid()->IsPackedData())
	{
		CopyDataFromSimpleInfo(RHInfo, opts->rh, false);
	}

	SwapTo(myTargetInfo, TInfo->Grid()->ScanningMode());

}

#endif