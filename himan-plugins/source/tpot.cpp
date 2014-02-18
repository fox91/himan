/**
 * @file tpot.cpp
 *
 * @brief Plugin to calculate potential temperature, pseudo-adiabatic
 * potential temperature or equivalent potential temperature.
 *
 * Created on: Nov 20, 2012
 * @author partio
 */

#include "tpot.h"
#include "plugin_factory.h"
#include "logger_factory.h"
#include "timer_factory.h"
#include <boost/lexical_cast.hpp>

#define HIMAN_AUXILIARY_INCLUDE

#include "fetcher.h"

#undef HIMAN_AUXILIARY_INCLUDE

using namespace std;
using namespace himan::plugin;

#include "tpot_cuda.h"
#include "cuda_helper.h"
#include "util.h"

tpot::tpot()
: itsThetaCalculation(false)
, itsThetaWCalculation(false)
, itsThetaECalculation(false)
{
	itsClearTextFormula = "TP = Tk * pow((1000/P), 0.286) ; TPW calculated with LCL ; TPE = X"; 
	itsCudaEnabledCalculation = true;

	itsLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("tpot"));

}

void tpot::Process(std::shared_ptr<const plugin_configuration> conf)
{
	Init(conf);

	/*
	 * Set target parameter to potential temperature
	 * - name TP-K
	 * - univ_id 8
	 * - grib2 descriptor 0'00'002
	 *
	 * We need to specify grib and querydata parameter information
	 * since we don't know which one will be the output format.
	 * (todo: we could check from conf but why bother?)
	 *
	 */

	vector<param> theParams;

	if (itsConfiguration->Exists("theta") && itsConfiguration->GetValue("theta") == "true")
	{
		itsThetaCalculation = true;

		itsLogger->Trace("Theta calculation requested");

		param p("TP-K", 8);

		p.GribDiscipline(0);
		p.GribCategory(0);
		p.GribParameter(2);		
		
		theParams.push_back(p);
	}

	if (itsConfiguration->Exists("thetaw") && itsConfiguration->GetValue("thetaw") == "true")
	{
		itsThetaWCalculation = true;

		itsLogger->Trace("ThetaW calculation requested");

		param p("TPW-K", 9);

		// Sharing number with thetae!

		p.GribDiscipline(0);
		p.GribCategory(0);
		p.GribParameter(3);

		theParams.push_back(p);
	}

	if (itsConfiguration->Exists("thetae") && itsConfiguration->GetValue("thetae") == "true")
	{
		itsThetaECalculation = true;

		itsLogger->Trace("ThetaE calculation requested");

		param p("TPE-K", 99999);

		// Sharing number with thetaw!

		p.GribDiscipline(0);
		p.GribCategory(0);
		p.GribParameter(3);

		theParams.push_back(p);

	}

	if (theParams.size() == 0)
	{
		// By default assume we'll calculate theta

		itsThetaCalculation = true;

		param p("TP-K", 8);

		p.GribDiscipline(0);
		p.GribCategory(0);
		p.GribParameter(2);

		theParams.push_back(p);
	}

	SetParams(theParams);

	Start();

}

/*
 * Calculate()
 *
 * This function does the actual calculation.
 */

void tpot::Calculate(shared_ptr<info> myTargetInfo, unsigned short threadIndex)
{

	shared_ptr<fetcher> theFetcher = dynamic_pointer_cast <fetcher> (plugin_factory::Instance()->Plugin("fetcher"));

	unique_ptr<logger> myThreadedLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("tpotThread #" + boost::lexical_cast<string> (threadIndex)));

	ResetNonLeadingDimension(myTargetInfo);

	myTargetInfo->FirstParam();

	params PParam = { param("P-PA"), param("P-HPA") };
	
	bool useCudaInThisThread = compiled_plugin_base::GetAndSetCuda(itsConfiguration, threadIndex);

	while (AdjustNonLeadingDimension(myTargetInfo))
	{

		myThreadedLogger->Debug("Calculating time " + myTargetInfo->Time().ValidDateTime()->String("%Y%m%d%H") +
								" level " + boost::lexical_cast<string> (myTargetInfo->Level().Value()));

		double PScale = 1;
		double TBase = 0;
		double TDBase  = 0;
		
		// Source infos
		shared_ptr<info> TInfo;
		shared_ptr<info> PInfo;
		shared_ptr<info> TDInfo;

		shared_ptr<NFmiGrid> PGrid;
		shared_ptr<NFmiGrid> TDGrid;

		bool isPressureLevel = (myTargetInfo->Level().Type() == kPressure);

		try
		{

			TInfo = theFetcher->Fetch(itsConfiguration,
										myTargetInfo->Time(),
										myTargetInfo->Level(),
										param("T-K"),
										itsConfiguration->UseCudaForPacking() && useCudaInThisThread);

			if (!isPressureLevel)
			{
				// Source info for P
				PInfo = theFetcher->Fetch(itsConfiguration,
											myTargetInfo->Time(),
											myTargetInfo->Level(),
											PParam,
											itsConfiguration->UseCudaForPacking() && useCudaInThisThread);

				if (PInfo->Param().Unit() == kHPa || PInfo->Param().Name() == "P-HPA")
				{
					PScale = 100;
				}

				PGrid = shared_ptr<NFmiGrid> (PInfo->Grid()->ToNewbaseGrid());
			}

			if (itsThetaWCalculation || itsThetaECalculation)
			{
				TDInfo = theFetcher->Fetch(itsConfiguration,
										myTargetInfo->Time(),
										myTargetInfo->Level(),
										param("TD-C"),
										itsConfiguration->UseCudaForPacking() && useCudaInThisThread);

				TDGrid = shared_ptr<NFmiGrid> (TDInfo->Grid()->ToNewbaseGrid());

			}
		}
		catch (HPExceptionType e)
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

		assert(isPressureLevel || ((PInfo->Grid()->AB() == TInfo->Grid()->AB() && PInfo->Grid()->AB() == TDInfo->Grid()->AB())));

		SetAB(myTargetInfo, TInfo);

		if (TInfo->Param().Unit() == kC)
		{
			TBase = himan::constants::kKelvin;
		}

		if (TDInfo && TDInfo->Param().Unit() == kC)
		{
			TBase = himan::constants::kKelvin;
		}

		shared_ptr<NFmiGrid> targetGrid(myTargetInfo->Grid()->ToNewbaseGrid());
		shared_ptr<NFmiGrid> TGrid(TInfo->Grid()->ToNewbaseGrid());

		size_t missingCount = 0;
		size_t count = 0;

		assert(targetGrid->Size() == myTargetInfo->Data()->Size());

		bool equalGrids = (*myTargetInfo->Grid() == *TInfo->Grid() && (isPressureLevel || *myTargetInfo->Grid() == *PInfo->Grid()));

		string deviceType;

#ifdef HAVE_CUDA

		// If we read packed data but grids are not equal we cannot use cuda
		// for calculations (our cuda routines do not know how to interpolate)

		if (!equalGrids && (TInfo->Grid()->IsPackedData() || TDInfo->Grid()->IsPackedData() || (PInfo && PInfo->Grid()->IsPackedData())))
		{
			myThreadedLogger->Debug("Unpacking for CPU calculation");

			Unpack({TInfo, TDInfo});

			if (PInfo)
			{
				Unpack({PInfo});
			}
		}

		if (useCudaInThisThread && equalGrids && (itsThetaCalculation && !itsThetaWCalculation && !itsThetaECalculation))
		{
			deviceType = "GPU";

			auto opts = CudaPrepare(myTargetInfo, TInfo, PInfo, TDInfo);

			tpot_cuda::Process(*opts);

			count = opts->N;
			missingCount = opts->missing;

			CudaFinish(move(opts), myTargetInfo, TInfo, PInfo, TDInfo);

		}
		else
#endif
		{

			deviceType = "CPU";

			myTargetInfo->ResetLocation();

			targetGrid->Reset();

			while (myTargetInfo->NextLocation() && targetGrid->Next())
			{
				count++;

				double T = kFloatMissing;
				double P = kFloatMissing;
				double TD = kFloatMissing;

				InterpolateToPoint(targetGrid, TGrid, equalGrids, T);

				if (isPressureLevel)
				{
					P = myTargetInfo->Level().Value() * 100;
				}
				else
				{
					InterpolateToPoint(targetGrid, PGrid, equalGrids, P);
				}

				if (itsThetaWCalculation || itsThetaECalculation)
				{
					InterpolateToPoint(targetGrid, TDGrid, equalGrids, TD);
				}

				if (T == kFloatMissing || P == kFloatMissing || ((itsThetaECalculation || itsThetaWCalculation) && TD == kFloatMissing))
				{
					missingCount++;

					myTargetInfo->Value(kFloatMissing);
					continue;
				}

				T -= TBase; // to Kelvin
				TD -= TDBase; // to Kelvin
				P *= PScale; // to Pa

				double value = kFloatMissing;
				double theta = kFloatMissing;

				if (itsThetaCalculation)
				{
					theta = Theta(P, T);

					myTargetInfo->Param(param("TP-K"));

					if (!myTargetInfo->Value(theta))
					{
						throw runtime_error(ClassName() + ": Failed to set value to matrix");
					}
				}
				
				if (itsThetaWCalculation)
				{
					
					value = ThetaW(P, T, TD);
					myTargetInfo->Param(param("TPW-K"));

					if (!myTargetInfo->Value(value))
					{
						throw runtime_error(ClassName() + ": Failed to set value to matrix");
					}
				}
				
				if (itsThetaECalculation)
				{
					value = ThetaE(P, T, TD, theta);
					
					myTargetInfo->Param(param("TPE-K"));

					if (!myTargetInfo->Value(value))
					{
						throw runtime_error(ClassName() + ": Failed to set value to matrix");
					}
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
		 *
		 */

		myThreadedLogger->Info("Missing values: " + boost::lexical_cast<string> (missingCount) + "/" + boost::lexical_cast<string> (count));

		if (itsConfiguration->FileWriteOption() != kSingleFile)
		{
			WriteToFile(myTargetInfo);
		}

	}
}

double tpot::Theta(double P, double T)
{

	double value = T * pow((1000 / (P*0.01)), 0.28586);

	return value;
}

double tpot::ThetaW(double P, double T, double TD)
{
   const double Pstep = 500; // Pa
   double value = kFloatMissing;

   // Search LCL level

   vector<double> LCL = util::LCL(P, T, TD);

   double Pint = LCL[0]; // Pa
   double Tint = LCL[1]; // K

   int i = 0;

   if (Tint == kFloatMissing || Pint == kFloatMissing)
   {
	   return kFloatMissing;
   }
   else
   {
	   /*
		* Units: Temperature in Kelvins, Pressure in Pascals
		*/

	   double T0 = Tint;

	   double Z = kFloatMissing;

	   while (++i < 500) // usually we don't reach this value
	   {
		   double TA = Tint;

		   if (i <= 2)
		   {
			   Z = i * Pstep/2;
		   }
		   else
		   {
			   Z = 2 * Pstep;
		   }

		   // Gammas() takes Pa
		   Tint = T0 + util::Gammas(Pint, Tint) * Z;

		   if (i > 2)
		   {
			   T0 = TA;
		   }

		   Pint += Pstep;

		   if (Pint >= 1e5)
		   {
			   value = Tint;
			   break;
		   }
	   }
   }

   assert(value == value); // check NaN

   return value;
}

double tpot::ThetaE(double P, double T, double TD, double theta)
{
	vector<double> LCL = util::LCL(P, T, TD);

	double TLCL = LCL[1];

	if (theta == kFloatMissing)
	{
		// theta was not calculated in this plugin session :(

		theta = Theta(P, T);
	}

	double Es = util::Es(TLCL);
	double ZQs = himan::constants::kEp * (Es / (P - Es));

	double value = theta * exp(himan::constants::kL * ZQs / himan::constants::kCp / (TLCL));

	return value;

}

#ifdef HAVE_CUDA

unique_ptr<tpot_cuda::options> tpot::CudaPrepare(shared_ptr<info> myTargetInfo, shared_ptr<info> TInfo, shared_ptr<info> PInfo, shared_ptr<info> TDInfo)
{
	unique_ptr<tpot_cuda::options> opts(new tpot_cuda::options);

	opts->is_constant_pressure = (myTargetInfo->Level().Type() == kPressure);

	opts->t = TInfo->ToSimple();
	opts->tp = myTargetInfo->ToSimple();

	opts->theta = itsThetaCalculation;
	opts->thetaw = itsThetaWCalculation;
	opts->thetae = itsThetaECalculation;

	if (!opts->is_constant_pressure)
	{
		opts->p = PInfo->ToSimple();

		if (PInfo->Param().Unit() == kPa || PInfo->Param().Name() == "P-PA")
		{
			opts->p_scale = 100;
		}
	}
	else
	{
		opts->p_const = myTargetInfo->Level().Value() * 100; // Pa
	}

	if (TDInfo)
	{
		opts->td = TDInfo->ToSimple();

		if (TDInfo->Param().Unit() == kC)
		{
			opts->td_base = himan::constants::kKelvin;
		}
	}

	opts->N = opts->t->size_x * opts->t->size_y;

	if (TInfo->Param().Unit() == kC)
	{
		opts->t_base = himan::constants::kKelvin;
	}

	return opts;
}

void tpot::CudaFinish(unique_ptr<tpot_cuda::options> opts, shared_ptr<info> myTargetInfo, shared_ptr<info> TInfo, shared_ptr<info> PInfo, shared_ptr<info> TDInfo)
{
	// Copy data back to infos

	myTargetInfo->Data()->Set(opts->tp->values, opts->N);
	opts->tp->free_values();

	SwapTo(myTargetInfo, TInfo->Grid()->ScanningMode());

	if (TInfo->Grid()->IsPackedData())
	{
		TInfo->Data()->Set(opts->t->values, opts->N);
		TInfo->Grid()->PackedData()->Clear();
		opts->t->free_values();
	}

	if (PInfo && PInfo->Grid()->IsPackedData())
	{
		PInfo->Data()->Set(opts->p->values, opts->N);
		PInfo->Grid()->PackedData()->Clear();
		opts->p->free_values();
	}

	if (TDInfo && TDInfo->Grid()->IsPackedData())
	{
		TDInfo->Data()->Set(opts->p->values, opts->N);
		TDInfo->Grid()->PackedData()->Clear();
		opts->td->free_values();
	}
	
	// opts is destroyed after leaving this function

}

#endif