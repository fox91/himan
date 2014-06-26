/**
 * @file hybrid_pressure.cpp
 *
 *  @date: Mar 23, 2013
 *  @author aaltom
 */

#include "hybrid_pressure.h"
#include "logger_factory.h"
#include <boost/lexical_cast.hpp>
#include "level.h"
#include "forecast_time.h"

using namespace std;
using namespace himan::plugin;

hybrid_pressure::hybrid_pressure()
{
	// Vertkoord_A and Vertkoord_B refer to full hybrid-level coefficients
	itsClearTextFormula = "P = Vertkoord_A + P0 * Vertkoord_B";

	itsLogger = logger_factory::Instance()->GetLog("hybrid_pressure");

}

void hybrid_pressure::Process(std::shared_ptr<const plugin_configuration> conf)
{
	Init(conf);

	SetParams({param("P-HPA", 1, 0, 3, 0)});

	Start();

}

/*
 * Calculate()
 *
 * This function does the actual calculation.
 */

void hybrid_pressure::Calculate(shared_ptr<info> myTargetInfo, unsigned short theThreadIndex)
{

	const param PParam("P-PA");
	const param QParam("Q-KGKG");
	const level PLevel(himan::kHeight, 0, "HEIGHT");

	auto myThreadedLogger = logger_factory::Instance()->GetLog("hybrid_pressureThread #" + boost::lexical_cast<string> (theThreadIndex));

	forecast_time forecastTime = myTargetInfo->Time();
	level forecastLevel = myTargetInfo->Level();

	myThreadedLogger->Info("Calculating time " + static_cast<string>(*forecastTime.ValidDateTime()) + " level " + static_cast<string> (forecastLevel));

	info_t PInfo = Fetch(forecastTime, PLevel, PParam, false);
	info_t QInfo = Fetch(forecastTime, PLevel, QParam, false);

	if (!PInfo || !QInfo)
	{
		myThreadedLogger->Warning("Skipping step " + boost::lexical_cast<string> (forecastTime.Step()) + ", level " + static_cast<string> (forecastLevel));
		return;
	}
	
	SetAB(myTargetInfo, QInfo);

	/* 
	 * Vertical coordinates for full hybrid levels.
	 * For Hirlam data, coefficients A and B are already interpolated to full level coefficients in the grib-file.
	 * For Harmonie and ECMWF interpolation is done, when reading data from the grib-file. (NFmiGribMessage::PV)
	 */

	std::vector<double> ab = QInfo->Grid()->AB();

   	double A = ab[0];
   	double B = ab[1];

	LOCKSTEP(myTargetInfo, PInfo)
	{
		double P = PInfo->Value();
		
		if (P == kFloatMissing)
		{
			continue;
		}

		double hybrid_pressure = 0.01 * (A + P * B);

		myTargetInfo->Value(hybrid_pressure);
	}


	myThreadedLogger->Info("[CPU] Missing values: " + boost::lexical_cast<string> (myTargetInfo->Data()->MissingCount()) + "/" + boost::lexical_cast<string> (myTargetInfo->Data()->Size()));

}
