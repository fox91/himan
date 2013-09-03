/**
 *
 * @file compiled_plugin_base.cpp
 *
 * @date Jan 15, 2013
 * @author partio
 */

#include "compiled_plugin_base.h"
#include <boost/thread.hpp>
#include "plugin_factory.h"

#define HIMAN_AUXILIARY_INCLUDE

#include "neons.h"
#include "writer.h"

#undef HIMAN_AUXILIARY_INCLUDE

using namespace std;
using namespace himan::plugin;

const unsigned int MAX_THREADS = 12; //<! Max number of threads we allow
const double kInterpolatedValueEpsilon = 0.00001; //<! Max difference between two grid points (if smaller, points are considered the same)

mutex itsAdjustDimensionMutex;

unsigned short compiled_plugin_base::ThreadCount(short userThreadCount) const
{
    unsigned int coreCount = boost::thread::hardware_concurrency(); // Number of cores

    unsigned short threadCount = MAX_THREADS;

    if (userThreadCount > 0)
    {
    	threadCount = userThreadCount;
    }
    else if (MAX_THREADS > coreCount)
    {
    	threadCount = static_cast<unsigned short> (coreCount);
    }

    return threadCount;
}


bool compiled_plugin_base::InterpolateToPoint(shared_ptr<const NFmiGrid> targetGrid, shared_ptr<NFmiGrid> sourceGrid, bool gridsAreEqual, double& value)
{

	/*
	 * Logic of interpolating values:
	 *
	 * 1) If source and target grids are equal, meaning that the grid AND the area
	 *    properties are effectively the same, do not interpolate. Instead return
	 *    the value of the source grid point that matches the ordering number of the
	 *    target grid point (ie. target grid point #1 --> source grid point #1 etc).
	 *
	 * 2) If actual interpolation is needed, first get the *grid* coordinates of the
	 *    latlon target point. Then check if those grid coordinates are very close
	 *    to a grid point -- if so, return the value of the grid point. This serves two
	 *    purposes:
	 *    - We don't need to interpolate if the distance between requested grid point
	 *      and actual grid point is small enough, saving some CPU cycles
	 *    - Sometimes when the requested grid point is close to grid edge, floating
	 *      point inaccuracies might move it outside the grid. If this happens, the
	 *      interpolation fails even though initially the grid point is valid.
	 *
	 * 3) If requested source grid point is not near and actual grid point, interpolate
	 *    the value of the point.
	 */

	// Step 1)

	if (gridsAreEqual)
	{
		value = sourceGrid->FloatValue(targetGrid->GridPoint());
		return true;
	}

	const NFmiPoint targetLatLonPoint = targetGrid->LatLon();
	const NFmiPoint sourceGridPoint = targetGrid->LatLonToGrid(targetLatLonPoint.X(), targetLatLonPoint.Y());

	// Step 2)

	bool noInterpolation = (fabs(sourceGridPoint.X() - round(sourceGridPoint.X())) < kInterpolatedValueEpsilon &&
		 fabs(sourceGridPoint.Y() - round(sourceGridPoint.Y())) < kInterpolatedValueEpsilon);

	if (noInterpolation)
	{
		value = sourceGrid->FloatValue(sourceGridPoint);
		//cout << value << endl;
		return true;
	}
	

	// Step 3)

	return sourceGrid->InterpolateToLatLonPoint(targetLatLonPoint, value);
	//return InterpolateToHilPP(sourceGrid, targetLatLonPoint, value);

}
/*
bool compiled_plugin_base::InterpolateToHilPP(shared_ptr<NFmiGrid> sourceGrid, const NFmiPoint &theLatLonPoint, double &value)
{
	double x = theLatLonPoint.X();
	double y = theLatLonPoint.Y();
	if(x > sourceGrid->OriginalXNumber()-1 || y > sourceGrid->OriginalYNumber()-1)
	{
		value = kFloatMissing;
		return false;
	}
  
  	if(!BiLinearInterpolation(sourceGrid, x, y, value))
	{
		value = kFloatMissing;
		return false;
	}
	return true;
}

bool compiled_plugin_base::BiLinearInterpolation(shared_ptr<NFmiGrid> sourceGrid, double x, double y, double &theValue)
{
	double X,Y;
	
	double kFmiEps = 0.0001;
	double dx = modf(x, &X);
	double dy = modf(y, &Y);
   	
   	if ( abs(x - X) < kFmiEps)
      	x = x - kFmiEps;

   	if ( abs(y - Y) < kFmiEps)
      	y = y - kFmiEps;

   	if ( abs(x - 1.0) < kFmiEps)
      	x = x - kFmiEps;

    if ( abs(y - 1.0) < kFmiEps)
      	y = y - kFmiEps;

	dx = x - float(int(x));
    dy = y - float(int(y));
	
	if( (dx != 0.) && (dx < kFmiEps)) { dx = 0;}
	if( (dy != 0.) && (dy < kFmiEps)) { dy = 0;}
	
	if(dx > (1-kFmiEps)) { dx = 0; X++; }
	if(dy > (1-kFmiEps)) { dy = 0; Y++; }
	
	double topLeftValue;
	double topRightValue;
	double bottomLeftValue;
	double bottomRightValue;
	
	if((X == sourceGrid->getPreviousX()) && (Y == sourceGrid->getPreviousY()))
	{
	  // Still sitting within the same grid cell -->
	  // previous cell corner values can be used
	  
	  topLeftValue = sourceGrid->getTopLeftValue();
	  topRightValue = sourceGrid->getTopRightValue();
	  bottomLeftValue = sourceGrid->getBottomLeftValue();
	  bottomRightValue = sourceGrid->getBottomRightValue();
	}
  	else
	{
		if(!sourceGrid->GridPoint(X,Y))
		{
			theValue = static_cast<double>(kFloatMissing);
			
			sourceGrid->setPreviousX(kMinDouble);
			sourceGrid->setPreviousY(kMinDouble);
			
			return false;						
		}
	  
	  	topLeftValue = sourceGrid->FloatValue(0, 1);
	  	topRightValue = sourceGrid->FloatValue(1, 1);
	  	bottomLeftValue = sourceGrid->FloatValue(0, 0);
	  	bottomRightValue = sourceGrid->FloatValue(1, 0);
	  
	  	sourceGrid->setTopLeftValue(topLeftValue);
	  	sourceGrid->setTopRightValue(topRightValue);
	  	sourceGrid->setBottomLeftValue(bottomLeftValue);
	  	sourceGrid->setBottomRightValue(bottomRightValue);
	  
	  	sourceGrid->setPreviousX(X);
	  	sourceGrid->setPreviousY(Y);
	}
  
  
  	if((topLeftValue == topRightValue) &&
	 	(topRightValue == bottomLeftValue) &&
	 	(bottomLeftValue == bottomRightValue))
	{
	  	theValue = topLeftValue;
	  	return true; // Each of four corner values are the same - nothing to interpolate!
	}
  
  
  	if((!sourceGrid->IsMissingValue(topLeftValue)) && (!sourceGrid->IsMissingValue(topRightValue))	&&
	 	(!sourceGrid->IsMissingValue(bottomLeftValue)) && (!sourceGrid->IsMissingValue(bottomRightValue)))
	{
		
		double pres = 0.0;
	  	pres += (1-dx)*(1-dy)*topRightValue;
	  	pres += dx*(1-dy)*bottomRightValue;
		pres += (1-dx)*dy*bottomLeftValue;
		pres += dx*dy*topLeftValue;
		
		theValue = pres;
		
		//newbase
		double interpolatedTopValue = dx*topRightValue + (1. - dx)*topLeftValue;
	  	double interpolatedBottomValue = dx*bottomRightValue + (1. - dx)*bottomLeftValue;
		

		//hil_pp
	  	//double interpolatedTopValue = dx*topLeftValue + (1. - dx)*bottomLeftValue;
	  	//double interpolatedBottomValue = dx*bottomRightValue + (1. - dx)*topRightValue;
	  
	  	theValue = dy*interpolatedTopValue + (1. - dy)*interpolatedBottomValue;
	  	return true;	
	}
  
  	// Missing value(s) found. Use nearest point interpolation method instead
  	sourceGrid->NearestPointInterpolation(x,y,theValue);
  	
  	return true;
}*/


bool compiled_plugin_base::AdjustLeadingDimension(shared_ptr<info> myTargetInfo)
{

    lock_guard<mutex> lock(itsAdjustDimensionMutex);

    // Leading dimension can be: time or level

    if (itsLeadingDimension == kTimeDimension)
    {
        if (!itsFeederInfo->NextTime())
        {
            return false;
        }

        myTargetInfo->Time(itsFeederInfo->Time());
    }
    else if (itsLeadingDimension == kLevelDimension)
    {
        if (!itsFeederInfo->NextLevel())
        {
            return false;
        }

        myTargetInfo->Level(itsFeederInfo->Level());
    }
    else
    {
        throw runtime_error(ClassName() + ": Invalid dimension type: " + boost::lexical_cast<string> (itsLeadingDimension));
    }

    return true;
}

bool compiled_plugin_base::AdjustNonLeadingDimension(shared_ptr<info> myTargetInfo)
{
    if (itsLeadingDimension == kTimeDimension)
    {
        return myTargetInfo->NextLevel();
    }
    else if (itsLeadingDimension == kLevelDimension)
    {
        return myTargetInfo->NextTime();
    }
    else
    {
        throw runtime_error(ClassName() + ": unsupported leading dimension: " + boost::lexical_cast<string> (itsLeadingDimension));
    }
}

void compiled_plugin_base::ResetNonLeadingDimension(shared_ptr<info> myTargetInfo)
{
    if (itsLeadingDimension == kTimeDimension)
    {
        myTargetInfo->ResetLevel();
    }
    else if (itsLeadingDimension == kLevelDimension)
    {
        myTargetInfo->ResetTime();
    }
    else
    {
        throw runtime_error(ClassName() + ": unsupported leading dimension: " + boost::lexical_cast<string> (itsLeadingDimension));
    }
}

himan::level compiled_plugin_base::LevelTransform(const himan::producer& sourceProducer,
													const himan::param& targetParam,
													const himan::level& targetLevel) const
{

	level sourceLevel;

	if (sourceProducer.TableVersion() != kHPMissingInt)
	{
		shared_ptr<neons> n = dynamic_pointer_cast <neons> (plugin_factory::Instance()->Plugin("neons"));

		string lvlName = n->NeonsDB().GetGridLevelName(targetParam.Name(), targetLevel.Type(), 204, sourceProducer.TableVersion());

		HPLevelType lvlType = kUnknownLevel;

		float lvlValue = targetLevel.Value();

		if (lvlName == "GROUND")
		{
			lvlType = kGround;
			lvlValue = 0;
		}
		else if (lvlName == "PRESSURE")
		{
			lvlType = kPressure;
		}
		else if (lvlName == "HYBRID")
		{
			lvlType = kHybrid;
		}
		else if (lvlName == "HEIGHT")
		{
			lvlType = kHeight;
		}
		else
		{
			throw runtime_error(ClassName() + ": Unknown level type: " + lvlName);
		}

		sourceLevel = level(lvlType, lvlValue, lvlName);
	}
	else
	{
		sourceLevel = targetLevel;
	}

	return sourceLevel;
}

bool compiled_plugin_base::SetAB(shared_ptr<info> myTargetInfo, shared_ptr<info> sourceInfo)
{
	if (myTargetInfo->Level().Type() == kHybrid)
	{
		int index = myTargetInfo->ParamIndex();

		myTargetInfo->Grid()->AB(sourceInfo->Grid()->AB());

		myTargetInfo->ParamIndex(index);
	}

	return true;
}

bool compiled_plugin_base::SwapTo(shared_ptr<info> myTargetInfo, HPScanningMode targetScanningMode)
{

	if (myTargetInfo->Grid()->ScanningMode() != targetScanningMode)
	{
		HPScanningMode originalMode = myTargetInfo->Grid()->ScanningMode();

		myTargetInfo->Grid()->ScanningMode(targetScanningMode);

		myTargetInfo->Grid()->Swap(originalMode);
	}

	return true;
}

void compiled_plugin_base::StoreGrib1ParameterDefinitions(vector<param> params, long table2Version)
{
	shared_ptr<neons> n = dynamic_pointer_cast<neons> (plugin_factory::Instance()->Plugin("neons"));

	for (unsigned int i = 0; i < params.size(); i++)
	{
		long parm_id = n->NeonsDB().GetGridParameterId(table2Version, params[i].Name());
		params[i].GribIndicatorOfParameter(parm_id);
		params[i].GribTableVersion(table2Version);
	}
}

void compiled_plugin_base::WriteToFile(shared_ptr<const plugin_configuration> conf, shared_ptr<const info> targetInfo)
{
	shared_ptr<writer> aWriter = dynamic_pointer_cast <writer> (plugin_factory::Instance()->Plugin("writer"));

	// writing might modify iterator positions --> create a copy

	shared_ptr<info> tempInfo(new info(*targetInfo));

	if (conf->FileWriteOption() == kNeons || conf->FileWriteOption() == kMultipleFiles)
	{
		// if info holds multiple parameters, we must loop over them all

		tempInfo->ResetParam();

		while (tempInfo->NextParam())
		{
			aWriter->ToFile(tempInfo, conf);
		}
	}
	else if (conf->FileWriteOption() == kSingleFile)
	{
		aWriter->ToFile(tempInfo, conf, conf->ConfigurationFile());
	}

	tempInfo.reset();
}
