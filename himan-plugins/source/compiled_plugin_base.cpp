/**
 *
 * @file compiled_plugin_base.cpp
 *
 * @date Jan 15, 2013
 * @author partio
 */

#include "compiled_plugin_base.h"
#include <boost/thread.hpp>

using namespace std;
using namespace himan::plugin;

const unsigned int MAX_THREADS = 12; //<! Max number of threads we allow
const double kInterpolatedValueEpsilon = 0.00001; //<! Max difference between two grid points (if smaller, points are considered the same)

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
	const NFmiPoint targetLatLonPoint = targetGrid->LatLon();
	const NFmiPoint targetGridPoint = targetGrid->GridPoint();

	if (gridsAreEqual)
	{
		value = sourceGrid->FloatValue(targetGridPoint);
		return true;
	}


	const NFmiPoint sourceGridPoint = sourceGrid->LatLonToGrid(targetLatLonPoint);

	bool noInterpolation = (fabs(targetGridPoint.X() - round(sourceGridPoint.X())) < kInterpolatedValueEpsilon &&
		 fabs(targetGridPoint.Y() - round(sourceGridPoint.Y())) < kInterpolatedValueEpsilon);

	if (noInterpolation)
	{
		value = sourceGrid->FloatValue();
		return true;
	}

	return sourceGrid->InterpolateToLatLonPoint(targetLatLonPoint, value);

}

bool compiled_plugin_base::AdjustLeadingDimension(shared_ptr<info> myTargetInfo)
{

    //lock_guard<mutex> lock(itsAdjustDimensionMutex);

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
