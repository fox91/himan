#include "generator.h"
#include "ensemble.h"
#include "plugin_factory.h"
#include <numeric>
#include <algorithm>

namespace himan
{
time_series::iterator::iterator(std::shared_ptr<const plugin_configuration> theConfiguration, himan::param theParam,
                                himan::level theLevel, himan::forecast_time theForecastTime,
                                HPTimeResolution theTimeResolution, int theStepSize)
    : itsConfiguration(theConfiguration),
      itsParam(theParam),
      itsLevel(theLevel),
      itsForecastTime(theForecastTime),
      itsTimeResolution(theTimeResolution),
      itsStepSize(theStepSize)
{
	f = GET_PLUGIN(fetcher);
	try
	{
		itsInfo = f->Fetch(itsConfiguration, itsForecastTime, itsLevel, itsParam);
	}
	catch (HPExceptionType& e)
	{
		if (e != kFileDataNotFound)
		{
			abort();
		}
		else
		{
			itsInfo = nullptr;
		}
	}
}

time_series::iterator& time_series::iterator::Next()
{
	itsForecastTime.ValidDateTime().Adjust(itsTimeResolution, itsStepSize);
	try
	{
		itsInfo = f->Fetch(itsConfiguration, itsForecastTime, itsLevel, itsParam);
	}
	catch (HPExceptionType& e)
	{
		if (e != kFileDataNotFound)
		{
			abort();
		}
		else
		{
			itsInfo = nullptr;
		}
	}
	return *this;
}

level_series::iterator::iterator(std::shared_ptr<const plugin_configuration> theConfiguration, himan::param theParam,
                                himan::level theLevel, himan::forecast_time theForecastTime, double theStepSize)
    : itsConfiguration(theConfiguration),
      itsParam(theParam),
      itsLevel(theLevel),
      itsForecastTime(theForecastTime),
      itsStepSize(theStepSize)
{
        f = GET_PLUGIN(fetcher);
        try
        {
                itsInfo = f->Fetch(itsConfiguration, itsForecastTime, itsLevel, itsParam);
        }
        catch (HPExceptionType& e)
        {
                if (e != kFileDataNotFound)
                {
                        abort();
                }
                else
                {
                        itsInfo = nullptr;
                }
        }
}

level_series::iterator& level_series::iterator::Next()
{
        itsLevel.Value(itsLevel.Value()+itsStepSize);
        try
        {
                itsInfo = f->Fetch(itsConfiguration, itsForecastTime, itsLevel, itsParam);
        }
        catch (HPExceptionType& e)
        {
                if (e != kFileDataNotFound)
                {
                        abort();
                }
                else
                {
                        itsInfo = nullptr;
                }
        }
        return *this;
}

/*
 * Finds the maximum values for a series of info_t on the interval [begin,end)
 * */

template <class RngExpr>
himan::matrix<double> Max(RngExpr value)
{
	auto begin = value.begin();
	auto end = value.end();

	// Set first field as first set of maximum values
	himan::matrix<double> ret = *begin;
	++begin;

	// Update set of maximum values
	for (; begin != end; ++begin)
	{

		auto in = (*begin).Values().begin();
		auto out = ret.Values().begin();
		auto inEnd = (*begin).Values().end();
		auto outEnd = ret.Values().end();

		for (; in != inEnd, out != outEnd; ++in, ++out)
		{
			*out = std::max(*in, *out);
		}
	}

	return ret;
}
template himan::matrix<double> Max<time_series>(time_series);
template himan::matrix<double> Max<level_series>(level_series);
template himan::matrix<double> Max<std::vector<himan::matrix<double>>>(std::vector<himan::matrix<double>>);

/*
 * Finds the minimum values for a series of info_t on the interval [begin,end)
 * */

template <class RngExpr>
himan::matrix<double>  Min(RngExpr value)
{
	auto begin = value.begin();
	auto end = value.end();

	himan::matrix<double> ret = *begin;	
	++begin;

	for (; begin != end; ++begin)
	{
		auto in = (*begin).Values().begin();
		auto out = ret.Values().begin();
		auto inEnd = (*begin).Values().end();
		auto outEnd = ret.Values().end();

		for (; in != inEnd, out != outEnd; ++in, ++out)
		{
			*out = std::min(*in, *out);
		}
	}

	return ret;
}
template himan::matrix<double> Min<time_series>(time_series);
template himan::matrix<double> Min<level_series>(level_series);

/*
 * Finds maximum value within given vertical gange.
 * */

template <class RngExpr>
himan::matrix<double> Max(RngExpr value, RngExpr height, double theLowerHeight, double theUpperHeight)
{
        auto beginValue = value.begin();
        auto endValue = value.end();
        auto beginHeight = height.begin();
        auto endHeight = height.end();

        auto ret = *beginValue;

        std::vector<double> lowerHeight(ret.Size(),theLowerHeight);
        std::vector<double> upperHeight(ret.Size(),theUpperHeight);

        auto prevHeightInfo = *beginHeight;
        auto prevValueInfo = *beginValue;

        for (; beginValue != endValue, beginHeight != endHeight; ++beginValue, ++beginHeight)
        {
                himan::matrix<double> heightInfo = *beginHeight;
                himan::matrix<double> valueInfo = *beginValue;

                auto Ret =  ret.Values().begin();
                auto RetEnd =  ret.Values().end();

                auto Value = valueInfo.Values().begin();
                auto PrevValue = prevValueInfo.Values().begin();

                auto Height = heightInfo.Values().begin();
                auto PrevHeight = prevHeightInfo.Values().begin();
                auto LowerHeight = lowerHeight.begin();
                auto UpperHeight = upperHeight.begin();

                for (; Ret != RetEnd; ++Value, ++Height, ++PrevValue, ++PrevHeight, ++LowerHeight, ++UpperHeight, ++Ret)
                {
                        register bool Mask = ((*Height <= *LowerHeight) && (*Height >= *UpperHeight)) || ((*Height >= *LowerHeight) && (*Height <= *UpperHeight));

                        *Ret = Mask ? std::max(*Value,*Ret) : *Ret;
                }

                // pass current height and value to previous height and value for next loop iteration
                prevHeightInfo = heightInfo;
                prevValueInfo = valueInfo;
        }

        return ret;
}
template himan::matrix<double> Max<level_series>(level_series, level_series, double, double);

/*
 * Finds value at given vertical coordinate. Interpolates between levels if needed.
 * */

template <class RngExpr>
himan::matrix<double> Value(RngExpr value, RngExpr height, double theHeight)
{

	auto beginValue = value.begin();
	auto endValue = value.end();
	auto beginHeight = height.begin();
	auto endHeight = height.end();


        himan::matrix<double> ret = *beginValue;

	std::vector<double> findHeight(ret.Size(),theHeight);

	himan::matrix<double> prevHeightInfo = *beginHeight;
	himan::matrix<double> prevValueInfo = *beginValue;

        for (; beginValue != endValue, beginHeight != endHeight; ++beginValue, ++beginHeight)
        {
		himan::matrix<double> heightInfo = *beginHeight;
                himan::matrix<double> valueInfo = *beginValue;

		auto Ret = ret.Values().begin();
                auto RetEnd = ret.Values().end();

                auto Value = valueInfo.Values().begin();
                auto PrevValue = prevValueInfo.Values().begin();

                auto Height = heightInfo.Values().begin();
                auto PrevHeight = prevHeightInfo.Values().begin();
                auto FindHeight = findHeight.begin();


                for (; Ret != RetEnd; ++Value, ++Height, ++PrevValue, ++PrevHeight, ++FindHeight, ++Ret)
                {
                        register bool Mask = ((*Height <= *FindHeight) && (*PrevHeight >= *FindHeight)) || ((*Height >= *FindHeight) && (*PrevHeight <= *FindHeight));
                        register double Interp = (*PrevValue*(*Height-*FindHeight)+*Value*(*FindHeight-*PrevHeight))/(*Height-*PrevHeight);

                        *Ret = Mask ? Interp : *Ret;
                }

		// pass current height and value to previous height and value for next loop iteration
		prevHeightInfo = heightInfo;
		prevValueInfo = valueInfo;
        }

        return ret;
}
template himan::matrix<double> Value<level_series>(level_series, level_series, double);

} // close namespace himan
