#include "generator.h"
#include "plugin_factory.h"
#include <numeric>

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
std::shared_ptr<himan::info> Max(RngExpr value)
{
	auto begin = value.begin();
	auto end = value.end();

	if (begin == end) return nullptr;

	// Find first field that contains data
	while (*begin == nullptr)
	{
		++begin;
		if (begin == end) return nullptr;
	}

	// Set first field as first set of maximum values
	auto maxInfo = *begin;
	maxInfo->ReGrid();
	++begin;

	// Update set of maximum values
	for (; begin != end; ++begin)
	{
		if (*begin == nullptr) continue;

		auto in = VEC((*begin)).begin();
		auto out = VEC(maxInfo).begin();
		auto inEnd = VEC((*begin)).end();
		auto outEnd = VEC(maxInfo).end();

		for (; in != inEnd, out != outEnd; ++in, ++out)
		{
			*out = std::max(*in, *out);
		}
	}

	return maxInfo;
}
template std::shared_ptr<himan::info> Max<time_series>(time_series);
template std::shared_ptr<himan::info> Max<level_series>(level_series);
template std::shared_ptr<himan::info> Max<std::vector<std::shared_ptr<himan::info>>>(std::vector<std::shared_ptr<himan::info>>);

/*
 * Finds the minimum values for a series of info_t on the interval [begin,end)
 * */

template <class RngExpr>
std::shared_ptr<himan::info> Min(RngExpr value)
{
	auto begin = value.begin();
	auto end = value.end();

	if (begin == end) return nullptr;

	while (*begin == nullptr)
	{
		++begin;
		if (begin == end) return nullptr;
	}

	auto minInfo = *begin;
	minInfo->ReGrid();
	++begin;

	for (; begin != end; ++begin)
	{
		if (*begin == nullptr) continue;

		auto in = VEC((*begin)).begin();
		auto out = VEC(minInfo).begin();
		auto inEnd = VEC((*begin)).end();
		auto outEnd = VEC(minInfo).end();

		for (; in != inEnd, out != outEnd; ++in, ++out)
		{
			*out = std::min(*in, *out);
		}
	}

	return minInfo;
}
template std::shared_ptr<himan::info> Min<time_series>(time_series);
template std::shared_ptr<himan::info> Min<level_series>(level_series);
template std::shared_ptr<himan::info> Min<std::vector<std::shared_ptr<himan::info>>>(std::vector<std::shared_ptr<himan::info>>);

/*
 * Finds value at given height. Interpolates between levels if needed.
 * */

template <class RngExpr>
std::shared_ptr<himan::info> Value(RngExpr value, RngExpr height, double theHeight)
{

	auto beginValue = value.begin();
	auto endValue = value.end();
	auto beginHeight = height.begin();
	auto endHeight = height.end();

        if (beginValue == endValue || beginHeight == endHeight) return nullptr;

        while (*beginValue == nullptr || *beginHeight == nullptr )
        {
                ++beginValue;
		++beginHeight;
                if (beginValue == endValue || beginHeight == endHeight) return nullptr;
        }

        auto retInfo = *beginValue;
        retInfo->ReGrid();
	retInfo->Data().Fill(kFloatMissing);

	std::vector<double> findHeight(retInfo->Data().Size(),theHeight);

	auto prevHeightInfo = *beginHeight;
	auto prevValueInfo = *beginValue;

        for (; beginValue != endValue, beginHeight != endHeight; ++beginValue, ++beginHeight)
        {
                if (*beginValue == nullptr || *beginHeight == nullptr) continue;

		auto Ret =  VEC((retInfo)).begin();
                auto RetEnd =  VEC((retInfo)).end();

                auto Value = VEC((*beginValue)).begin();
                auto PrevValue = VEC(prevValueInfo).begin();

                auto Height = VEC((*beginHeight)).begin();
                auto PrevHeight = VEC(prevHeightInfo).begin();
                auto FindHeight = findHeight.begin();


                for (; Ret != RetEnd; ++Value, ++Height, ++PrevValue, ++PrevHeight, ++FindHeight, ++Ret)
                {
                        register bool Mask = ((*Height <= *FindHeight) && (*PrevHeight >= *FindHeight)) || ((*Height >= *FindHeight) && (*PrevHeight <= *FindHeight));
                        register double Interp = (*PrevValue*(*Height-*FindHeight)+*Value*(*FindHeight-*PrevHeight))/(*Height-*PrevHeight);

                        *Ret = Mask ? Interp : *Ret;
                }

		// pass current height and value to previous height and value for next loop iteration
		prevHeightInfo = *beginHeight;
		prevValueInfo = *beginValue;
        }

        return retInfo;
}
template std::shared_ptr<himan::info> Value<level_series>(level_series, level_series, double);
template std::shared_ptr<himan::info> Value<std::vector<std::shared_ptr<himan::info>>>(std::vector<std::shared_ptr<himan::info>>, std::vector<std::shared_ptr<himan::info>>, double);

} // close namespace himan
