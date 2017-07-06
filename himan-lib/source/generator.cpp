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

template <class InputIt>
std::shared_ptr<himan::info> Max(InputIt begin, InputIt end)
{
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
template std::shared_ptr<himan::info> Max<time_series::iterator>(time_series::iterator, time_series::iterator);
template std::shared_ptr<himan::info> Max<level_series::iterator>(level_series::iterator, level_series::iterator);
template std::shared_ptr<himan::info> Max<std::vector<std::shared_ptr<himan::info>>::iterator>(
    std::vector<std::shared_ptr<himan::info>>::iterator, std::vector<std::shared_ptr<himan::info>>::iterator);

/*
 * Finds the minimum values for a series of info_t on the interval [begin,end)
 * */

template <class InputIt>
std::shared_ptr<himan::info> Min(InputIt begin, InputIt end)
{
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
template std::shared_ptr<himan::info> Min<time_series::iterator>(time_series::iterator, time_series::iterator);
template std::shared_ptr<himan::info> Min<level_series::iterator>(level_series::iterator, level_series::iterator);
template std::shared_ptr<himan::info> Min<std::vector<std::shared_ptr<himan::info>>::iterator>(
    std::vector<std::shared_ptr<himan::info>>::iterator, std::vector<std::shared_ptr<himan::info>>::iterator);

/*
 * Finds value at given height. Interpolates between levels if needed.
 * */

template <class InputIt>
std::shared_ptr<himan::info> Value(InputIt beginValue, InputIt endValue, InputIt beginHeight, InputIt endHeight, double theHeight)
{
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
	std::vector<bool> mask(retInfo->Data().Size());
	std::vector<double> interp(findHeight.size());

	auto prevHeightInfo = *beginHeight;
	auto prevValueInfo = *beginValue;

        for (; beginValue != endValue, beginHeight != endHeight; ++beginValue, ++beginHeight)
        {
                if (*beginValue == nullptr || *beginHeight == nullptr) continue;

		// compute mask for the interpolation step
                auto Height = VEC((*beginHeight)).begin();
                auto PrevHeight = VEC(prevHeightInfo).begin();
		auto FindHeight = findHeight.begin();

		auto Mask = mask.begin();
		auto MaskEnd = mask.end();

		for (; Mask != MaskEnd; ++Mask, ++Height, ++FindHeight, ++PrevHeight)
		{
			*Mask = ((*Height <= *FindHeight) && (*PrevHeight >= *FindHeight)) || ((*Height >= *FindHeight) && (*PrevHeight <= *FindHeight));
		}

		// proceed to next level if no work is to be done for all grid points
		if (!std::accumulate(mask.begin(),mask.end(),false))
			continue;

		// compute interpolated value for all grid points
                auto Value = VEC((*beginValue)).begin();
                auto PrevValue = VEC(prevValueInfo).begin();

                Height = VEC((*beginHeight)).begin();
                PrevHeight = VEC(prevHeightInfo).begin();
                FindHeight = findHeight.begin();

                auto Interp = interp.begin();
		auto InterpEnd = interp.end();


                for (;  Interp != InterpEnd; ++Value, ++Height, ++PrevValue, ++PrevHeight, ++FindHeight, ++Interp)
                {
                        *Interp = (*PrevValue*(*Height-*FindHeight)+*Value*(*FindHeight-*PrevHeight))/(*Height-*PrevHeight);
                }

		// write interpolated value to return vector for grid points with mask value == true 
                Interp = interp.begin();

                Mask = mask.begin();
                MaskEnd = mask.end();

		auto Ret =  VEC((retInfo)).begin();

                for (; Mask != MaskEnd; ++Mask, ++Interp, ++Ret)
                {
                        *Ret = *Mask ? *Interp : *Ret;
                }

		// pass current height and value to previous height and value for next loop iteration
		prevHeightInfo = *beginHeight;
		prevValueInfo = *beginValue;
        }

        return retInfo;
}
template std::shared_ptr<himan::info> Value<level_series::iterator>(level_series::iterator, level_series::iterator, level_series::iterator, level_series::iterator, double);
template std::shared_ptr<himan::info> Value<std::vector<std::shared_ptr<himan::info>>::iterator>(
    std::vector<std::shared_ptr<himan::info>>::iterator, std::vector<std::shared_ptr<himan::info>>::iterator, std::vector<std::shared_ptr<himan::info>>::iterator, std::vector<std::shared_ptr<himan::info>>::iterator, double);

} // close namespace himan
