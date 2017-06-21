#include "generator.h"
#include "plugin_factory.h"

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

	for (; begin != end; ++begin)
	{
		if (*begin == nullptr) continue;

		auto in = VEC((*begin)).begin();
		auto out = VEC(minInfo).begin();
		auto inEnd = VEC((*begin)).end();
		auto outEnd = VEC(minInfo).end();

		for (; out != inEnd, out != outEnd; ++in, ++out)
		{
			*out = std::min(*in, *out);
		}
	}

	return minInfo;
}
template std::shared_ptr<himan::info> Min<time_series::iterator>(time_series::iterator, time_series::iterator);
template std::shared_ptr<himan::info> Min<std::vector<std::shared_ptr<himan::info>>::iterator>(
    std::vector<std::shared_ptr<himan::info>>::iterator, std::vector<std::shared_ptr<himan::info>>::iterator);
}
