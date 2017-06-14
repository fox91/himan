#ifndef GENERATOR_H
#define GENERATOR_H

#include "fetcher.h"
#include "himan_common.h"
#include <iterator>

namespace himan
{
/*
 *  This is a generator type object that allows to fetch a time series of data on the fly.
 *  It is build upon the input generator model. Incrementing the iterator fetches the next
 *  time step of the desired parameter.
 *
 *  TODO A range check for the begin and end time, i.e. begin-end = n*time resolution.
 * */

class time_series
{
   public:
	class iterator : public std::iterator<std::input_iterator_tag, std::shared_ptr<himan::info>>
	{
	   public:
		explicit iterator(std::shared_ptr<const plugin_configuration> theConfiguration, himan::param theParam,
		                  himan::level theLevel, himan::forecast_time theForecastTime,
		                  HPTimeResolution theTimeResolution, int theStepSize);
		iterator& operator++() { return Next(); }
		iterator operator++(int)
		{
			iterator retval = *this;
			++(*this);
			return retval;
		}
		bool operator==(iterator other) const
		{
			return (itsForecastTime == other.itsForecastTime && itsLevel == other.itsLevel);
		}
		bool operator!=(iterator other) const { return !(*this == other); }
		const std::shared_ptr<himan::info>& operator*() const { return itsInfo; }
	   private:
		iterator& Next();
		std::shared_ptr<himan::plugin::fetcher> f;
		std::shared_ptr<const plugin_configuration> itsConfiguration;
		himan::param itsParam;
		himan::level itsLevel;
		himan::forecast_time itsForecastTime;
		HPTimeResolution itsTimeResolution;
		int itsStepSize;

		std::shared_ptr<himan::info> itsInfo;
	};

	time_series(std::shared_ptr<const plugin_configuration> theConfiguration, himan::param theParam,
	            himan::level theLevel, himan::forecast_time theStartTime, himan::forecast_time theEndTime,
	            HPTimeResolution theTimeResolution, int theStepSize)
	    : itsBegin(theConfiguration, theParam, theLevel, theStartTime, theTimeResolution, theStepSize),
	      itsEnd(theConfiguration, theParam, theLevel, theEndTime, theTimeResolution, theStepSize)
	{
	}

	iterator begin() const { return itsBegin; }
	iterator end() const { return itsEnd; }
   private:
	iterator itsBegin;
	iterator itsEnd;
};

// TODO these functions are independent of the implementation of the generator class but are example use cases.
// They should be moved to some other namespace, e.g. numerical functions, util, modifier, etc.
template <class InputIt>
std::shared_ptr<himan::info> Max(InputIt, InputIt);

template <class InputIt>
std::shared_ptr<himan::info> Min(InputIt, InputIt);

}  // namespace himan

#endif  // GENERATOR_H
