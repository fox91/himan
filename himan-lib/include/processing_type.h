/**
 * @file processing_type.h
 *
 * @brief simple class to describe parameter processing type metadata,
 * like for example probability related metadata
 */

#pragma once

#include "himan_common.h"
#include "serialization.h"

namespace himan
{
enum HPProcessingType
{
	kUnknownProcessingType = 0,
	kProbabilityGreaterThan,
	kProbabilityLessThan,
	kProbabilityBetween,
	kProbabilityEquals,
	kProbabilityNotEquals,
	kProbabilityEqualsIn,
	kFractile,
	kEnsembleMean,  // to differ from aggregation 'average'
	kSpread,
	kStandardDeviation,
	kEFI
};

const boost::unordered_map<HPProcessingType, std::string> HPProcessingTypeToString =
    ba::map_list_of(kUnknownProcessingType, "unknown")(kProbabilityGreaterThan, "probability greater than")(
        kProbabilityLessThan, "probability less than")(kProbabilityBetween, "probability between")(
        kProbabilityEquals, "probability equals")(kProbabilityNotEquals, "probability is not equal to")(
        kProbabilityEqualsIn, "probability equal to value in list")(kFractile, "fractile")(
        kEnsembleMean, "ensemble mean")(kSpread, "spread")(kStandardDeviation, "standard deviation")(kEFI, "efi");

class processing_type
{
   public:
	processing_type() = default;
	processing_type(HPProcessingType theType);
	processing_type(HPProcessingType theType, double theValue, double theValue2);
	processing_type(HPProcessingType theType, double theValue, double theValue2, int theNumberOfEnsembleMembers);
	~processing_type() = default;

	processing_type(const processing_type& other) = default;
	processing_type& operator=(const processing_type& other) = default;

	bool operator==(const processing_type& other) const;
	bool operator!=(const processing_type& other) const;
	operator std::string() const;

	std::string ClassName() const
	{
		return "himan::processing_type";
	}

	HPProcessingType Type() const;
	void Type(HPProcessingType theType);

	double Value() const;
	void Value(double theValue);

	double Value2() const;
	void Value2(double theValue2);

	/**
	 * @brief Return how many ensemble members were used in the processing (if any)
	 */

	int NumberOfEnsembleMembers() const;
	void NumberOfEnsembleMembers(int theNumberOfEnsembleMembers);

	std::ostream& Write(std::ostream& file) const;

   private:
	HPProcessingType itsType = kUnknownProcessingType;
	double itsValue = kHPMissingValue;
	double itsValue2 = kHPMissingValue;

	int itsNumberOfEnsembleMembers = kHPMissingInt;

#ifdef SERIALIZATION
	friend class cereal::access;

	template <class Archive>
	void serialize(Archive& ar)
	{
		ar(CEREAL_NVP(itsType), CEREAL_NVP(itsFirstValue), CEREAL_NVP(itsSecondValue),
		   CEREAL_NVP(itsNumberOfEnsembleMembers));
	}
#endif
};

inline std::ostream& operator<<(std::ostream& file, const processing_type& ob)
{
	return ob.Write(file);
}
}  // namespace himan
