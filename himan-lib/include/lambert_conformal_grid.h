/**
 * @file   lambert_conformal_grid.h
 *
 */

#ifndef LAMBERT_CONFORMAL_GRID_H
#define LAMBERT_CONFORMAL_GRID_H

#include "grid.h"
#include "logger.h"
#include "point.h"
#include "serialization.h"
#include <string>

class OGRCoordinateTransformation;

namespace himan
{
class lambert_conformal_grid : public regular_grid
{
   public:
	lambert_conformal_grid(HPScanningMode theScanningMode, const point& theFirstPoint, size_t ni, size_t nj, double di,
	                       double dj, double theOrientation, double theStandardParallel1, double theStandardParallel2,
	                       const earth_shape<double>& earthShape, bool firstPointIsProjected = false);
	lambert_conformal_grid(HPScanningMode theScanningMode, const point& theFirstPoint, size_t ni, size_t nj, double di,
	                       double dj, std::unique_ptr<OGRSpatialReference> spRef, bool firstPointIsProjected = false);

	virtual ~lambert_conformal_grid() = default;
	lambert_conformal_grid(const lambert_conformal_grid& other);
	lambert_conformal_grid& operator=(const lambert_conformal_grid& other) = delete;

	virtual std::string ClassName() const
	{
		return "himan::lambert_conformal_grid";
	}
	virtual std::ostream& Write(std::ostream& file) const;

	bool operator==(const grid& other) const;
	bool operator!=(const grid& other) const;

	std::unique_ptr<grid> Clone() const override;

	double Orientation() const;
	double StandardParallel1() const;
	double StandardParallel2() const;

	OGRSpatialReference SpatialReference() const;

	size_t Hash() const override;
	double Cone() const;

   private:
	bool EqualsTo(const lambert_conformal_grid& other) const;
	void CreateCoordinateTransformations(const point& firstPoint, bool isProjected);

#ifdef SERIALIZATION
	friend class cereal::access;

	template <class Archive>
	void serialize(Archive& ar)
	{
		ar(cereal::base_class<regular_grid>(this));
	}
#endif
};

inline std::ostream& operator<<(std::ostream& file, const lambert_conformal_grid& ob)
{
	return ob.Write(file);
}
}  // namespace himan

#ifdef SERIALIZATION
CEREAL_REGISTER_TYPE(himan::lambert_conformal_grid);
#endif

#endif /* LAMBERT_CONFORMAL_GRID_H */
