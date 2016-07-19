/**
 * @file reduced_gaussian_grid.cpp
 *
 *  Created on: Jun 18, 2016
 *  @author partio
 */

#include "reduced_gaussian_grid.h"
#include "logger_factory.h"
#include <boost/lexical_cast.hpp>

using namespace himan;

reduced_gaussian_grid::reduced_gaussian_grid()
	: grid(kIrregularGrid, kReducedGaussian)
	, itsN(kHPMissingInt)
	, itsNumberOfLongitudesAlongParallels()
	, itsNj(kHPMissingInt)
	, itsBottomLeft()
	, itsTopRight()
	, itsBottomRight()
	, itsTopLeft()
	, itsDj(kHPMissingValue)
{
	itsLogger = logger_factory::Instance()->GetLog("reduced_gaussian_grid");
}

reduced_gaussian_grid::reduced_gaussian_grid(const reduced_gaussian_grid& other) 
	: grid(other)
	, itsN(other.itsN)
	, itsNumberOfLongitudesAlongParallels(other.itsNumberOfLongitudesAlongParallels)
	, itsNj(other.itsNj)
	, itsBottomLeft(other.itsBottomLeft)
	, itsTopRight(other.itsTopRight)
	, itsBottomRight(other.itsBottomRight)
	, itsTopLeft(other.itsTopLeft)
	, itsDj(kHPMissingValue)
{
	itsLogger = logger_factory::Instance()->GetLog("reduced_gaussian_grid");
}

int reduced_gaussian_grid::N() const
{
	return itsN;
}

void reduced_gaussian_grid::N(int theN)
{
	itsN = theN;
}

size_t reduced_gaussian_grid::Size() const
{
	assert(itsNumberOfLongitudesAlongParallels.size() == 2 * static_cast<size_t> (itsN));

	return std::accumulate(itsNumberOfLongitudesAlongParallels.begin(), itsNumberOfLongitudesAlongParallels.end(), 0);
}

std::vector<int> reduced_gaussian_grid::NumberOfLongitudesAlongParallels() const
{
	return itsNumberOfLongitudesAlongParallels;
}

void reduced_gaussian_grid::NumberOfLongitudesAlongParallels(std::vector<int> theNumberOfLongitudesAlongParallels)
{
	assert((itsN == kHPMissingInt && itsNumberOfLongitudesAlongParallels.size() == 0) || static_cast<size_t> (itsN*2) == theNumberOfLongitudesAlongParallels.size());
	itsNumberOfLongitudesAlongParallels = theNumberOfLongitudesAlongParallels;
}

point reduced_gaussian_grid::FirstPoint() const
{
	switch (itsScanningMode)
	{
		case kTopLeft:
			return itsTopLeft;
		case kBottomLeft:
			return itsBottomLeft;
		default:
			throw std::runtime_error("Scanning mode not supported: " + HPScanningModeToString.at(itsScanningMode));
	}
}

point reduced_gaussian_grid::LastPoint() const
{
	switch (itsScanningMode)
	{
		case kTopLeft:
			return itsBottomRight;
		case kBottomLeft:
			return itsTopRight;
		default:
			throw std::runtime_error("Scanning mode not supported: " + HPScanningModeToString.at(itsScanningMode));
	}
}

std::ostream& reduced_gaussian_grid::Write(std::ostream& file) const
{
	grid::Write(file);
	
	file << "__itsN__ " << itsN << std::endl;
	file << "__itsNj__ " << itsNj << std::endl;
	file << "__itsNumberOfLongitudesAlongParallels__";
	
	for (auto& num : itsNumberOfLongitudesAlongParallels)
	{
		file << " " << num;
	}

	file << std::endl;
	
	file << itsBottomLeft;
	file << itsTopRight;
	file << itsBottomRight;
	file << itsTopLeft;

	return file;
}

point reduced_gaussian_grid::LatLon(size_t x, size_t y) const
{
	assert(itsScanningMode == kTopLeft);
	assert(TopLeft() != point());
	assert(BottomRight() != point());

	double lonspan = (BottomRight().X() - TopLeft().X()); // longitude span of the whole area in degrees
	lonspan = (lonspan < 0) ? lonspan + 360 : lonspan;
	assert(lonspan >= 0 && lonspan <= 360);

	const size_t currentNumOfLongitudes = itsNumberOfLongitudesAlongParallels[y];
	const double di = (lonspan / (static_cast<double> (currentNumOfLongitudes)-1.));
	const double dj = Dj();

	double latitude = TopLeft().Y() - static_cast<double> (y) * dj;
	double longitude = TopLeft().X() + static_cast<double> (x) * di;

	if (longitude > 360)
	{
		longitude -= 360;
	}

	return point(longitude, latitude);
}

point reduced_gaussian_grid::LatLon(size_t theLocationIndex) const
{
	assert(itsNj > 0);
	assert(itsNumberOfLongitudesAlongParallels.size() > 0);
	
	if (theLocationIndex < static_cast<size_t> (itsNumberOfLongitudesAlongParallels[0]))
	{
		return LatLon(theLocationIndex, 0);
	}
	
	size_t sum = 0, y = 0;

	for (size_t i = 0; i < itsNj; i++)
	{
		int numLongitudes = itsNumberOfLongitudesAlongParallels[i];
		
		if (sum + numLongitudes > theLocationIndex)
		{
			y = i;
			break;
		}
		
		sum += numLongitudes;
	}

	assert(theLocationIndex >= sum);
	
	size_t x = theLocationIndex - sum;
	assert(y > 0);

	return LatLon(x, y);
}

point reduced_gaussian_grid::LatLonToGridPoint(const himan::point& latlon) const
{
	double offset = 0;

	if (TopLeft().X() == 0 && (BottomRight().X() < 0 || BottomRight().X() > 180))
	{
		offset = 360;
	}

	const double dj = Dj();
	assert(dj > 0.);

	const double y = static_cast<int> ((TopLeft().Y() - latlon.Y()) / dj); // grid y [0 .. Nj-1] 

	assert(itsNj > 0);

	if (y < 0. || y > (static_cast<double> (itsNj)-1.))
	{
		// out of bounds
		return point();
	}

	const int numCurrentLongitudes = itsNumberOfLongitudesAlongParallels[static_cast<size_t> (rint(y))]; // number of longitudes for the nearest parallel

	assert(numCurrentLongitudes > 0);

	double lonspan = (BottomRight().X() - TopLeft().X()); // longitude span of the whole area in degrees
	lonspan = (lonspan < 0) ? lonspan + 360 : lonspan;
	assert(lonspan >= 0 && lonspan <= 360);
	
	const double di = (lonspan / (numCurrentLongitudes-1)); // longitude distance between two points in degrees for the current parallel
	double x = (latlon.X() - (TopLeft().X() - offset)) / di; //grid x in current parallel

	if (offset != 0)
	{
		x = fmod(x, numCurrentLongitudes); // wrap around with global fields
	}
	
	return point(x, y);
}

double reduced_gaussian_grid::Value(size_t x, size_t y) const
{
	if (y == 0)
	{
		return itsData.At(x);
	}
	
	size_t sum = 0, i = 0;
	
	while (i < itsNj && i <= (y-1))
	{
		sum += itsNumberOfLongitudesAlongParallels[i];
		i++;
	}

	return itsData.At(sum + x);
}

bool reduced_gaussian_grid::Swap(HPScanningMode newScanningMode)
{
	// not implemented yet
	return false;
}

reduced_gaussian_grid* reduced_gaussian_grid::Clone() const
{
	return new reduced_gaussian_grid(*this);
}

double reduced_gaussian_grid::Di() const
{
	return kHPMissingValue;
}

double reduced_gaussian_grid::Dj() const
{
	if (itsDj == kHPMissingValue)
	{	
		assert(itsNj != kHPMissingInt);
		itsDj = (TopLeft().Y() - BottomRight().Y()) / (static_cast<double> (itsNj)-1.);
	}

	return itsDj;
}

size_t reduced_gaussian_grid::Ni() const
{
	// no bound checking ...
	return itsNumberOfLongitudesAlongParallels[itsN];
}

void reduced_gaussian_grid::Nj(size_t theNj)
{
	itsNj = theNj;
}

size_t reduced_gaussian_grid::Nj() const
{
	return itsNj;
}

point reduced_gaussian_grid::BottomLeft() const
{
	return itsBottomLeft;
}

point reduced_gaussian_grid::TopRight() const
{
	return itsTopRight;
}

void reduced_gaussian_grid::BottomLeft(const point& theBottomLeft)
{
	itsBottomLeft = theBottomLeft;
	UpdateCoordinates();
}

void reduced_gaussian_grid::TopRight(const point& theTopRight)
{
	itsTopRight = theTopRight;
	UpdateCoordinates();
}

point reduced_gaussian_grid::BottomRight() const
{
	return itsBottomRight;
}

point reduced_gaussian_grid::TopLeft() const
{
	return itsTopLeft;
}

void reduced_gaussian_grid::BottomRight(const point& theBottomRight)
{
	itsBottomRight = theBottomRight;
	UpdateCoordinates();
}

void reduced_gaussian_grid::TopLeft(const point& theTopLeft)
{
	itsTopLeft = theTopLeft;
	UpdateCoordinates();
}

void reduced_gaussian_grid::UpdateCoordinates() const
{
	const point missing;
	
	if ((itsBottomLeft != missing && itsTopRight != missing) && (itsTopLeft == missing || itsBottomRight == missing))
	{
		itsTopLeft = point(itsBottomLeft.X(), itsTopRight.Y());
		itsBottomRight = point(itsTopRight.X(), itsBottomLeft.Y());
	}
	else if ((itsBottomLeft == missing || itsTopRight == missing) && (itsTopLeft != missing && itsBottomRight != missing))
	{
		itsBottomLeft = point(itsTopLeft.X(), itsBottomRight.Y());
		itsTopRight = point(itsBottomRight.X(), itsTopLeft.Y());	
	}
}

bool reduced_gaussian_grid::operator!=(const grid& other) const
{
	return !(other == *this);
}

bool reduced_gaussian_grid::operator==(const grid& other) const
{
	const reduced_gaussian_grid* g = dynamic_cast<const reduced_gaussian_grid*> (&other);

	if (g)
	{
		return EqualsTo(*g);
	}

	return false;
}

bool reduced_gaussian_grid::EqualsTo(const reduced_gaussian_grid& other) const
{
	if (!grid::EqualsTo(other))
	{
		return false;
	}

	if (itsBottomLeft != other.BottomLeft())
	{
		itsLogger->Trace("BottomLeft does not match: X " + boost::lexical_cast<std::string> (itsBottomLeft.X()) + " vs " + boost::lexical_cast<std::string> (other.BottomLeft().X()));
		itsLogger->Trace("BottomLeft does not match: Y " + boost::lexical_cast<std::string> (itsBottomLeft.Y()) + " vs " + boost::lexical_cast<std::string> (other.BottomLeft().Y()));
		return false;
	}

	if (itsTopRight != other.TopRight())
	{
		itsLogger->Trace("TopRight does not match: X " + boost::lexical_cast<std::string> (itsTopRight.X()) + " vs " + boost::lexical_cast<std::string> (other.TopRight().X()));
		itsLogger->Trace("TopRight does not match: Y " + boost::lexical_cast<std::string> (itsTopRight.Y()) + " vs " + boost::lexical_cast<std::string> (other.TopRight().Y()));
		return false;
	}

	if (itsN != other.N())
	{
		itsLogger->Trace("N does not match: " + boost::lexical_cast<std::string> (itsN) + " vs " + boost::lexical_cast<std::string> (other.N()));
		return false;
	}

	if (itsNj != other.Nj())
	{
		itsLogger->Trace("Nj does not match: " + boost::lexical_cast<std::string> (itsNj) + " vs " + boost::lexical_cast<std::string> (other.Nj()));
		return false;
	}

	if (itsNumberOfLongitudesAlongParallels.size() != other.NumberOfLongitudesAlongParallels().size())
	{
		itsLogger->Trace("NumberOfLongitudesAlongParallels size does not match: " + boost::lexical_cast<std::string> (itsNumberOfLongitudesAlongParallels.size()) 
			+ " vs " + boost::lexical_cast<std::string> (other.NumberOfLongitudesAlongParallels().size()));
		return false;
	}
	else
	{
		for (size_t i = 0; i < itsNumberOfLongitudesAlongParallels.size(); i++)
		{
			if (itsNumberOfLongitudesAlongParallels[i] != other.NumberOfLongitudesAlongParallels()[i])
			{
				return false;
			}
		}
	}
	
	return true;
}