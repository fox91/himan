#include "himan_unit.h"
#include "util.h"

#define BOOST_TEST_MODULE util

using namespace std;
using namespace himan;

const double kEpsilon = 1e-3;

BOOST_AUTO_TEST_CASE(UV_TO_GEOGRAPHICAL)
{
	// Transform grid coordinates to lat and lon in stereographic projection

	himan::point stereoUV(8.484046, 3.804569);
	double lon = 72.79;

	himan::point latlon = util::UVToGeographical(lon, stereoUV);

	BOOST_CHECK_CLOSE(latlon.X(), 6.144442, kEpsilon); 
	BOOST_CHECK_CLOSE(latlon.Y(), -6.978511, kEpsilon); 

	stereoUV.X(-0.2453410);
	stereoUV.Y(0.5808838);
	lon = 23.39;

	latlon = util::UVToGeographical(lon, stereoUV);

	BOOST_CHECK_CLOSE(latlon.X(), 5.4238806e-03, kEpsilon); 
	BOOST_CHECK_CLOSE(latlon.Y(), 0.6305464, kEpsilon); 

}

BOOST_AUTO_TEST_CASE(FILTER2D)
{
	// Filter a plane with given filter kernel
	// Declare matrices
	himan::matrix<double> A(11,8,1);
	himan::matrix<double> B(3,3,1);
	himan::matrix<double> C;
	himan::matrix<double> D(11,8,1);

	// Fill matrix A that will be smoothend with checker-board pattern
	for(size_t i=0; i < A.Size(); ++i)
	{
		if(i % 2 == 0)
      		{
			A.Set(i, 0);
      		}
      		else
      		{
			A.Set(i, 36);
      		}

	}

	// Fill matrix D with solution of the smoothened matrix Filter2D(A,B) 
	for(size_t i=0; i < D.SizeX(); ++i)
	{
		for(size_t j=0; j < D.SizeY(); ++j)
		{
			if(i == 0 || i == A.SizeX()-1 || j == 0 || j == A.SizeY()-1)
			{
				D.Set(i,j,0,18);
			}
			else if ((i % 2 != 0 && j % 2 != 0) || (i % 2 == 0 && j % 2 == 0))
			{
				D.Set(i,j,0,16);
			}
			else
			{
				D.Set(i,j,0,20);
			}
		}
	}

	// Fill matrix B (filter kernel) with constant values 1/9 so that sum(B) == 1
	double filter_coeff(1.0/9.0);
	B.Fill(filter_coeff);

	// Compute smoothened matrix
	C = himan::util::Filter2D(A,B);
	
	// Compare results
	for(size_t i=0; i < C.Size(); ++i)
	{
		BOOST_CHECK_CLOSE(C.At(i),D.At(i),kEpsilon);
	}

	// computed filtered matrix
	std::cout << "Matrix C computed with Filter2D:" << std::endl;
	for (size_t i=0; i < C.SizeX();++i){
    		for (size_t j=0; j < C.SizeY();++j){
      			std::cout << C.At(i,j,0) << " ";
    		}
    		std::cout << std::endl;
  	}

	std::cout << std::endl << "Matrix D as reference case for Filter2D computation:" << std::endl; 

	for (size_t i=0; i < D.SizeX();++i){
    		for (size_t j=0; j < D.SizeY();++j){
      			std::cout << D.At(i,j,0) << " ";
    		}
    		std::cout << std::endl;
  	}

}
