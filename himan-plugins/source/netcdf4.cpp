/**
 * @file netcdf4.cpp
 *
 */

#include "fminc4.h"
#include "variable.h"
#include "group.h"
#include "dimension.h"
#include "common.h"

#include "grid.h"

#include "producer.h"
#include "latitude_longitude_grid.h"
#include "logger.h"
#include <fstream>

#include "netcdf4.h"

#include "plugin_factory.h"
#include "radon.h"

#ifdef __clang__

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wsign-conversion"
#pragma clang diagnostic ignored "-Winvalid-source-encoding"
#pragma clang diagnostic ignored "-Wshorten-64-to-32"

#endif

#ifdef __clang__

#pragma clang diagnostic pop

#endif

using namespace std;
using namespace himan;
using namespace himan::plugin;
using namespace fminc4;

//--------------------------------------------------------------------------------
// TODO Maybe should be member functions of grid_latitude_longitude
vector<double> ComputeLatsFromLatLonGrid(const shared_ptr<himan::grid>& grid)
{
	auto rg = dynamic_pointer_cast<latitude_longitude_grid>(grid);

        himan::point firstGridPoint = rg->FirstPoint();
        himan::point lastGridPoint = rg->LastPoint();

	switch (rg->ScanningMode())
        {
        	case kTopLeft:
                {
	                std::vector<double> w(rg->Nj());
                        for(size_t j = 0; j < rg->Nj(); ++j)
                        {
        	                w[j] = firstGridPoint.Y() - double(j) * rg->Dj();
                        }
			return w;
                }
		default:
			cout << "not implemented\n";
	}
	return vector<double>{};

}

vector<double> ComputeLonsFromLatLonGrid(const shared_ptr<himan::grid>& grid)
{
	auto rg = dynamic_pointer_cast<latitude_longitude_grid>(grid);

        himan::point firstGridPoint = rg->FirstPoint();
        himan::point lastGridPoint = rg->LastPoint();

	switch (rg->ScanningMode())
        {
                case kTopLeft:
                {       
                        std::vector<double> v(rg->Ni());
                        for(size_t i = 0; i < rg->Ni(); ++i)
                        {       
                                v[i] = firstGridPoint.X() + double(i) * rg->Di();
                        }
			return v;
		}
		default:
			cout << "not implemented\n";
	}
	return vector<double>{};
}
//--------------------------------------------------------------------------------

bool CreateMetadataTime(nc_group& file, const forecast_time& ftime, const producer& p)
{
	// create netcdf dimension and according variable
	try
	{
		nc_dim d = file.AddDim("time", NC_UNLIMITED);
		file.AddVar<long>("time",{d});
	}
	catch(const int& err)
	{
		cout << nc_strerror(err) << '\n';
		return false;
	}

        auto r = GET_PLUGIN(radon);
        r->RadonDB().Query("SELECT units FROM cf_time where producer_id = " + to_string(p.Id()));
        const auto row = r->RadonDB().FetchRow();

        string tMask = row[0];
        
	try
	{
        	nc_var<long> t = file.GetVar<long>("time");
		t.AddTextAtt("long_name","time");
        	t.AddTextAtt("units",ftime.OriginDateTime().String(tMask));
	}
	catch(const int& err)
        {
                cout << nc_strerror(err) << '\n';
                return false;
        }

	return true;
}

bool CreateMetadataGrid(nc_group& file, const shared_ptr<himan::grid>& grid, const producer& p)
{
	switch (grid->Type())
        {
                case kLatitudeLongitude:
                {
                        auto rg = dynamic_pointer_cast<latitude_longitude_grid>(grid);

			try
			{
                        	//auto lon = file.AddDim("longitude", rg->Ni());
				//auto longvar = file.AddVar<double>("longitude",{lon});
			}
			catch(const int& err)
			{
				return false;
			}

			/*try
			{
				auto longvar = file.GetVar<double>("longitude");
				longvar.Write(ComputeLonsFromLatLonGrid(grid));
			}
			catch(const int& err)
                        {
                                return false;
                        }*/

			break;
                }
		case kRotatedLatitudeLongitude:
		{
                        auto rg = dynamic_pointer_cast<rotated_latitude_longitude_grid>(grid);

                        auto r = GET_PLUGIN(radon);
                        r->RadonDB().Query("SELECT units_rlon,units_rlat,grid_mapping_name FROM cf_geom_rotated_latitude_longitude where producer_id = " + to_string(p.Id()));
                        const auto row = r->RadonDB().FetchRow();
			const string units_rlon = row[0];
			const string units_rlat = row[1];
                        const string grid_mapping = row[2];

			// create lat, lon axis
                        try
                        {
                                // add latitude dimension
                                auto lat = file.AddDim("rlat", rg->Nj());
                                file.AddVar<double>("rlat",{lat});

				// add longitude dimension
                                auto lon = file.AddDim("rlon", rg->Ni());
                                file.AddVar<double>("rlon",{lon});
                        }
                        catch(const int& err)
                        {
                                return false;
                        }

			auto lats = ComputeLatsFromLatLonGrid(grid);
			auto lons = ComputeLonsFromLatLonGrid(grid);

			// add attributes to axis
                        try
                        {
				// write longitude values
                                auto longvar = file.GetVar<double>("rlon");
				longvar.AddTextAtt("units",units_rlon);
                                longvar.Write(lons);

				// write latitude values
                                auto latvar = file.GetVar<double>("rlat");
                                latvar.AddTextAtt("units",units_rlat);
                                latvar.Write(lats);
                        }
                        catch(const int& err)
                        {
                                return false;
                        }

			// create grid_mapping variable
			try
			{
				// create rotation metadata
				file.AddVar<signed char>("rotated_pole",{});
			}
			catch(const int& err)
                        {
                                return false;
                        }

			const double SPlat = rg->SouthPole().Y();
			const double SPlon = rg->SouthPole().X();

			// add projection metadata as grid_mapping attributes
			try
			{
				nc_var<signed char> rotvar = file.GetVar<signed char>("rotated_pole");
				rotvar.AddTextAtt("grid_mapping_name",grid_mapping);
				rotvar.AddAtt<double>("grid_south_pole_latitude",SPlat);
				rotvar.AddAtt<double>("grid_south_pole_longitude",SPlon);
			}
			catch(const int& err)
                        {
				cout << "error" << nc_strerror(err) << '\n';
                                return false;
                        }

                        break;
                }
                default:
                        cout << "error\n";
	}	

	return true;
}

bool CreateMetadataLevel(nc_group& file, const level& lev, const producer& p)
{
        try
        {
                nc_dim d = file.AddDim("level", NC_UNLIMITED);
		file.AddVar<double>("level",{d});
        }
        catch(const int& err)
        {
                cout << nc_strerror(err) << '\n';
                return false;
        }

        auto r = GET_PLUGIN(radon);
        r->RadonDB().Query("SELECT units,positive FROM cf_level where producer_id = " + to_string(p.Id()) + " and level_id = " + to_string(6));
        const auto row = r->RadonDB().FetchRow();

        string levelUnits = row[0];
	string levelPositive = row[1];

        try
        {
                nc_var<double> l = file.GetVar<double>("level");
                l.AddTextAtt("units",levelUnits);
		if(!levelPositive.empty())
			l.AddTextAtt("positive",levelPositive);
        }
        catch(const int& err)
        {
                cout << nc_strerror(err) << '\n';
                return false;
        }

	return true;
}

template <typename VARTYPE>
bool CreateParam(nc_group& file, const param& par, const producer& p)
{
	try
	{
		// Create variable with dimension order as they appear in the header, i.e. in order they are created.
                file.AddVar<VARTYPE>(par.Name(),file.ListDims());
        }
	catch(const int& err)
        {
                cout << nc_strerror(err) << '\n';
                return false;
        }

	string s = "_FillValue";
	VARTYPE miss = 32700;//NC_FILL_DOUBLE;

        try
        {
                nc_var<VARTYPE> parameter = file.GetVar<VARTYPE>(par.Name());
                parameter.AddTextAtt("long_name",par.Name());
		parameter.AddAtt(s,miss);
        }
        catch(...)
        {
                std::cout << "fetch failed\n";
        }

	return true;
}

template <typename VARTYPE>
bool WriteParam(nc_group& file, const param& par, const producer& p, const matrix<VARTYPE>& data, size_t tidx)
{
	const vector<VARTYPE>& vals = data.Values();

	nc_var<VARTYPE> parameter;

	VARTYPE missing;

	if(std::is_same<VARTYPE, double>::value)
		missing = 32700;//NC_FILL_DOUBLE;

        vector<VARTYPE> arr(vals.size());
        replace_copy_if(vals.begin(), vals.end(), arr.begin(), [](const VARTYPE& val) { return himan::IsMissing(val); },
                        missing);

        try
        {
                parameter = file.GetVar<VARTYPE>(par.Name());
        }
        catch(...)
        {
                std::cout << "fetch failed\n";
        }

        try
        {
                parameter.Write(arr, {tidx,0,0,0}, {1,data.sizeZ(),data.SizeY(),data.SizeX()});
        }
        catch (...)
        {
	}

	return true;             
}

netcdf4::netcdf4()
{
        itsLogger = logger("netcdf4");
}

bool netcdf4::ToFile(info<double>& anInfo, string& outputFile, bool appendToFile)
{
        return ToFile<double>(anInfo, outputFile, appendToFile);
}

template <typename T>
bool netcdf4::ToFile(info<T>& anInfo, string& outputFile, bool appendToFile)
{
        // Write only that data which is currently set at descriptors
        nc_group theFile;

        if(appendToFile)
        {
		try
		{
               		theFile = Open(outputFile);
		}
		catch(const int& err)
		{
			//no such file or directory
			if(err == 2)
			{
				theFile = Create(outputFile);
			}
			else
			{
				itsLogger.Error("Could not open or create netcdf file");
				exit(1);
			}
		}
        }
        else
        {
		try
		{
               		theFile = Create(outputFile);
		}
		catch(const int& err)
		{
			itsLogger.Error("Could not create netcdf file");
                        exit(1);
		}
        }

	CreateMetadataTime(theFile, anInfo.Time(), anInfo.Producer());
        CreateMetadataLevel(theFile, anInfo.Level(), anInfo.Producer());
        CreateMetadataGrid(theFile, anInfo.Grid(), anInfo.Producer());

	CreateParam<T>(theFile, anInfo.Param(), anInfo.Producer());

	WriteParam<T>(theFile, anInfo.Param(), anInfo.Producer(), anInfo.Data(), anInfo.template Index<forecast_time>());

        string verb = (appendToFile ? "Appended to " : "Wrote ");
        itsLogger.Info(verb + "file '" + outputFile);

        return true;
}
template bool netcdf4::ToFile<double>(info<double>&, string&, bool);
template bool netcdf4::ToFile<float>(info<float>&, string&, bool);
