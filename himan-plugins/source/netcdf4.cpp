/**
 * @file netcdf4.cpp
 *
 **/

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

bool CreateMetadataTime(nc_group& file, const vector<forecast_time>& ftimes, const producer& p)
{
        // create netcdf dimension and according variable
        /*try
        {
		const auto row = r->RadonDB().FetchRow();
		string tDim = row[0];

                nc_dim d = file.AddDim("T", ftimes.size());
                file.AddVar<long>("T",{d});
        }
        catch(const int& err)
        {
                return false;
        }*/

        auto r = GET_PLUGIN(radon);
        r->RadonDB().Query("SELECT units,time_dimension_name FROM cf_time where producer_id = " + to_string(p.Id()));
        const auto row = r->RadonDB().FetchRow();

        string tMask = row[0];
	string tDim = row[1];

        try
        {
		nc_dim d = file.AddDim(tDim, ftimes.size());
                nc_var<long> t = file.AddVar<long>(tDim,{d});
                t.AddTextAtt("long_name","time");
                t.AddTextAtt("units",ftimes.front().OriginDateTime().String(tMask));
        }
        catch(const int& err)
        {
                return false;
        }

	vector<long> steps;
	steps.reserve(ftimes.size());

	for(auto ftime : ftimes)
	{
		steps.push_back(ftime.Step().Hours());
	}

        try
        {
                nc_var<long> t = file.GetVar<long>(tDim);
		t.Write(steps);
        }
        catch(const int& err)
        {
                return false;
        }

        return true;
}

bool CreateMetadataGrid(nc_group& file, const shared_ptr<himan::grid>& grid, const producer& p)
{
	switch (grid->Type())
        {
		case kRotatedLatitudeLongitude:
		{
                        auto rg = dynamic_pointer_cast<rotated_latitude_longitude_grid>(grid);

                        auto r = GET_PLUGIN(radon);
                        r->RadonDB().Query("SELECT units_rlon,units_rlat,grid_mapping_name,lon_dimension_name,lat_dimension_name FROM cf_geom_rotated_latitude_longitude where producer_id = " + to_string(p.Id()));
                        const auto row = r->RadonDB().FetchRow();
			const string units_rlon = row[0];
			const string units_rlat = row[1];
                        const string grid_mapping = row[2];
			const string lon_dim_name = row[3];
			const string lat_dim_name = row[4];

			// create lat, lon axis
                        try
                        {
                                // add latitude dimension
                                auto lat = file.AddDim(lat_dim_name, rg->Nj());
                                file.AddVar<double>(lat_dim_name,{lat});

				// add longitude dimension
                                auto lon = file.AddDim(lon_dim_name, rg->Ni());
                                file.AddVar<double>(lon_dim_name,{lon});
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
                                auto longvar = file.GetVar<double>(lon_dim_name);
				longvar.AddTextAtt("units",units_rlon);
                                longvar.Write(lons);

				// write latitude values
                                auto latvar = file.GetVar<double>(lat_dim_name);
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
                                return false;
                        }

			// add coordinate reference system
		        try
        		{
               			// create rotation metadata
               			file.AddVar<int>("crs",{});
        		}
        		catch(const int& err)
        		{
                		return false;
        		}

                        try
                        {
                                nc_var<int> crs = file.GetVar<int>("crs");
                                crs.AddTextAtt("grid_mapping_name",grid_mapping);
                                crs.AddAtt<double>("semi_major_axis",rg->EarthShape().A());
                                crs.AddAtt<double>("inverse_flattening",rg->EarthShape().F());
                        }
			catch(const int& err)
                        {
                                return false;
                        }

                        break;
                }
                default:
                        cout << "error\n";
	}	

	return true;
}

bool CreateMetadataLevel(nc_group& file, const vector<level>& lev, const producer& p)
{
	nc_dim d;
        try
        {
                d = file.AddDim("level", lev.size());
        }
        catch(const int& err)
        {
                return false;
        }

	switch(lev.front().Type())
	{
		case kHeight: 
		{
        		auto r = GET_PLUGIN(radon);
        		r->RadonDB().Query("SELECT units,positive,dimension_name FROM cf_level where producer_id = " + to_string(p.Id()) + " and level_id = " + to_string(6));
        		const auto row = r->RadonDB().FetchRow();

        		string levelUnits = row[0];
			string levelPositive = row[1];
			string levelName = row[2];

        		try
        		{
				d.Name(levelName);
                		nc_var<double> l = file.AddVar<double>(levelName,{d});;
                		l.AddTextAtt("units",levelUnits);
				if(!levelPositive.empty())
					l.AddTextAtt("positive",levelPositive);
        		}
        		catch(const int& err)
        		{
                		return false;
        		}

                        vector<double> lvlvalue;
                        for(auto value : lev)
                        {
                                lvlvalue.push_back(value.Value());
                        }

                        try
                        {
                                nc_var<double> t = file.GetVar<double>(levelName);
                                t.Write(lvlvalue);
                        }
                        catch(const int& err)
                        {
                                return false;
                        }


			break;
		}
		case kHybrid:
		{
        		try
        		{
                		nc_var<double> l = file.AddVar<double>("level",{d});
                		l.AddTextAtt("standard_name","atmosphere_hybrid_sigma_pressure_coordinate");
				l.AddTextAtt("formula_terms","ak: A b: B ps: Ps");
			
				file.AddVar<double>("A",{d});
				file.AddVar<double>("B",{d});
        		}
        		catch(const int& err)
        		{
                		return false;
        		}

			vector<double> AB = lev.front().AB();
			vector<double> lvlvalue;
			vector<double> A;
			vector<double> B;
			for(auto value : lev)
			{
				lvlvalue.push_back(value.Value());
				//A.push_back(AB[static_cast<size_t>(value.Value())]);
				//B.push_back(AB[static_cast<size_t>(value.Value()) + 137]);
			}

                        try
                        {
                                nc_var<double> t = file.GetVar<double>("level");
                                t.Write(lvlvalue);

                                //nc_var<double> a = file.GetVar<double>("A");
                                //nc_var<double> b = file.GetVar<double>("B");

                               	//a.Write(A);
                               	//b.Write(B);
                        }
                        catch(const int& err)
                        {
                                return false;
                        }

		        break;	
		}
		default:
			break;
	}

	return true;
}

template <typename VARTYPE>
bool CreateParam(nc_group& file, const param& par, const producer& p)
{
        auto r = GET_PLUGIN(radon);

        r->RadonDB().Query("SELECT netcdf_name FROM cf_param WHERE producer_id = " + to_string(p.Id()) + " AND param_id = ( SELECT id FROM param WHERE name = '" + par.Name() + "')");
        auto row = r->RadonDB().FetchRow();
        string parm = row[0];

	try
	{
		// Create variable with dimension order as they appear in the header, i.e. in order they are created.
                file.AddVar<VARTYPE>(parm,file.ListDims());
        }
	catch(const int& err)
        {
                return false;
        }

	string s = "_FillValue";
	VARTYPE missing;
	if(std::is_same<VARTYPE, double>::value)
        {
                missing = NC_FILL_DOUBLE;
        }
        else if(std::is_same<VARTYPE, float>::value)
        {
                missing = NC_FILL_FLOAT;
        }

	string crs = "crs";

        try
        {
                nc_var<VARTYPE> parameter = file.GetVar<VARTYPE>(parm);
                parameter.AddTextAtt("long_name",par.Name());
		parameter.AddAtt(s,missing);
		parameter.AddTextAtt("grid_mapping",crs);
        }
        catch(...)
        {
                std::cout << "fetch failed\n";
        }

	return true;
}

template <typename VARTYPE>
bool WriteParam(nc_group& file, const param& par, const producer& p, const matrix<VARTYPE>& data, size_t tidx, size_t zidx)
{
        auto r = GET_PLUGIN(radon);

        r->RadonDB().Query("SELECT netcdf_name FROM cf_param WHERE producer_id = " + to_string(p.Id()) + " AND param_id = ( SELECT id FROM param WHERE name = '" + par.Name() + "')");
        auto row = r->RadonDB().FetchRow();
        string parm = row[0];

	const vector<VARTYPE>& vals = data.Values();

	nc_var<VARTYPE> parameter;

	VARTYPE missing;

	if(std::is_same<VARTYPE, double>::value)
	{
		missing = NC_FILL_DOUBLE;
	}
	else if(std::is_same<VARTYPE, float>::value)
	{
		missing = NC_FILL_FLOAT;
	}

        vector<VARTYPE> arr(vals.size());
        replace_copy_if(vals.begin(), vals.end(), arr.begin(), [](const VARTYPE& val) { return himan::IsMissing(val); },
                        missing);

        try
        {
                parameter = file.GetVar<VARTYPE>(parm);
        }
        catch(const int& err)
        {
		return false;
        }

        try
        {
                parameter.Write(arr, {tidx,zidx,0,0}, {1,data.SizeZ(),data.SizeY(),data.SizeX()});
        }
        catch (const int& err)
        {
		return false;
	}

	return true;             
}

bool WriteLevel(nc_group& file, const level& lev, const producer& p, size_t lidx)
{
	switch(lev.Type())
	{
		case kHeight:
		{
        		try
        		{
                		nc_var<double> t = file.GetVar<double>("level");
                		t.Write(lev.Value(), {lidx});
        		}
        		catch(const int& err)
        		{
                		return false;
        		}
			break;
		}
		case kHybrid:
		{
			try
        		{
                		nc_var<double> t = file.GetVar<double>("level");
                		t.Write(lev.Value(), {lidx});

				nc_var<double> a = file.GetVar<double>("A");
				nc_var<double> b = file.GetVar<double>("B");

				a.Write(lev.AB()[static_cast<size_t>(lev.Value())], {lidx});
				b.Write(lev.AB()[static_cast<size_t>(lev.Value()) + 137], {lidx});
        		}
        		catch(const int& err)
        		{
                		return false;
        		}	
        		break;
		}
		default:
			break;
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

        try
        {
	        theFile = Open(outputFile);
        }
        catch(const int& err)
        {
                itsLogger.Error(nc_strerror(err));
                exit(1);
        }

	WriteParam<T>(theFile, anInfo.Param(), anInfo.Producer(), anInfo.Data(), anInfo.template Index<forecast_time>(), anInfo.template Index<level>());

        string verb = (appendToFile ? "Appended to " : "Wrote ");
        itsLogger.Info(verb + "file '" + outputFile);

        return true;
}
template bool netcdf4::ToFile<double>(info<double>&, string&, bool);
template bool netcdf4::ToFile<float>(info<float>&, string&, bool);

template <typename T>
bool netcdf4::InitFile(const info<T>& anInfo, const string& outputFile)
{
        nc_group theFile = Create(outputFile);
	CreateMetadataTime(theFile, anInfo.template Iterator<forecast_time>().Values(), anInfo.Producer());
        CreateMetadataLevel(theFile, anInfo.template Iterator<level>().Values(), anInfo.Producer());
        CreateMetadataGrid(theFile, anInfo.Grid(), anInfo.Producer());

	for(auto parameter : anInfo.template Iterator<param>().Values())
		CreateParam<T>(theFile, parameter, anInfo.Producer());

        return true;
}
template bool netcdf4::InitFile<double>(const info<double>&, const string&);
template bool netcdf4::InitFile<float>(const info<float>&, const string&);

//###############################################################################
// READING PART
//###############################################################################

template <typename T>
vector<shared_ptr<himan::info<T>>> netcdf4::FromFile(const string& theInputFile, const search_options& options) const
{
        vector<shared_ptr<himan::info<T>>> infos;

	nc_group theFile;
	try
	{
		theFile = Open(theInputFile);
		itsLogger.Info("Opening file " + theInputFile);
	}
	catch(const int& err)
	{
                itsLogger.Error("Opening file '" + theInputFile + "' failed: " + nc_strerror(err) );
                return infos;
        }

        auto r = GET_PLUGIN(radon);
        r->RadonDB().Query("SELECT units,time_dimension_name FROM cf_time where producer_id = " + to_string(options.prod.Id()));
        auto row = r->RadonDB().FetchRow();

        string tMask = row[0];
        string tDim = row[1];

	auto theTime = theFile.GetVar<long>(tDim);
	auto times = theTime.Read();

	raw_time origin_time(theTime.GetAtt<string>("units").front(),tMask);
	auto ti = find_if(times.begin(),times.end(),[=](long l){return options.time == forecast_time(origin_time,time_duration(kHourResolution,l));});
	size_t timeIdx = ti - times.begin();

        r->RadonDB().Query("SELECT dimension_name FROM cf_level WHERE producer_id = " + to_string(options.prod.Id()) + " AND level_id = (SELECT level_id FROM level_grib1 WHERE producer_id = " + to_string(options.prod.Id()) + " AND grib_level_id = " + to_string(static_cast<int>(options.level.Type())) + ")");
	row = r->RadonDB().FetchRow();

	string lDim = row[0];

	auto theLevel = theFile.GetVar<double>(lDim);
	auto levels = theLevel.Read();
	
	auto li = find(levels.begin(),levels.end(),options.level.Value());
	size_t lvlIdx = li - levels.begin();

	r->RadonDB().Query("SELECT netcdf_name FROM cf_param WHERE producer_id = " + to_string(options.prod.Id()) + " AND param_id = ( SELECT id FROM param WHERE name = '" + options.param.Name() + "')");
	row = r->RadonDB().FetchRow();
	string par = row[0];

	auto pot = theFile.GetVar<T>(par);

	auto vals = pot.Read({timeIdx,lvlIdx,0,0},{1,1,576,661});


	himan::info<T> ifo;
	ifo.Producer(options.prod);
	ifo.template Set<param>({options.param});
	ifo.template Set<forecast_time>({options.time});
	ifo.template Set<level>({options.level});
	ifo.template Set<forecast_type>({options.ftype});
	rotated_latitude_longitude_grid newGrid;

	auto pole = theFile.GetVar<signed char>("rotated_pole");
	double lat = pole.GetAtt<double>("grid_south_pole_latitude").front();
        double lon = pole.GetAtt<double>("grid_south_pole_longitude").front();
			newGrid.Ni(661);
			newGrid.Nj(576);

			newGrid.Di(0.1);
			newGrid.Dj(0.1);

			newGrid.SouthPole(himan::point(lon, lat));
			newGrid.UVRelativeToGrid(true);

			newGrid.ScanningMode(kTopLeft);
			newGrid.FirstPoint(himan::point(-26, 22.5));
  			newGrid.LastPoint(himan::point(40, -35));

	auto b = make_shared<base<T>>();
	b->grid = shared_ptr<grid>(make_shared<rotated_latitude_longitude_grid>(newGrid));

	ifo.Create(b,true);

	ifo.template Find<param>(options.param);
	ifo.template Find<forecast_time>(options.time);
	ifo.template Find<level>(options.level);
	ifo.template Find<forecast_type>(options.ftype);
	ifo.Data().Set(vals);

	infos.push_back(make_shared<himan::info<T>>(ifo));

        return infos;
}
template vector<shared_ptr<himan::info<double>>> netcdf4::FromFile<double>(const string&, const search_options&) const;
template vector<shared_ptr<himan::info<float>>> netcdf4::FromFile<float>(const string&, const search_options&) const;


