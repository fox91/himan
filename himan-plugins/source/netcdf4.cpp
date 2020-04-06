/**
 * @file netcdf4.cpp
 *
 **/

#include "common.h"
#include "dimension.h"
#include "fminc4.h"
#include "group.h"
#include "variable.h"

#include "grid.h"

#include "latitude_longitude_grid.h"
#include "logger.h"
#include "producer.h"
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
			for (size_t j = 0; j < rg->Nj(); ++j)
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
			for (size_t i = 0; i < rg->Ni(); ++i)
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
//###############################################################################
// WRITING PART
//###############################################################################

// Might scrap this and rebuild. Don't look too closely here.
#if 0
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

bool netcdf4::ToFile(info<double>& anInfo)
{
        return ToFile<double>(anInfo);
}

template <typename T>
bool netcdf4::ToFile(info<T>& anInfo)
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
template bool netcdf4::ToFile<double>(info<double>&);
template bool netcdf4::ToFile<float>(info<float>&);

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
#endif
//###############################################################################
// READING PART
//###############################################################################

// TODO Move these functions to radon db
string TimeVarname(const producer& prod)
{
	/*auto r = GET_PLUGIN(radon);
	r->RadonDB().Query("SELECT time_dimension_name FROM cf_time where producer_id = " + to_string(prod.Id()));
	auto row = r->RadonDB().FetchRow();

return row[0];*/
	return "time";
}

string TimeMask(const producer& prod)
{
	/*auto r = GET_PLUGIN(radon);
	r->RadonDB().Query("SELECT units FROM cf_time where producer_id = " + to_string(prod.Id()));
	auto row = r->RadonDB().FetchRow();

	return row[0];*/
	return "seconds since %Y-%m-%d %H:%M:%S";
}

string LevelVarname(const producer& prod, const level& level)
{
	/*auto r = GET_PLUGIN(radon);
	    r->RadonDB().Query("SELECT dimension_name FROM cf_level WHERE producer_id = "
	            + to_string(prod.Id()) + " AND level_id = (SELECT level_id FROM level_grib1 WHERE producer_id = "
	            + to_string(prod.Id()) + " AND grib_level_id = " + to_string(static_cast<int>(level.Type())) + ")");
	    auto row = r->RadonDB().FetchRow();

	    return row[0];*/
	return "";
}

string Varname(const producer& prod, const param& param)
{
	/*auto r = GET_PLUGIN(radon);
	    r->RadonDB().Query("SELECT netcdf_name FROM cf_param WHERE producer_id = " + to_string(prod.Id()) + " AND
	   param_id = ( SELECT id FROM param WHERE name = '" + param.Name() + "')"); auto row = r->RadonDB().FetchRow();

	    return row[0];*/

	return "lwe_thickness_of_surface_snow_amount";
}
// TODO END

/*
 * Read a full netcdf array into a himan matrix. This is limited to max three netcdf dimensions.
 */

template <typename T>
matrix<T> ReadData(nc_var theVar)
{
	vector<nc_dim> dims = theVar.GetDims();
	assert(dims.size() <= 3);

	size_t width = 1;
	size_t height = 1;
	size_t depth = 1;

	switch (dims.size())
	{
		case 1:
			width = dims[0].Size();
			break;
		case 2:
			width = dims[1].Size();
			height = dims[0].Size();
			break;
		case 3:
			width = dims[2].Size();
			height = dims[1].Size();
			depth = dims[0].Size();
			break;
		default:
			throw 1;  // throw some error;
	}

	matrix<T> ret(width, height, depth, MissingValue<T>());

	switch (theVar.Type())
	{
		case NC_DOUBLE:
		{
			std::vector<double> raw_data = theVar.Read<double>();
			std::transform(raw_data.begin(), raw_data.end(), ret.Values().begin(),
			               [](double c) -> T { return (c == NC_FILL_DOUBLE) ? MissingValue<T>() : static_cast<T>(c); });
			break;
		}

		case NC_FLOAT:
		{
			std::vector<float> raw_data = theVar.Read<float>();
			std::transform(raw_data.begin(), raw_data.end(), ret.Values().begin(),
			               [](float c) -> T { return (c == NC_FILL_FLOAT) ? MissingValue<T>() : static_cast<T>(c); });
			break;
		}

		default:
			break;
	}

	return ret;
}

/*
 * Read a slice of netcdf array into a himan matrix. Slice is defined by start and steps.
 */

template <typename T>
matrix<T> ReadData(nc_var theVar, vector<size_t> theStart, vector<size_t> theSteps)
{
	return matrix<T>();
}

std::unique_ptr<himan::grid> ReadGrid(nc_group theGroup, nc_var theVar)
{
	HPGridType gType = kUnknownGridType;

	auto mapping = theVar.GetAtt<std::string>("grid_mapping").front();

	nc_var mapVar = theGroup.GetVar(mapping);

	if (mapping == "latitude_longitude")
		gType = kLatitudeLongitude;

	// add other grid types

	nc_var xVar;
	nc_var yVar;

	for (auto v : theGroup.ListVars())
	{
		string axis;
		try
		{
			axis = v.GetAtt<string>("axis").front();
		}
		catch (...)
		{
			continue;
		}

		if (axis == "X")
			xVar = v;
		if (axis == "Y")
			yVar = v;
	}

	std::vector<double> xCoord = ReadData<double>(xVar).Values();
	std::vector<double> yCoord = ReadData<double>(yVar).Values();

	bool iNegative = xCoord.front() > xCoord.back();
	bool jPositive = yCoord.front() < yCoord.back();

	HPScanningMode m = kUnknownScanningMode;

	if (!iNegative && !jPositive)
	{
		m = kTopLeft;
	}
	else if (iNegative && !jPositive)
	{
		m = kTopRight;
	}
	else if (iNegative && jPositive)
	{
		m = kBottomRight;
	}
	else if (!iNegative && jPositive)
	{
		m = kBottomLeft;
	}
	else
	{
		throw runtime_error("WHAT?");
	}

	double X0 = xCoord.front();
	double Y0 = yCoord.front();

	himan::point firstPoint(X0, Y0);

	unique_ptr<grid> newGrid;

	switch (gType)
	{
		case kLatitudeLongitude:
		{
			newGrid = unique_ptr<latitude_longitude_grid>(new latitude_longitude_grid);
			latitude_longitude_grid* const rg = dynamic_cast<latitude_longitude_grid*>(newGrid.get());

			size_t ni = xCoord.size();
			size_t nj = yCoord.size();

			rg->Ni(ni);
			rg->Nj(nj);

			rg->Di(abs(xCoord.front() - xCoord.back()) / static_cast<double>(ni - 1));
			rg->Dj(abs(yCoord.front() - yCoord.back()) / static_cast<double>(nj - 1));

			rg->ScanningMode(m);

			rg->FirstPoint(firstPoint);

			double X1 = xCoord.back();
			double Y1 = yCoord.back();

			rg->LastPoint(point(X1, Y1));

			break;
		}
		default:
			// TODO throw error message of unsupported grid
			break;
	}

	// TODO probably need to check 1. does attribute exist, 2. what is type of the attribute, 3. convert type to
	// double, 4. other shapes
	double a = mapVar.GetAtt<double>("earth_radius").front();
	double b = a;
	newGrid->EarthShape(himan::earth_shape<double>(a, b));

	return newGrid;
}

/*
 * Read time from dimension attribute of 2D data
 */

forecast_time ReadTime2D(nc_group theGroup, nc_var theVar, string theMask)
{
	nc_var forecast_period;
	nc_var origin_time;

	// 2D case, time is not a dimension of the variable
	// find origin time
	vector<string> dimensions;

	string s = theVar.GetAtt<string>("coordinates").front();
	string delimiter = " ";

	size_t pos = 0;
	std::string token;
	while ((pos = s.find(delimiter)) != string::npos)
	{
		token = s.substr(0, pos);
		dimensions.push_back(token);
		s.erase(0, pos + delimiter.length());
	}

	for (auto dim : dimensions)
	{
		nc_var dimension;
		try
		{
			dimension = theGroup.GetVar(dim);
		}
		catch (...)
		{
			continue;
		}
		if (dimension.GetAtt<string>("standard_name").front() == "forecast_reference_time")
			origin_time = dimension;
		if (dimension.GetAtt<string>("standard_name").front() == "forecast_period")
			forecast_period = dimension;
	}

	string originDateTimeStr = origin_time.GetAtt<string>("units").front();
	raw_time originDateTime(originDateTimeStr, theMask);

	forecast_time t(originDateTime, originDateTime);

	uint64_t stepOrigin = origin_time.Read<uint64_t>().front();
	int stepForecast = forecast_period.Read<int>().front();

	// TODO adjust time based one units attribute of netcdf variable; kSecondResolution required?
	t.OriginDateTime().Adjust(kMinuteResolution, static_cast<int>(stepOrigin / 60));
	t.ValidDateTime().Adjust(kMinuteResolution, static_cast<int>(stepOrigin / 60) + stepForecast / 60);

	return t;
}

/*
 * Read time for given start index from time dimension of the 3/4-dimensional variable
 *
 * Not yet done
 */

forecast_time ReadTime(nc_group theGroup, nc_var theVar, vector<size_t> theIndex, string theMask)
{
	forecast_time t;
	return t;
}

/*
 * Search options return a vector of indexes to the starting points of the current variables to for example a certain
 * Level or Time step. Let's say we have a 4D variable in netcdf with 10 timesteps and 5 levels in addition to the
 * horizontal directions. Using search options matching index positions for Time and Level are returned if they are
 * present as variable dimensions. TODO !under construction!
 */

std::pair<std::vector<size_t>, std::vector<size_t>> ApplySearchOptions(nc_group& theFile, const search_options& options)
{
	nc_var var;
	try
	{
		var = theFile.GetVar(Varname(options.prod, options.param));
	}
	catch (const int& err)
	{
		// itsLogger.Error("Opening param '" + options.param.Name() + "' failed: " + nc_strerror(err) );
		return std::make_pair(std::vector<size_t>{}, std::vector<size_t>{});
	}

	auto dims = var.GetDims();
	size_t nDims = dims.size();
	std::vector<size_t> ret{nDims};

	// limit to 4 Dims
	if (nDims > 4)
	{
		// itsLogger.Error("Opening param '" + options.param.Name() + "' failed: " + nc_strerror(err) );
		return std::make_pair(std::vector<size_t>{}, std::vector<size_t>{});
	}

	// find time dim
	nc_var timevar;
	std::pair<int, size_t> dimnum_time{0, 0};
	bool time_found = false;
	for (auto dim : dims)
	{
		if (dim.Name() == TimeVarname(options.prod))
		{
			timevar = theFile.GetVar(dim.Name());
			auto times = timevar.Read<long>();

			raw_time origin_time(timevar.GetAtt<string>("units").front(), TimeMask(options.prod));
			auto ti = find_if(times.begin(), times.end(), [=](long l) {
				return options.time == forecast_time(origin_time, time_duration(kHourResolution, l));
			});
			dimnum_time.second = ti - times.begin();

			time_found = true;
			break;
		}
		++dimnum_time.first;
	}

	// find level dim
	nc_var levelvar;
	std::pair<int, size_t> dimnum_level{0, 0};
	bool level_found = false;
	for (auto dim : dims)
	{
		if (dim.Name() == LevelVarname(options.prod, options.level))
		{
			levelvar = theFile.GetVar(dim.Name());

			// TODO check level type from units attribute

			auto levels = levelvar.Read<double>();

			auto li = find_if(levels.begin(), levels.end(), [=](double l) { return options.level.Value() == l; });
			dimnum_level.second = li - levels.begin();

			level_found = true;
			break;
		}
		++dimnum_level.first;
	}

	if (time_found)
		ret[dimnum_time.first] = dimnum_time.second;
	if (level_found)
		ret[dimnum_level.first] = dimnum_level.second;

	// TODO, set correct values for ny and nx
	size_t ny = 1;
	size_t nx = 1;
	return std::make_pair(ret, std::vector<size_t>{1, 1, ny, nx});
}

template <typename T>
vector<shared_ptr<himan::info<T>>> netcdf4::FromFile(const file_information& theInputFile,
                                                     const search_options& options) const
{
	vector<shared_ptr<himan::info<T>>> infos;

	// open the file
	nc_group theFile;
	try
	{
		theFile = Open(theInputFile.file_location);
		itsLogger.Info("Opening file " + theInputFile.file_location);
	}
	catch (const int& err)
	{
		itsLogger.Error("Opening file '" + theInputFile.file_location + "' failed: " + nc_strerror(err));
		return infos;
	}

	// try open variable
	nc_var var;
	try
	{
		var = theFile.GetVar(Varname(options.prod, options.param));
	}
	catch (const int& err)
	{
		itsLogger.Error("Opening param '" + options.param.Name() + "' failed: " + nc_strerror(err));
		return infos;
	}

	// read variable data
	matrix<T> data;
	try
	{
		// if var has only two dimensions assume is single level, single time forecast field.
		if (var.GetDims().size() == 2)
		{
			data = ReadData<T>(var);
		}
		else
		{
			// Data is contains multiple time steps/levels. Find starting point for 2D field from search options and
			// download.
			// TODO Have to verify this later. Comment out for now
			// auto access_point = ApplySearchOptions(theFile, options);
			// data = ReadData(var, access_point.first, access_point.second);
		}
	}
	catch (const int& err)
	{
		itsLogger.Error("Reading data '" + options.param.Name() + "' failed: " + nc_strerror(err));
		return infos;
	}

	// validate origin time & forecast time
	// TODO Read time for non-2D data
	if (ReadTime2D(theFile, var, TimeMask(options.prod)) != options.time)
	{
		itsLogger.Error("Reading data time failed");
		return infos;
	}

	// TODO validate level

	// set info metadata
	auto newInfo = std::make_shared<himan::info<T>>();

	newInfo->Producer(options.prod);
	newInfo->template Set<param>({options.param});
	newInfo->template Set<forecast_time>({options.time});
	newInfo->template Set<level>({options.level});
	newInfo->template Set<forecast_type>({options.ftype});

	unique_ptr<grid> newGrid = ReadGrid(theFile, var);
	// latitude_longitude_grid* const rg = dynamic_cast<latitude_longitude_grid*>(newGrid.get());

	auto b = make_shared<base<T>>();
	b->grid = shared_ptr<grid>(newGrid->Clone());
	b->data = data;

	newInfo->Create(b, true);

	// Set descriptors

	newInfo->template Find<param>(options.param);
	newInfo->template Find<forecast_time>(options.time);
	newInfo->template Find<level>(options.level);
	newInfo->template Find<forecast_type>(options.ftype);

	infos.push_back(newInfo);

	return infos;
}
template vector<shared_ptr<himan::info<double>>> netcdf4::FromFile<double>(const file_information&,
                                                                           const search_options&) const;
template vector<shared_ptr<himan::info<float>>> netcdf4::FromFile<float>(const file_information&,
                                                                         const search_options&) const;
