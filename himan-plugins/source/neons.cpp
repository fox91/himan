/**
 * @file neons.cpp
 *
 * @date Nov 20, 2012
 * @author partio
 */

#include "neons.h"
#include "logger_factory.h"
#include "plugin_factory.h"
#include <thread>
#include <sstream>
#include "util.h"

using namespace std;
using namespace himan::plugin;

const int MAX_WORKERS = 16;
once_flag oflag;

neons::neons() : itsInit(false), itsNeonsDB()
{
	itsLogger = unique_ptr<logger> (logger_factory::Instance()->GetLog("neons"));
	
	// no lambda functions for gcc 4.4 :(
	// call_once(oflag, [](){ NFmiNeonsDBPool::MaxWorkers(MAX_WORKERS); });

	call_once(oflag, &himan::plugin::neons::InitPool, this);
}

void neons::InitPool()
{
	NFmiNeonsDBPool::Instance()->MaxWorkers(MAX_WORKERS);

	char* host;

	host = getenv("HOSTNAME");

	// TODO: a smarter way to check if we are in production server

	if (host != NULL)
	{
		if (string(host) == "gogol.fmi.fi" || string(host) == "tolstoi.fmi.fi")
		{
			NFmiNeonsDBPool::Instance()->ExternalAuthentication(true);
			NFmiNeonsDBPool::Instance()->ReadWriteTransaction(true);
		}
		else if (string(host) == "jeera.fmi.fi" || string(host) == "sahrami.fmi.fi")
		{
			NFmiNeonsDBPool::Instance()->ReadWriteTransaction(true);
			NFmiNeonsDBPool::Instance()->Username("wetodb");
			NFmiNeonsDBPool::Instance()->Password("3loHRgdio");
		}
	}
}

vector<string> neons::Files(const search_options& options)
{

	Init();

	vector<string> files;

	string analtime = options.time.OriginDateTime()->String("%Y%m%d%H%M%S");
	string levelvalue = boost::lexical_cast<string> (options.level.Value());

	string ref_prod = options.configuration->SourceProducer().Name(); //producerInfo["ref_prod"];
//	long proddef = options.configuration->SourceProducer().Id(); // producerInfo["producer_id"];
	long no_vers = options.configuration->SourceProducer().TableVersion(); // producerInfo["no_vers"];

	string level_name = options.level.Name();

	vector<vector<string> > gridgeoms = itsNeonsDB->GetGridGeoms(ref_prod, analtime, options.configuration->SourceGeomName());

	if (gridgeoms.empty())
	{
		itsLogger->Warning("No geometries found for given search options");
		return files;
	}

	for (size_t i = 0; i < gridgeoms.size(); i++)
	{
		string tablename = gridgeoms[i][1];
		string dset = gridgeoms[i][2];

		/// @todo GFS (or in fact codetable 2) has wrong temperature parameter defined

		string parm_name = options.param.Name();

		if (parm_name == "T-K" && no_vers == 2)
		{
		 	parm_name = "T-C";
		}

		string query = "SELECT parm_name, lvl_type, lvl1_lvl2, fcst_per, file_location, file_server "
				   "FROM "+tablename+" "
				   "WHERE dset_id = "+dset+" "
				   "AND parm_name = upper('"+parm_name+"') "
				   "AND lvl_type = upper('"+level_name+"') "
				   "AND lvl1_lvl2 = " +levelvalue+" "
				   "AND fcst_per = "+boost::lexical_cast<string> (options.time.Step())+" "
				   "ORDER BY dset_id, fcst_per, lvl_type, lvl1_lvl2";

		itsNeonsDB->Query(query);

		vector<string> values = itsNeonsDB->FetchRow();

		if (values.empty())
		{
			continue;
		}

		itsLogger->Trace("Found data for parameter " + parm_name + " from neons geometry " + gridgeoms[i][0]);
		
		files.push_back(values[4]);

	}

	return files;

}

bool neons::Save(shared_ptr<const info> resultInfo, const string& theFileName)
{
	Init();

	stringstream query;

	/*
	 * 1. Get grid information
	 * 2. Get model information
	 * 3. Get data set information (ie model run)
	 * 4. Insert or update
	 */

	himan::point firstGridPoint = resultInfo->Grid()->FirstGridPoint();

	/*
	 * pas_latitude and pas_longitude cannot be checked programmatically
	 * since f.ex. in the case for GFS in neons we have value 500 and
	 * by calculating we have value 498. But not check these columns should
	 * not matter as long as row_cnt, col_cnt, lat_orig and lon_orig match
	 * (since pas_latitude and pas_longitude are derived from these anyway)
	 */

	query 	<< "SELECT geom_name "
			<< "FROM grid_reg_geom "
			<< "WHERE row_cnt = " << resultInfo->Nj()
			<< " AND col_cnt = " << resultInfo->Ni()
			<< " AND lat_orig = " << (firstGridPoint.Y() * 1e3)
			<< " AND long_orig = " << (firstGridPoint.X() * 1e3);
//			<< " AND pas_latitude = " << static_cast<long> (resultInfo->Dj() * 1e3)
//			<< " AND pas_longitude = " << static_cast<long> (resultInfo->Di() * 1e3);

	itsNeonsDB->Query(query.str());

	vector<string> row;

	row = itsNeonsDB->FetchRow();

	if (row.empty())
	{
		itsLogger->Warning("Grid geometry not found from neons");
		return false;
	}

	string geom_name = row[0];

	query.str("");

	query 	<< "SELECT "
			<< "nu.model_id AS process, "
			<< "nu.ident_id AS centre, "
			<< "m.model_name, "
			<< "m.model_type, "
			<< "type_smt "
			<< "FROM "
			<< "grid_num_model_grib nu, "
			<< "grid_model m, "
			<< "grid_model_name na, "
			<< "fmi_producers f "
			<< "WHERE f.producer_id = " << resultInfo->Producer().Id()
			<< " AND m.model_type = f.ref_prod "
			<< " AND nu.model_name = m.model_name "
			<< " AND m.flag_mod = 0 "
			<< " AND nu.model_name = na.model_name "
			<< " AND m.model_name = na.model_name ";

	itsNeonsDB->Query(query.str());

	row = itsNeonsDB->FetchRow();

	if (row.empty())
	{
		itsLogger->Warning("Producer definition not found from neons (id: " + boost::lexical_cast<string> (resultInfo->Producer().Id()) + ")");
		return false;
	}

	string process = row[0];
	string centre = row[1];
	string model_name = row[2];
	string model_type = row[3];

	/*
		query 	<< "SELECT "
				<< "m.model_name, "
				<< "model_type, "
				<< "type_smt "
				<< "FROM grid_num_model_grib nu, "
				<< "grid_model m, "
				<< "grid_model_name na "
				<< "WHERE nu.model_id = " << info.process
				<< " AND nu.ident_id = " << info.centre
				<< " AND m.flag_mod = 0 "
				<< " AND nu.model_name = na.model_name "
				<< " AND m.model_name = na.model_name";

	*/

	query.str("");

	query	<< "SELECT "
			<< "dset_id, "
			<< "table_name, "
			<< "rec_cnt_dset "
			<< "FROM as_grid "
			<< "WHERE "
			<< "model_type = '" << model_type << "'"
			<< " AND geom_name = '" << geom_name << "'"
			<< " AND dset_name = 'AF'"
			<< " AND base_date = '" << resultInfo->OriginDateTime().String("%Y%m%d%H%M") << "'";

	itsNeonsDB->Query(query.str());

	row = itsNeonsDB->FetchRow();

	if (row.empty())
	{
		itsLogger->Warning("Data set definition not found from neons");
		return false;
	}

	string table_name = row[1];
	string dset_id = row[0];

	string eps_specifier = "0";

	query.str("");

	query << "UPDATE as_grid "
		   << "SET rec_cnt_dset = "
		   << "rec_cnt_dset + 1, "
		   << "date_maj_dset = sysdate "
		   << "WHERE dset_id = " << dset_id;

	try
	{
		itsNeonsDB->Execute(query.str());
	}
	catch (int e)
	{
	   itsLogger->Error("Error code: " + boost::lexical_cast<string> (e));
	   itsLogger->Error("Query: " + query.str());
   	   itsNeonsDB->Rollback();
   	   return false;
	}

	query.str("");

	string host = "himan_test_host";

	char* kone = getenv("HOST");

	if (kone != NULL)
	{
		host = string(kone);
	}

	query  << "INSERT INTO " << table_name
		   << " (dset_id, parm_name, lvl_type, lvl1_lvl2, fcst_per, eps_specifier, file_location, file_server) "
		   << "VALUES ("
		   << dset_id << ", "
		   << "'" << resultInfo->Param().Name() << "', "
		   << "upper('" << resultInfo->Level().Name() << "'), "
		   << resultInfo->Level().Value() << ", "
		   << resultInfo->Time().Step() << ", "
		   << "'" << eps_specifier << "', "
		   << "'" << theFileName << "', "
		   << "'" << host << "')";


	try
	{
		itsNeonsDB->Execute(query.str());
		itsNeonsDB->Commit();
	}
	catch (int e)
	{
		itsLogger->Error("Error code: " + boost::lexical_cast<string> (e));
		itsLogger->Error("Query: " + query.str());
		itsNeonsDB->Rollback();
		return false;
	}

	itsLogger->Info("Saved information on file '" + theFileName + "' to neons");

	return true;
}

string neons::GribParameterName(const long fmiParameterId, const long codeTableVersion)
{	
	Init();
	
	string paramName = itsNeonsDB->GetGridParameterName(fmiParameterId, codeTableVersion, codeTableVersion);
	return paramName; 
	
}

string neons::GribParameterName(const long fmiParameterId, const long category, const long discipline, const long producer)
{
	Init();
	
	string paramName = itsNeonsDB->GetGridParameterName(fmiParameterId, category, discipline, producer);
	return paramName;   
}
