/**
 * @file json_parser.cpp
 *
 */
#include "interpolate.h"
#include "json_parser.h"
#include "lambert_conformal_grid.h"
#include "latitude_longitude_grid.h"
#include "plugin_factory.h"
#include "point.h"
#include "point_list.h"
#include "reduced_gaussian_grid.h"
#include "stereographic_grid.h"
#include "util.h"
#include <boost/algorithm/string/trim.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <map>
#include <stdexcept>
#include <utility>

#define HIMAN_AUXILIARY_INCLUDE

#include "cache.h"
#include "radon.h"

#undef HIMAN_AUXILIARY_INCLUDE

using namespace himan;
using namespace std;

unique_ptr<grid> ParseAreaAndGridFromPoints(configuration& conf, const boost::property_tree::ptree& pt);
unique_ptr<grid> ParseAreaAndGridFromDatabase(configuration& conf, const boost::property_tree::ptree& pt);
void ParseLevels(shared_ptr<info> anInfo, const boost::property_tree::ptree& pt);
void ParseSourceProducer(shared_ptr<configuration> conf, shared_ptr<info> anInfo,
                         const boost::property_tree::ptree& pt);
void ParseTargetProducer(shared_ptr<configuration> conf, shared_ptr<info> anInfo,
                         const boost::property_tree::ptree& pt);
vector<forecast_type> ParseForecastTypes(const boost::property_tree::ptree& pt);

vector<level> LevelsFromString(const string& levelType, const string& levelValues);

static logger itsLogger;

/*
 * Parse()
 *
 * Read command line options and create info instance.
 *
 * Steps taken:
 *
 * 1) Read command line options. Options specified in command line will
 *	override those in the conf file.
 *
 * 2) Read configuration file (if specified).
 *
 * 3) Create configuration instance.
 *
 * Some of the required information is missing, this function will not
 * behave nicely and will throw an error.
 *
 */

json_parser::json_parser()
{
	itsLogger = logger("json_parser");
}
vector<shared_ptr<plugin_configuration>> json_parser::Parse(shared_ptr<configuration> conf)
{
	if (conf->ConfigurationFile().empty())
	{
		throw runtime_error("Configuration file not defined");
	}

	vector<shared_ptr<plugin_configuration>> plugins = ParseConfigurationFile(conf);

	if (plugins.size() == 0)
	{
		throw runtime_error("Empty processqueue");
	}

	return plugins;
}

vector<shared_ptr<plugin_configuration>> json_parser::ParseConfigurationFile(shared_ptr<configuration> conf)
{
	itsLogger.Trace("Parsing configuration file '" + conf->ConfigurationFile() + "'");

	boost::property_tree::ptree pt;

	try
	{
		boost::property_tree::json_parser::read_json(conf->ConfigurationFile(), pt);
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error reading configuration file: ") + e.what());
	}

	vector<shared_ptr<plugin_configuration>> pluginContainer;
	/* Create our base info */

	auto baseInfo = make_shared<info>();

	/* Check producers */

	ParseSourceProducer(conf, baseInfo, pt);
	ParseTargetProducer(conf, baseInfo, pt);

	/* Check area definitions */

	auto g = ParseAreaAndGrid(conf, pt);

	baseInfo->itsBaseGrid = move(g);

	/* Check time definitions */

	conf->FirstSourceProducer();
	ParseTime(conf, baseInfo, pt);

	/* Check levels */

	// ParseLevels(baseInfo, pt);

	/* Check file_write */

	try
	{
		string theFileWriteOption = pt.get<string>("file_write");

		if (theFileWriteOption == "database")
		{
			conf->FileWriteOption(kDatabase);
		}
		else if (theFileWriteOption == "single")
		{
			conf->FileWriteOption(kSingleFile);
		}
		else if (theFileWriteOption == "multiple")
		{
			conf->FileWriteOption(kMultipleFiles);
		}
		else if (theFileWriteOption == "cache only")
		{
			conf->FileWriteOption(kCacheOnly);
		}
		else
		{
			throw runtime_error("Invalid value for file_write: " + theFileWriteOption);
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing key file_write: ") + e.what());
	}

	/* Check file_compression */

	try
	{
		string theFileCompression = pt.get<string>("file_compression");

		if (theFileCompression == "gzip")
		{
			conf->FileCompression(kGZIP);
		}
		else if (theFileCompression == "bzip2")
		{
			conf->FileCompression(kBZIP2);
		}
		else
		{
			conf->FileCompression(kNoCompression);
		}

		if (conf->FileCompression() != kNoCompression && conf->FileWriteOption() == kSingleFile)
		{
			itsLogger.Warning(
			    "file_write_option value 'single' conflicts with file_compression, using 'multiple' instead");
			conf->FileWriteOption(kMultipleFiles);
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing key file_write: ") + e.what());
	}

	/* Check read_data_from_database */

	try
	{
		string theReadDataFromDatabase = pt.get<string>("read_data_from_database");

		if (!util::ParseBoolean(theReadDataFromDatabase) || conf->DatabaseType() == kNoDatabase)
		{
			conf->ReadDataFromDatabase(false);
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing key read_data_from_database: ") + e.what());
	}

	// Check global use_cache option

	try
	{
		string theUseCache = pt.get<string>("use_cache");

		if (!util::ParseBoolean(theUseCache))
		{
			conf->UseCache(false);
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing key use_cache: ") + e.what());
	}

	// Check global cache_limit option

	try
	{
		int theCacheLimit = pt.get<int>("cache_limit");

		if (theCacheLimit < 1)
		{
			itsLogger.Warning("cache_limit must be larger than 0");
		}
		else
		{
			conf->CacheLimit(theCacheLimit);
			plugin::cache_pool::Instance()->CacheLimit(theCacheLimit);
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing key use_cache: ") + e.what());
	}

	// Check global file_type option

	try
	{
		string theFileType = boost::to_upper_copy(pt.get<string>("file_type"));

		if (theFileType == "GRIB")
		{
			conf->itsOutputFileType = kGRIB;
		}
		else if (theFileType == "GRIB1")
		{
			conf->itsOutputFileType = kGRIB1;
		}
		else if (theFileType == "GRIB2")
		{
			conf->itsOutputFileType = kGRIB2;
		}
		else if (theFileType == "FQD" || theFileType == "QUERYDATA")
		{
			conf->itsOutputFileType = kQueryData;
		}
		else
		{
			throw runtime_error("Invalid option for 'file_type': " + theFileType);
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing key file_type: ") + e.what());
	}

	// Check global forecast_type option

	auto forecastTypes = ParseForecastTypes(pt);

	if (forecastTypes.empty())
	{
		// Default to deterministic
		forecastTypes.push_back(forecast_type(kDeterministic));
	}

	baseInfo->Set<forecast_type>(forecastTypes);

	/* Check dynamic_memory_allocation */

	try
	{
		string theUseDynamicMemoryAllocation = pt.get<string>("dynamic_memory_allocation");

		if (util::ParseBoolean(theUseDynamicMemoryAllocation))
		{
			conf->UseDynamicMemoryAllocation(true);
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing key file_write: ") + e.what());
	}

	/*
	 * Check processqueue.
	 *
	 * Some configuration elements might be replicated here; if so they will overwrite
	 * those specified in the upper level.
	 */

	boost::property_tree::ptree& pq = pt.get_child("processqueue");

	if (pq.empty())
	{
		throw runtime_error(ClassName() + ": processqueue missing");
	}

	for (boost::property_tree::ptree::value_type& element : pq)
	{
		auto anInfo = make_shared<info>(*baseInfo);

		try
		{
			ParseTime(conf, anInfo, element.second);
		}
		catch (...)
		{
			// do nothing
		}

		try
		{
			g = ParseAreaAndGrid(conf, element.second);

			anInfo->itsBaseGrid = move(g);
		}
		catch (...)
		{
			// do nothing
		}

		try
		{
			ParseLevels(anInfo, element.second);
		}
		catch (exception& e)
		{
			throw runtime_error(string("Error parsing level information: ") + e.what());
		}

		// Check local use_cache option

		bool delayedUseCache = conf->UseCache();

		try
		{
			string theUseCache = element.second.get<string>("use_cache");

			delayedUseCache = util::ParseBoolean(theUseCache);
		}
		catch (boost::property_tree::ptree_bad_path& e)
		{
			// Something was not found; do nothing
		}
		catch (exception& e)
		{
			throw runtime_error(string("Error parsing use_cache key: ") + e.what());
		}

		// Check local file_type option

		HPFileType delayedFileType = conf->itsOutputFileType;

		// Check local async options

		try
		{
			string async = element.second.get<string>("async");
			if (async == "true")
			{
				conf->AsyncExecution(true);
			}
		}
		catch (boost::property_tree::ptree_bad_path& e)
		{
			// Something was not found; do nothing
		}
		catch (exception& e)
		{
			throw runtime_error(string("Error parsing async key: ") + e.what());
		}

		try
		{
			string theFileType = boost::to_upper_copy(element.second.get<string>("file_type"));

			if (theFileType == "GRIB")
			{
				delayedFileType = kGRIB;
			}
			else if (theFileType == "GRIB1")
			{
				delayedFileType = kGRIB1;
			}
			else if (theFileType == "GRIB2")
			{
				delayedFileType = kGRIB2;
			}
			else if (theFileType == "FQD" || theFileType == "QUERYDATA")
			{
				delayedFileType = kQueryData;
			}
			else
			{
				throw runtime_error("Invalid option for 'file_type': " + theFileType);
			}
		}
		catch (boost::property_tree::ptree_bad_path& e)
		{
			// Something was not found; do nothing
		}
		catch (exception& e)
		{
			throw runtime_error(string("Error parsing meta information: ") + e.what());
		}

		// Check local file_write option

		HPFileWriteOption delayedFileWrite = conf->FileWriteOption();

		try
		{
			string theFileWriteOption = element.second.get<string>("file_write");

			if (theFileWriteOption == "database")
			{
				delayedFileWrite = kDatabase;
			}
			else if (theFileWriteOption == "single")
			{
				delayedFileWrite = kSingleFile;
			}
			else if (theFileWriteOption == "multiple")
			{
				delayedFileWrite = kMultipleFiles;
			}
			else if (theFileWriteOption == "cache only")
			{
				delayedFileWrite = kCacheOnly;
			}
			else
			{
				throw runtime_error("Invalid value for file_write: " + theFileWriteOption);
			}
		}
		catch (boost::property_tree::ptree_bad_path& e)
		{
			// Something was not found; do nothing
		}
		catch (exception& e)
		{
			throw runtime_error(string("Error parsing meta information: ") + e.what());
		}

		// Check local forecast_type option

		forecastTypes = ParseForecastTypes(element.second);

		if (!forecastTypes.empty())
		{
			anInfo->Set<forecast_type>(forecastTypes);
		}

		boost::property_tree::ptree& plugins = element.second.get_child("plugins");

		// Check local producer option

		try
		{
			ParseSourceProducer(conf, anInfo, element.second);
		}
		catch (...)
		{
		}

		try
		{
			ParseTargetProducer(conf, anInfo, element.second);
		}
		catch (...)
		{
		}

		if (plugins.empty())
		{
			throw runtime_error(ClassName() + ": plugin definitions not found");
		}

		for (boost::property_tree::ptree::value_type& plugin : plugins)
		{
			shared_ptr<plugin_configuration> pc = make_shared<plugin_configuration>(*conf);

			pc->UseCache(delayedUseCache);
			pc->itsOutputFileType = delayedFileType;
			pc->FileWriteOption(delayedFileWrite);

			if (plugin.second.empty())
			{
				throw runtime_error(ClassName() + ": plugin definition is empty");
			}

			for (boost::property_tree::ptree::value_type& kv : plugin.second)
			{
				string key = kv.first;
				string value;

				try
				{
					value = kv.second.get<string>("");
				}
				catch (...)
				{
					continue;
				}

				if (key == "name")
				{
					pc->Name(value);
				}
				else if (key == "param_list" || key == "options")
				{
					boost::property_tree::ptree params = plugin.second.get_child(key);

					if (params.empty())
					{
						throw runtime_error(ClassName() + ": param_list definition is empty");
					}

					for (boost::property_tree::ptree::value_type& param : params)
					{
						string name;
						std::vector<std::pair<std::string, std::string>> opts;

						for (boost::property_tree::ptree::value_type& paramOpt : param.second)
						{
							string paramOptName = paramOpt.first;
							string paramOptValue = paramOpt.second.get<string>("");

							if (paramOptName.empty())
							{
								throw runtime_error(ClassName() + ": param_list parameter option name is empty");
							}

							if (paramOptValue.empty())
							{
								throw runtime_error(ClassName() + ": param_list parameter option '" + paramOptName +
								                    "' value is empty");
							}

							if (paramOptName == "name" || paramOptName == "producer")
							{
								name = paramOptValue;
							}
							else
							{
								opts.push_back(std::make_pair(paramOptName, paramOptValue));
							}
						}
						pc->AddParameter(name, opts);
					}
				}
				else if (key == "async")
				{
					pc->AsyncExecution(util::ParseBoolean(value));
				}
				else
				{
					if (value.empty())
					{
						for (boost::property_tree::ptree::value_type& listval : kv.second)
						{
							// pc->AddOption(key, value);
							pc->AddOption(key, himan::util::Expand(listval.second.get<string>("")));
						}
					}
					else
					{
						pc->AddOption(key, himan::util::Expand(value));
					}
				}
			}

			if (pc->Name().empty())
			{
				throw runtime_error(ClassName() + ": plugin name not found from configuration");
			}

			pc->Info(make_shared<info>(*anInfo));  // We have to have a copy for all configs.
			                                       // Each plugin will later on create a data backend.

			ASSERT(pc.unique());

			pluginContainer.push_back(pc);
		}

	}  // END for

	return pluginContainer;
}

raw_time GetLatestOriginDateTime(const shared_ptr<configuration> conf, const string& latest)
{
	using namespace himan;

	auto strlist = himan::util::Split(latest, "-", false);

	unsigned int offset = 0;

	if (strlist.size() == 2)
	{
		// will throw if origintime is not in the form "latest-X", where X : integer >= 0
		offset = static_cast<unsigned>(stoi(strlist[1]));
	}

	HPDatabaseType dbtype = conf->DatabaseType();
	producer sourceProducer = conf->SourceProducer();

	raw_time latestOriginDateTime;

	auto r = GET_PLUGIN(radon);

	auto latestFromDatabase = r->RadonDB().GetLatestTime(static_cast<int>(sourceProducer.Id()), "", offset);

	if (!latestFromDatabase.empty())
	{
		return raw_time(latestFromDatabase, "%Y-%m-%d %H:%M:%S");
	}

	throw runtime_error("Latest time not found from " + HPDatabaseTypeToString.at(dbtype) + " for producer " +
	                    to_string(sourceProducer.Id()));
}

void json_parser::ParseTime(shared_ptr<configuration> conf, std::shared_ptr<info> anInfo,
                            const boost::property_tree::ptree& pt)
{
	/* Check origin time */
	const string mask = "%Y-%m-%d %H:%M:%S";

	std::vector<raw_time> originDateTimes;

	try
	{
		auto originDateTime = pt.get<string>("origintime");

		boost::algorithm::to_lower(originDateTime);

		if (originDateTime.find("latest") != string::npos)
		{
			if (conf->DatabaseType() == kNoDatabase)
			{
				throw std::invalid_argument("Unable to get latest time from database when no database mode is enabled");
			}
			originDateTimes.push_back(GetLatestOriginDateTime(conf, originDateTime));
		}
		else
		{
			originDateTimes.push_back(raw_time(originDateTime, mask));
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		try
		{
			auto datesList = himan::util::Split(pt.get<string>("origintimes"), ",", false);

			for (const auto& dateString : datesList)
			{
				originDateTimes.push_back(raw_time(dateString, mask));
			}
		}
		catch (boost::property_tree::ptree_bad_path& ee)
		{
			throw runtime_error("Origin datetime not found with keys 'origintime' or 'origindatetimes'");
		}
	}
	catch (exception& e)
	{
		throw runtime_error(ClassName() + ": " + string("Error parsing origin time information: ") + e.what());
	}

	ASSERT(!originDateTimes.empty());

	/* Check time steps */

	/*
	 * Three ways of providing information on steps:
	 * - hours
	 * - start_hour + stop_hour + step
	 * - start_minute + stop_minute + step
	 */

	try
	{
		string hours = pt.get<string>("hours");
		vector<string> timesStr = himan::util::Split(hours, ",", true);

		vector<int> times;

		for (size_t i = 0; i < timesStr.size(); i++)
		{
			times.push_back(stoi(timesStr[i]));
		}

		sort(times.begin(), times.end());

		vector<forecast_time> theTimes;

		// Create forecast_time with both times origintime, then adjust the validtime

		for (const auto& originDateTime : originDateTimes)
		{
			for (int hour : times)
			{
				forecast_time theTime(originDateTime, originDateTime);

				theTime.ValidDateTime().Adjust(kHourResolution, hour);

				theTimes.push_back(theTime);
			}
		}

		anInfo->Set<forecast_time>(theTimes);

		return;
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
	}
	catch (exception& e)
	{
		throw runtime_error(ClassName() + ": " + string("Error parsing time information from 'times': ") + e.what());
	}

	// hours was not specified
	// check if start/stop times are

	// First check step_unit which is deprecated and issue warning
	try
	{
		string stepUnit = pt.get<string>("step_unit");

		if (!stepUnit.empty())
		{
			itsLogger.Warning("Key 'step_unit' is deprecated");
		}
	}
	catch (exception& e)
	{
	}

	try
	{
		int start = pt.get<int>("start_hour");
		int stop = pt.get<int>("stop_hour");
		int step = pt.get<int>("step");

		if (step <= 0)
		{
			throw runtime_error("step size must be > 0");
		}

		conf->itsForecastStep = step;

		HPTimeResolution stepResolution = kHourResolution;

		vector<forecast_time> theTimes;

		for (const auto& originDateTime : originDateTimes)
		{
			int curtime = start;

			do
			{
				forecast_time theTime(originDateTime, originDateTime);

				theTime.ValidDateTime().Adjust(stepResolution, curtime);

				theTime.StepResolution(stepResolution);

				theTimes.push_back(theTime);

				curtime += step;

			} while (curtime <= stop);
		}

		anInfo->Set<forecast_time>(theTimes);

		return;
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
	}
	catch (exception& e)
	{
		throw runtime_error(ClassName() + ": " + string("Error parsing time information from 'start_hour': ") +
		                    e.what());
	}

	try
	{
		// try start_minute/stop_minute

		int start = pt.get<int>("start_minute");
		int stop = pt.get<int>("stop_minute");
		int step = pt.get<int>("step");

		conf->itsForecastStep = step;

		HPTimeResolution stepResolution = kMinuteResolution;

		int curtime = start;

		vector<forecast_time> theTimes;

		for (const auto& originDateTime : originDateTimes)
		{
			do
			{
				forecast_time theTime(originDateTime, originDateTime);

				theTime.ValidDateTime().Adjust(stepResolution, curtime);

				theTime.StepResolution(stepResolution);

				theTimes.push_back(theTime);

				curtime += step;

			} while (curtime <= stop);
		}

		anInfo->Set<forecast_time>(theTimes);
	}
	catch (exception& e)
	{
		throw runtime_error(ClassName() + ": " + string("Error parsing time information: ") + e.what());
	}
}

unique_ptr<grid> ParseAreaAndGridFromDatabase(configuration& conf, const boost::property_tree::ptree& pt)
{
	using himan::kBottomLeft;
	using himan::kTopLeft;

	unique_ptr<grid> g;

	logger log("json_parser");

	try
	{
		string geom = pt.get<string>("target_geom_name");

		conf.TargetGeomName(geom);

		g = util::GridFromDatabase(geom);
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		log.Fatal(string("Error parsing area information: ") + e.what());
		himan::Abort();
	}

	return g;
}

unique_ptr<grid> ParseAreaAndGridFromPoints(configuration& conf, const boost::property_tree::ptree& pt)
{
	unique_ptr<grid> g;

	// check points

	try
	{
		vector<string> stations = himan::util::Split(pt.get<string>("points"), ",", false);

		g = unique_ptr<point_list>(new point_list());

		vector<station> theStations;

		int i = 1;

		for (const string& line : stations)
		{
			vector<string> point = himan::util::Split(line, " ", false);

			if (point.size() != 2)
			{
				cout << "Error::json_parser Line " + line + " is invalid" << endl;
				continue;
			}

			string lon = point[0];
			string lat = point[1];

			boost::algorithm::trim(lon);
			boost::trim(lat);

			theStations.push_back(station(i, "point_" + to_string(i), stod(lon), stod(lat)));

			i++;
		}

		if (theStations.size())
		{
			dynamic_cast<point_list*>(g.get())->Stations(theStations);
			return g;
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (const exception& e)
	{
		throw runtime_error(string("Fatal::json_parser Error parsing points: ") + e.what());
	}

	// Check stations

	try
	{
		vector<string> stations = himan::util::Split(pt.get<string>("stations"), ",", false);

		g = unique_ptr<point_list>(new point_list);

		vector<station> theStations;

		auto r = GET_PLUGIN(radon);

		for (const string& str : stations)
		{
			unsigned long fmisid;

			try
			{
				fmisid = static_cast<unsigned long>(stol(str));
			}
			catch (const invalid_argument& e)
			{
				cout << "Error::json_parser Invalid fmisid: " << str << endl;
				continue;
			}

			auto stationinfo = r->RadonDB().GetStationDefinition(kFmiSIDNetwork, fmisid);

			if (stationinfo.empty())
			{
				cout << "Error::json_parser Station " << str << " not found from database" << endl;
				continue;
			}

			theStations.push_back(station(static_cast<int>(fmisid), stationinfo["station_name"],
			                              stod(stationinfo["longitude"]), stod(stationinfo["latitude"])));
		}

		dynamic_cast<point_list*>(g.get())->Stations(theStations);
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (const exception& e)
	{
		throw runtime_error(string("Fatal::json_parser Error parsing stations: ") + e.what());
	}

	if (g && dynamic_cast<point_list*>(g.get())->Stations().empty())
	{
		throw runtime_error("Fatal::json_parser No valid points or stations found");
	}

	return g;
}

unique_ptr<grid> json_parser::ParseAreaAndGrid(shared_ptr<configuration> conf, const boost::property_tree::ptree& pt)
{
	/*
	 * Parse area and grid from different possible options.
	 * Order or parsing:
	 *
	 * 1. 'source_geom_name': this is used in fetching data, it's not used to create an area instance
	 * 2. radon style geom name: 'target_geom_name'
	 * 3. irregular grid: 'points' and 'stations'
	 * 4. bounding box: 'bbox'
	 * 5. manual definition:
	 * -> 'projection',
	 * -> 'bottom_left_longitude', 'bottom_left_latitude',
	 * -> 'top_right_longitude', 'top_right_latitude'
	 * -> 'orientation'
	 * -> 'south_pole_longitude', 'south_pole_latitude'
	 * -> 'ni', 'nj'
	 * -> 'scanning_mode'
	 *
	 */

	// 1. Check for source geom name

	try
	{
		conf->SourceGeomNames(himan::util::Split(pt.get<string>("source_geom_name"), ",", false));
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing area information: ") + e.what());
	}

	// 2. radon-style geom_name

	auto g = ParseAreaAndGridFromDatabase(*conf, pt);

	if (g)
	{
		return g;
	}

	// 3. Points

	auto ig = ParseAreaAndGridFromPoints(*conf, pt);

	if (ig)
	{
		// Disable cuda interpolation (too inefficienct for single points)
		itsLogger.Trace("Disabling cuda interpolation for single point data");
		conf->UseCudaForInterpolation(false);
		return ig;
	}

	// 4. Target geometry is still not set, check for bbox

	unique_ptr<grid> rg;

	try
	{
		vector<string> coordinates = himan::util::Split(pt.get<string>("bbox"), ",", false);

		rg = unique_ptr<latitude_longitude_grid>(new latitude_longitude_grid);

		dynamic_cast<latitude_longitude_grid*>(rg.get())->BottomLeft(point(stod(coordinates[0]), stod(coordinates[1])));
		dynamic_cast<latitude_longitude_grid*>(rg.get())->TopRight(point(stod(coordinates[2]), stod(coordinates[3])));

		dynamic_cast<latitude_longitude_grid*>(rg.get())->Ni(pt.get<size_t>("ni"));
		dynamic_cast<latitude_longitude_grid*>(rg.get())->Nj(pt.get<size_t>("nj"));

		dynamic_cast<latitude_longitude_grid*>(rg.get())->ScanningMode(
		    HPScanningModeFromString.at(pt.get<string>("scanning_mode")));

		return rg;
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing bbox: ") + e.what());
	}

	// 5. Check for manual definition of area

	try
	{
		HPScanningMode mode = HPScanningModeFromString.at(pt.get<string>("scanning_mode"));

		if (mode != kTopLeft && mode != kBottomLeft)
		{
			throw runtime_error(ClassName() + ": scanning mode " + HPScanningModeToString.at(mode) +
			                    " not supported (ever)");
		}

		string projection = pt.get<string>("projection");

		if (projection == "latlon")
		{
			rg = unique_ptr<latitude_longitude_grid>(new latitude_longitude_grid);
			dynamic_cast<latitude_longitude_grid*>(rg.get())->ScanningMode(mode);

			dynamic_cast<latitude_longitude_grid*>(rg.get())->BottomLeft(
			    point(pt.get<double>("bottom_left_longitude"), pt.get<double>("bottom_left_latitude")));
			dynamic_cast<latitude_longitude_grid*>(rg.get())->TopRight(
			    point(pt.get<double>("top_right_longitude"), pt.get<double>("top_right_latitude")));
			dynamic_cast<latitude_longitude_grid*>(rg.get())->Ni(pt.get<size_t>("ni"));
			dynamic_cast<latitude_longitude_grid*>(rg.get())->Nj(pt.get<size_t>("nj"));
		}
		else if (projection == "rotated_latlon")
		{
			rg = unique_ptr<rotated_latitude_longitude_grid>(new rotated_latitude_longitude_grid);
			dynamic_cast<rotated_latitude_longitude_grid*>(rg.get())->ScanningMode(mode);

			dynamic_cast<rotated_latitude_longitude_grid*>(rg.get())->BottomLeft(
			    point(pt.get<double>("bottom_left_longitude"), pt.get<double>("bottom_left_latitude")));
			dynamic_cast<rotated_latitude_longitude_grid*>(rg.get())->TopRight(
			    point(pt.get<double>("top_right_longitude"), pt.get<double>("top_right_latitude")));
			dynamic_cast<rotated_latitude_longitude_grid*>(rg.get())->SouthPole(
			    point(pt.get<double>("south_pole_longitude"), pt.get<double>("south_pole_latitude")));
			dynamic_cast<rotated_latitude_longitude_grid*>(rg.get())->Ni(pt.get<size_t>("ni"));
			dynamic_cast<rotated_latitude_longitude_grid*>(rg.get())->Nj(pt.get<size_t>("nj"));
		}
		else if (projection == "stereographic")
		{
			rg = unique_ptr<stereographic_grid>(new stereographic_grid);
			dynamic_cast<stereographic_grid*>(rg.get())->ScanningMode(mode);

			dynamic_cast<stereographic_grid*>(rg.get())->FirstPoint(
			    point(pt.get<double>("first_point_longitude"), pt.get<double>("first_point_latitude")));
			dynamic_cast<stereographic_grid*>(rg.get())->Di(pt.get<double>("di"));
			dynamic_cast<stereographic_grid*>(rg.get())->Dj(pt.get<double>("dj"));
			dynamic_cast<stereographic_grid*>(rg.get())->Orientation(pt.get<double>("orientation"));
			dynamic_cast<stereographic_grid*>(rg.get())->Ni(pt.get<size_t>("ni"));
			dynamic_cast<stereographic_grid*>(rg.get())->Nj(pt.get<size_t>("nj"));
		}
		else
		{
			throw runtime_error(ClassName() + ": Unknown type: " + projection);
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		throw runtime_error(e.what());
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing area: ") + e.what());
	}

	rg->EarthShape(earth_shape<double>(6371220.));

	return rg;
}

void ParseSourceProducer(shared_ptr<configuration> conf, shared_ptr<info> anInfo, const boost::property_tree::ptree& pt)
{
	std::vector<producer> sourceProducers;
	vector<string> sourceProducersStr = himan::util::Split(pt.get<string>("source_producer"), ",", false);

	const HPDatabaseType dbtype = conf->DatabaseType();

	if (dbtype == kRadon)
	{
		auto r = GET_PLUGIN(radon);

		for (const auto& prodstr : sourceProducersStr)
		{
			long pid = stol(prodstr);

			producer prod(pid);

			map<string, string> prodInfo = r->RadonDB().GetProducerDefinition(static_cast<unsigned long>(pid));

			if (!prodInfo.empty())
			{
				prod.Name(prodInfo["ref_prod"]);

				if (!prodInfo["ident_id"].empty())
				{
					prod.Centre(stol(prodInfo["ident_id"]));
					prod.Process(stol(prodInfo["model_id"]));
				}

				prod.Class(static_cast<HPProducerClass>(stoi(prodInfo["producer_class"])));

				sourceProducers.push_back(prod);
			}
			else
			{
				itsLogger.Warning("Failed to find source producer from Radon: " + prodstr);
			}
		}
	}
	else if (dbtype != kNoDatabase && sourceProducers.size() == 0)
	{
		itsLogger.Fatal("Source producer information was not found from database");
		himan::Abort();
	}
	else if (dbtype == kNoDatabase)
	{
		for (const auto& prodstr : sourceProducersStr)
		{
			sourceProducers.push_back(producer(stoi(prodstr)));
		}
	}

	conf->SourceProducers(sourceProducers);
}

void ParseTargetProducer(shared_ptr<configuration> conf, shared_ptr<info> anInfo, const boost::property_tree::ptree& pt)
{
	const HPDatabaseType dbtype = conf->DatabaseType();

	/*
	 * Target producer is also set to target info; source infos (and producers) are created
	 * as data is fetched from files.
	 */

	long pid = stol(pt.get<string>("target_producer"));
	producer prod(pid);

	auto r = GET_PLUGIN(radon);
	auto prodInfo = r->RadonDB().GetProducerDefinition(static_cast<unsigned long>(pid));

	if (!prodInfo.empty())
	{
		if (prodInfo["ident_id"].empty() || prodInfo["model_id"].empty())
		{
			itsLogger.Warning("Centre or ident information not found for producer " + prodInfo["ref_prod"]);
		}
		else
		{
			prod.Centre(stol(prodInfo["ident_id"]));
			prod.Process(stol(prodInfo["model_id"]));
		}

		prod.Name(prodInfo["ref_prod"]);

		if (prodInfo["producer_class"].empty())
		{
			prod.Class(kGridClass);
		}
		else
		{
			prod.Class(static_cast<HPProducerClass>(stoi(prodInfo["producer_class"])));
		}
	}
	else if (dbtype != kNoDatabase)
	{
		itsLogger.Warning("Unknown target producer: " + pt.get<string>("target_producer"));
	}

	conf->TargetProducer(prod);
	anInfo->Producer(prod);
}

void ParseLevels(shared_ptr<info> anInfo, const boost::property_tree::ptree& pt)
{
	try
	{
		string levelTypeStr = pt.get<string>("leveltype");
		string levelValuesStr = pt.get<string>("levels");

		vector<level> levels = LevelsFromString(levelTypeStr, levelValuesStr);

		anInfo->Set<level>(levels);
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
		throw runtime_error(e.what());
		// Something was not found; do nothing
	}
	catch (exception& e)
	{
		throw runtime_error(e.what());
	}
}

vector<level> LevelsFromString(const string& levelType, const string& levelValues)
{
	HPLevelType theLevelType = HPStringToLevelType.at(boost::to_lower_copy(levelType));
	vector<level> levels;

	if (theLevelType == kHeightLayer || theLevelType == kGroundDepth || theLevelType == kPressureDelta)
	{
		const vector<string> levelsStr = himan::util::Split(levelValues, ",", false);
		for (size_t i = 0; i < levelsStr.size(); i++)
		{
			const vector<string> levelIntervals = himan::util::Split(levelsStr[i], "_", false);

			if (levelIntervals.size() != 2)
			{
				throw runtime_error(
				    "height_layer, ground_depth and pressure delta requires two level values per definition (lx1_ly1, "
				    "lx2_ly2, ..., "
				    "lxN_lyN)");
			}

			levels.push_back(level(theLevelType, stof(levelIntervals[0]), stof(levelIntervals[1])));
		}
	}
	else
	{
		const vector<string> levelsStr = himan::util::Split(levelValues, ",", true);
		for (size_t i = 0; i < levelsStr.size(); i++)
		{
			levels.push_back(level(theLevelType, stof(levelsStr[i]), levelType));
		}
	}

	ASSERT(!levels.empty());

	return levels;
}

vector<forecast_type> ParseForecastTypes(const boost::property_tree::ptree& pt)
{
	vector<forecast_type> forecastTypes;

	try
	{
		vector<string> types = himan::util::Split(pt.get<string>("forecast_type"), ",", false);

		for (string& type : types)
		{
			boost::algorithm::to_lower(type);
			HPForecastType forecastType;

			if (type.find("pf") != string::npos)
			{
				forecastType = kEpsPerturbation;
				string list = "";
				for (size_t i = 2; i < type.size(); i++)
					list += type[i];

				vector<string> range = himan::util::Split(list, "-", false);

				if (range.size() == 1)
				{
					forecastTypes.push_back(forecast_type(forecastType, stod(range[0])));
				}
				else
				{
					ASSERT(range.size() == 2);

					int start = stoi(range[0]);
					int stop = stoi(range[1]);

					while (start <= stop)
					{
						forecastTypes.push_back(forecast_type(forecastType, start));
						start++;
					}
				}
			}
			else
			{
				if (type == "cf")
				{
					forecastTypes.push_back(forecast_type(kEpsControl, 0));
				}
				else if (type == "det" || type == "deterministic")
				{
					forecastTypes.push_back(forecast_type(kDeterministic));
				}
				else if (type == "an" || type == "analysis")
				{
					forecastTypes.push_back(forecast_type(kAnalysis));
				}
				else
				{
					throw runtime_error("Invalid forecast_type: " + type);
				}
			}
		}
	}
	catch (boost::property_tree::ptree_bad_path& e)
	{
	}
	catch (const boost::exception& e)
	{
		throw runtime_error(string("Invalid forecast_type value"));
	}
	catch (exception& e)
	{
		throw runtime_error(string("Error parsing key forecast_type: ") + e.what());
	}

	return forecastTypes;
}
