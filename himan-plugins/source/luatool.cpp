#include "luatool.h"
#include <boost/filesystem.hpp>
#include "logger_factory.h"
#include "forecast_time.h"
#include "util.h"
#include "metutil.h"
#include "regular_grid.h"
#include "plugin_factory.h"
#include <boost/foreach.hpp>

#define HIMAN_AUXILIARY_INCLUDE

#include "hitool.h"

#undef HIMAN_AUXILIARY_INCLUDE

extern "C"
{
#include <lualib.h>
}

namespace luabind
{
namespace detail
{
namespace has_get_pointer_
{
	template<class T>
	T * get_pointer(std::shared_ptr<T> const& p) { return p.get(); }
}
}
}

#include <luabind/luabind.hpp>
#include <luabind/object.hpp>
#include <luabind/return_reference_to_policy.hpp>
#include <luabind/iterator_policy.hpp>
#include <luabind/adopt_policy.hpp>
#include <luabind/operator.hpp>

#include <boost/get_pointer.hpp>

using namespace himan;
using namespace himan::plugin;
using namespace luabind;

#define LUA_MEMFN(r, t,  m, ...) static_cast<r(t::*)(__VA_ARGS__)>(&t::m)
#define LUA_CMEMFN(r, t, m, ...) static_cast<r(t::*)(__VA_ARGS__) const>(&t::m)

void BindLib(lua_State* L);
void BindPlugins(lua_State* L);
void BindEnum(lua_State* L);
int BindErrorHandler(lua_State* L);

luatool::luatool() : L(0)
{
	itsClearTextFormula = "<interpreted>";
	itsLogger = logger_factory::Instance()->GetLog("luatool");
}

luatool::~luatool()
{
	if (L)
	{
		lua_close(L);
	}
}

void luatool::Process(std::shared_ptr<const plugin_configuration> conf)
{
	Init(conf);

	SetParams({param("DUMMY")});

	Start();

}


void luatool::Calculate(std::shared_ptr<info> myTargetInfo, unsigned short threadIndex)
{

	InitLua(myTargetInfo);

	auto myThreadedLogger = logger_factory::Instance()->GetLog("luatoolThread #" + boost::lexical_cast<std::string> (threadIndex));

	globals(L)["logger"] = myThreadedLogger.get();

	BOOST_FOREACH(const std::string& luaFile, itsConfiguration->GetValueList("luafile"))
	{
		if (luaFile.empty())
		{
			continue;
		}

		itsLogger->Info("Starting script " + luaFile);

		ReadFile(luaFile);
	}
}

void luatool::InitLua(std::shared_ptr<info> myTargetInfo)
{

	L = luaL_newstate();
	luaL_openlibs(L);

	assert(L);

	open(L);

	set_pcall_callback(&BindErrorHandler);

    BindEnum(L);
    BindLib(L);
    BindPlugins(L);

	itsLogger->Trace("luabind finished");

	// Set some variable that are needed in luatool calculations
	// but are too hard or complicated to create in the lua side

	globals(L)["luatool"] = boost::ref(this);
	globals(L)["result"] = myTargetInfo;
	globals(L)["current_time"] = myTargetInfo->Time();
	globals(L)["current_level"] = myTargetInfo->Level();

	auto h = std::dynamic_pointer_cast<hitool> (plugin_factory::Instance()->Plugin("hitool"));

	h->Configuration(itsConfiguration);
	h->Time(myTargetInfo->Time());

	globals(L)["hitool"] = h;

	// Define a nice iterator for multiple infos

	const char* nextvalue = R"(
function nextvalue(...)
	local arg = {...}
	local tmp = {...}

	if #tmp == 0 then
		return function() end, nil, nil
	end

	local function _nextvalue(arg, i)
		i = i+1

		local gridsize = arg[1]:GetGrid():Size()

		if i > gridsize then return nil end

		for j=1,#arg do
			local f = arg[j]
			local val = f:GetIndexValue(i)
			if val == nil then return val end

			tmp[j] = val
		end

		return i, unpack(tmp)
	end

	return _nextvalue, arg, 0
end
)";

	luaL_dostring(L, nextvalue);
}

bool luatool::ReadFile(const std::string& luaFile)
{
	if (!boost::filesystem::exists(luaFile))
	{
		std::cerr << "Error: script " << luaFile << " does not exist\n";
		return false;
	}
	
	try
	{
		if (luaL_dofile(L, luaFile.c_str()))
		{
			itsLogger->Error(lua_tostring(L, -1));
			return false;
		}
	}
	catch (const error& e)
	{
		return false;
	}
	catch (const std::exception& e)
	{
		itsLogger->Error(e.what());
		return false;
	}

	return true;

}

int BindErrorHandler(lua_State* L)
{
    // log the error message
    luabind::object msg(luabind::from_stack( L, -1 ));
    std::ostringstream str;
    str << "lua> run-time error: " << msg;
    std::cout << str.str() << std::endl;

    // log the callstack
    std::string traceback = luabind::call_function<std::string>( luabind::globals(L)["debug"]["traceback"] );
    traceback = std::string("lua> ") + traceback;
    std::cout << traceback.c_str() << std::endl;

    // return unmodified error object
    return 1;
}

void BindEnum(lua_State* L)
{

	module(L)
	[
		class_<HPLevelType>("HPLevelType")
			.enum_("constants")
		[
			value("kUnknownLevel", kUnknownLevel),
			value("kGround", kHeight),
			value("kTopOfAtmosphere", kHeight),
			value("kPressure", kHeight),
			value("kMeanSea", kHeight),
			value("kAltitude", kHeight),
			value("kHeight", kHeight),
			value("kHybrid", kHybrid),
			value("kGndLayer", kGndLayer),
			value("kDepth", kDepth),
			value("kEntireAtmosphere", kEntireAtmosphere),
			value("kEntireOcean", kEntireOcean)
		]
		,
		class_<HPTimeResolution>("HPTimeResolution")
			.enum_("constants")
		[
			value("kUnknownTimeResolution", kUnknownTimeResolution),
			value("kHourResolution", kHourResolution),
			value("kMinuteResolution", kMinuteResolution)
		]
		,
		class_<HPFileType>("HPFileType")
			.enum_("constants")
		[
			value("kUnknownFile", kUnknownFile),
			value("kGRIB1", kGRIB1),
			value("kGRIB2", kGRIB2),
			value("kGRIB", kGRIB),
			value("kQueryData", kQueryData),
			value("kNetCDF", kNetCDF)
		]
		,
		class_<HPProjectionType>("HPProjectionType")
			.enum_("constants")
		[
			value("kUnknownProjection", kUnknownProjection),
			value("kLatLonProjection", kLatLonProjection),
			value("kRotatedLatLonProjection", kRotatedLatLonProjection),
			value("kStereographicProjection", kStereographicProjection)
		]
		,
		class_<HPScanningMode>("HPScanningMode")
			.enum_("constants")
		[
			value("kUnknownScanningMode", kUnknownScanningMode),
			value("kTopLeft", kTopLeft),
			value("kTopRight", kTopRight),
			value("kBottomLeft", kBottomLeft),
			value("kBottomRight", kBottomRight)
		]
		,
		class_<HPAggregationType>("HPAggregationType")
			.enum_("constants")
		[
			value("kUnknownAggregationType", kUnknownAggregationType),
			value("kAverage", kAverage),
			value("kAccumulation", kAccumulation),
			value("kMaximum", kMaximum),
			value("kMinimum", kMinimum),
			value("kDifference", kDifference)
		],
		class_<HPModifierType>("HPModifierType")
			.enum_("constants")
		[
			value("kUnknownModifierType", kUnknownModifierType),
			value("kAverageModifier", kAverageModifier),
			value("kAccumulationModifier", kAccumulationModifier),
			value("kMaximumModifier", kMaximumModifier),
			value("kMinimumModifier", kMinimumModifier),
			value("kDifferenceModifier", kDifferenceModifier),
			value("kMaximumMinimumModifier", kMaximumMinimumModifier),
			value("kCountModifier", kCountModifier),
			value("kFindHeightModifier", kFindHeightModifier),
			value("kFindValueModifier", kFindValueModifier),
			value("kIntegralModifier", kIntegralModifier),
			value("kPlusMinusAreaModifier", kPlusMinusAreaModifier)
		],
		class_<HPGridType>("HPGridType")
			.enum_("constants")
		[
			value("kUnknownGridType", kUnknownGridType),
			value("kRegularGrid", kRegularGrid),
			value("kIrregularGrid", kIrregularGrid)
		]
	];
}

namespace info_wrapper
{
// These are convenience functions for accessing info class contents

const std::vector<double>& Locations(std::shared_ptr<info> anInfo)
{
	return anInfo->Grid()->Data().Values();
}

bool SetValue(std::shared_ptr<info> anInfo, int index, double value)
{
	return anInfo->Grid()->Value(--index, value);
}

double GetValue(std::shared_ptr<info> anInfo, int index)
{
	return anInfo->Grid()->Value(--index);
}

bool SetData(std::shared_ptr<info> anInfo, const std::vector<double>& values)
{
	return anInfo->Grid()->Data().Set(values);
}

} // namespace info_wrapper

namespace hitool_wrapper
{
// The following functions are all wrappers for hitool:
// we cannot specify hitool functions directly in the lua binding
// because that leads to undefined symbols in plugins and that
// forces us to link luatool with hitool which is not nice!

std::vector<double> VerticalMaximumMultiParamGrid(std::shared_ptr<hitool> h, const params& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->VerticalMaximum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalMaximumMultiParam(std::shared_ptr<hitool> h, const params& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->VerticalMaximum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalMaximumGrid(std::shared_ptr<hitool> h, const param& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->VerticalMaximum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalMaximum(std::shared_ptr<hitool> h, const param& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->VerticalMaximum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalMinimumMultiParamGrid(std::shared_ptr<hitool> h, const params& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->VerticalMinimum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalMinimumMultiParam(std::shared_ptr<hitool> h, const params& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->VerticalMinimum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalMinimumGrid(std::shared_ptr<hitool> h, const param& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->VerticalMinimum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalMinimum(std::shared_ptr<hitool> h, const param& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->VerticalMinimum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalSumMultiParamGrid(std::shared_ptr<hitool> h, const params& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->VerticalSum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalSumMultiParam(std::shared_ptr<hitool> h, const params& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->VerticalSum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalSumGrid(std::shared_ptr<hitool> h, const param& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->VerticalSum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalSum(std::shared_ptr<hitool> h, const param& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->VerticalSum(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalAverageMultiParamGrid(std::shared_ptr<hitool> h, const params& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->VerticalAverage(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalAverageMultiParam(std::shared_ptr<hitool> h, const params& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->VerticalAverage(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalAverageGrid(std::shared_ptr<hitool> h, const param& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->VerticalAverage(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalAverage(std::shared_ptr<hitool> h, const param& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->VerticalAverage(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalCountMultiParamGrid(std::shared_ptr<hitool> h, const params& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue, const std::vector<double>& findValue)
{
	return h->VerticalCount(theParams, firstLevelValue, lastLevelValue, findValue);
}

std::vector<double> VerticalCountMultiParam(std::shared_ptr<hitool> h, const params& theParams, double firstLevelValue, double lastLevelValue, double findValue)
{
	return h->VerticalCount(theParams, firstLevelValue, lastLevelValue, findValue);
}

std::vector<double> VerticalCountGrid(std::shared_ptr<hitool> h, const param& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue, const std::vector<double>& findValue)
{
	return h->VerticalCount(theParams, firstLevelValue, lastLevelValue, findValue);
}

std::vector<double> VerticalCount(std::shared_ptr<hitool> h, const param& theParams, double firstLevelValue, double lastLevelValue, double findValue)
{
	return h->VerticalCount(theParams, firstLevelValue, lastLevelValue, findValue);
}

std::vector<double> VerticalHeightMultiParamGrid(std::shared_ptr<hitool> h, const params& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue, const std::vector<double>& findValue, size_t findNth)
{
	return h->VerticalHeight(theParams, firstLevelValue, lastLevelValue, findValue, findNth);
}

std::vector<double> VerticalHeightMultiParam(std::shared_ptr<hitool> h, const params& theParams, double firstLevelValue, double lastLevelValue, double findValue, size_t findNth)
{
	return h->VerticalHeight(theParams, firstLevelValue, lastLevelValue, findValue, findNth);
}

std::vector<double> VerticalHeightGrid(std::shared_ptr<hitool> h, const param& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue, const std::vector<double>& findValue, size_t findNth)
{
	return h->VerticalHeight(theParams, firstLevelValue, lastLevelValue, findValue, findNth);
}

std::vector<double> VerticalHeight(std::shared_ptr<hitool> h, const param& theParams, double firstLevelValue, double lastLevelValue, double findValue, size_t findNth)
{
	return h->VerticalHeight(theParams, firstLevelValue, lastLevelValue, findValue, findNth);
}

std::vector<double> VerticalValueMultiParamGrid(std::shared_ptr<hitool> h, const params& theParams, const std::vector<double>& findValue)
{
	return h->VerticalValue(theParams, findValue);
}

std::vector<double> VerticalValueMultiParam(std::shared_ptr<hitool> h, const params& theParams, double findValue)
{
	return h->VerticalValue(theParams, findValue);
}

std::vector<double> VerticalValueGrid(std::shared_ptr<hitool> h, const param& theParams, const std::vector<double>& findValue)
{
	return h->VerticalValue(theParams, findValue);
}

std::vector<double> VerticalValue(std::shared_ptr<hitool> h, const param& theParams, double findValue)
{
	return h->VerticalValue(theParams, findValue);
}

std::vector<double> VerticalPlusMinusAreaMultiParamGrid(std::shared_ptr<hitool> h, const params& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->PlusMinusArea(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalPlusMinusAreaMultiParam(std::shared_ptr<hitool> h, const params& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->PlusMinusArea(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalPlusMinusAreaGrid(std::shared_ptr<hitool> h, const param& theParams, const std::vector<double>& firstLevelValue, const std::vector<double>& lastLevelValue)
{
	return h->PlusMinusArea(theParams, firstLevelValue, lastLevelValue);
}

std::vector<double> VerticalPlusMinusArea(std::shared_ptr<hitool> h, const param& theParams, double firstLevelValue, double lastLevelValue)
{
	return h->PlusMinusArea(theParams, firstLevelValue, lastLevelValue);
}

void Time(std::shared_ptr<hitool> h, const forecast_time& theTime)
{
	h->Time(theTime);
}
} // namespace hitool_wrapper

void BindLib(lua_State* L)
{
	module(L)
	[
		class_<himan::info, std::shared_ptr<himan::info>>("info")
			.def(constructor<>())
			.def("ClassName", &info::ClassName)
			.def("First", &info::First)
			.def("ResetLocation", &info::ResetLocation)
			.def("FirstLocation", &info::FirstLocation)
			.def("NextLocation", &info::NextLocation)
			.def("ResetParam", &info::ResetParam)
			.def("FirstParam", &info::FirstParam)
			.def("NextParam", &info::NextParam)
			.def("ResetLevel", &info::ResetLevel)
			.def("FirstLevel", &info::FirstLevel)
			.def("NextLevel", &info::NextLevel)
			.def("ResetTime", &info::ResetTime)
			.def("FirstTime", &info::FirstTime)
			.def("NextTime", &info::NextTime)
			.def("GetValue", LUA_CMEMFN(double, info, Value, void))
			.def("SetValue", LUA_MEMFN(bool, info, Value, double))
			.def("SizeLocations", LUA_CMEMFN(size_t, info, SizeLocations, void))
			.def("SizeTimes", LUA_CMEMFN(size_t, info, SizeTimes, void))
			.def("SizeParams", LUA_CMEMFN(size_t, info, SizeParams, void))
			.def("SizeLevels", LUA_CMEMFN(size_t, info, SizeLevels, void))
			.def("LocationIndex", LUA_CMEMFN(size_t, info, LocationIndex, void))
			.def("GetTimeIndex", LUA_CMEMFN(size_t, info, TimeIndex, void))
			.def("GetParamIndex", LUA_CMEMFN(size_t, info, ParamIndex, void))
			.def("GetLevelIndex", LUA_CMEMFN(size_t, info, LevelIndex, void))
			.def("GetLevel", LUA_CMEMFN(level, info, Level, void))
			.def("GetTime", LUA_CMEMFN(forecast_time, info, Time, void))
			.def("GetGrid", LUA_CMEMFN(grid*, info, Grid, void))
			.def("SetParam", LUA_MEMFN(void, info, SetParam, const param&))
			.def("SetTime", LUA_MEMFN(void, info, SetTime, const forecast_time&))
			.def("SetLevel", LUA_MEMFN(void, info, SetLevel, const level&))
			.def("LatLon", &info::LatLon)
			// These are local functions to luatool
			.def("Locations", &info_wrapper::Locations, return_stl_iterator)
			.def("SetIndexValue", &info_wrapper::SetValue)
			.def("GetIndexValue", &info_wrapper::GetValue)
			.def("SetData", &info_wrapper::SetData)
		,
		class_<grid, std::shared_ptr<grid>>("grid")
		,
		class_<regular_grid, grid, std::shared_ptr<regular_grid>>("regular_grid")
			.def(constructor<>())
			.def("ClassName", &grid::ClassName)
			.def("Size", &grid::Size)
			.def("GetData", LUA_MEMFN(matrix<double>&, grid, Data, void))
			.def("GetNi", LUA_CMEMFN(size_t, regular_grid, Ni, void))
			.def("GetNj", LUA_CMEMFN(size_t, regular_grid, Nj, void))
			.def("GetDi", LUA_CMEMFN(double, regular_grid, Di, void))
			.def("GetDj", LUA_CMEMFN(double, regular_grid, Dj, void))
			.def("GetScanningMode", LUA_CMEMFN(HPScanningMode, regular_grid, ScanningMode, void))
			.def("GetProjection", LUA_CMEMFN(HPProjectionType, regular_grid, Projection, void))
			.def("GetAB", LUA_CMEMFN(std::vector<double>, regular_grid, AB, void))
			.def("SetAB", LUA_MEMFN(void, regular_grid, AB, const std::vector<double>&))
			.def("GetBottomLeft", LUA_CMEMFN(point, regular_grid, BottomLeft, void))
			.def("SetBottomLeft", LUA_MEMFN(void, regular_grid, BottomLeft, const point&))
			.def("GetTopRight", LUA_CMEMFN(point, regular_grid, TopRight, void))
			.def("SetTopRight", LUA_MEMFN(void, regular_grid, BottomLeft, const point&))
			.def("GetFirstGridPoint", LUA_CMEMFN(point, regular_grid, FirstGridPoint, void))
			.def("GetLastGridPoint", LUA_CMEMFN(point, regular_grid, LastGridPoint, void))
			.def("Swap", &regular_grid::Swap)
			.def("LatLon", &regular_grid::LatLon)
		,
		class_<matrix<double>>("data")
			.def(constructor<size_t,size_t,size_t>())
			.def("Size", &matrix<double>::Size)
			.def("ClassName", &matrix<double>::ClassName)
			.def("GetSizeX", LUA_CMEMFN(size_t, matrix<double>, SizeX, void))
			.def("GetSizeY", LUA_CMEMFN(size_t, matrix<double>, SizeY, void))
			.def("GetSizeZ", LUA_CMEMFN(size_t, matrix<double>, SizeZ, void))
			.def("SetSizeX", LUA_MEMFN(void, matrix<double>, SizeX, size_t))
			.def("SetSizeY", LUA_MEMFN(void, matrix<double>, SizeY, size_t))
			.def("SetSizeZ", LUA_MEMFN(void, matrix<double>, SizeZ, size_t))
			.def("Resize", LUA_MEMFN(void, matrix<double>, Resize, size_t, size_t, size_t))
			.def("GetValue", LUA_CMEMFN(double, matrix<double>, At, size_t))
			.def("GetValues", &matrix<double>::Values)
			.def("SetValue", LUA_MEMFN(bool, matrix<double>, Set, size_t, double))
			.def("SetValues", LUA_MEMFN(bool, matrix<double>, Set, const std::vector<double>&))
			.def("Fill", &matrix<double>::Fill)
			.def("GetMissingValue", LUA_CMEMFN(double, matrix<double>, MissingValue, void))
			.def("SetMissingValue", LUA_MEMFN(void, matrix<double>, MissingValue, double))
			.def("Clear", &matrix<double>::Clear)
			.def("IsMissing", LUA_CMEMFN(bool, matrix<double>, IsMissing, size_t))
			.def("MissingCount", &matrix<double>::MissingCount)
		,
		class_<param>("param")
			.def(constructor<const std::string&>())
			.def("ClassName", &param::ClassName)
			.def("GetName", LUA_CMEMFN(std::string, param, Name, void))
			.def("SetName", LUA_MEMFN(void, param, Name, const std::string&))
			.def("GetGrib2Number", LUA_CMEMFN(long, param, GribParameter, void))
			.def("SetGrib2Number", LUA_MEMFN(void, param, GribParameter, long))
			.def("GetGrib2Discipline", LUA_CMEMFN(long, param, GribDiscipline, void))
			.def("SetGrib2Discipline", LUA_MEMFN(void, param, GribDiscipline, long))
			.def("GetGrib2Category", LUA_CMEMFN(long, param, GribCategory, void))
			.def("SetGrib2Category", LUA_MEMFN(void, param, GribCategory, long))
			.def("GetGrib1Parameter", LUA_CMEMFN(long, param, GribIndicatorOfParameter, void))
			.def("SetGrib1Parameter", LUA_MEMFN(void, param, GribIndicatorOfParameter, long))
			.def("GetGrib1TableVersion", LUA_CMEMFN(long, param, GribTableVersion, void))
			.def("SetGrib1TableVersion", LUA_MEMFN(void, param, GribTableVersion, long))
			.def("GetUnivId", LUA_CMEMFN(unsigned long, param, UnivId, void))
			.def("SetUnivId", LUA_MEMFN(void, param, UnivId, unsigned long))
			.def("GetAggregation", LUA_CMEMFN(const aggregation&, param, Aggregation, void))
			.def("SetAggregation", LUA_MEMFN(void, param, Aggregation, const aggregation&))
		,
		class_<level>("level")
			.def(constructor<HPLevelType, double>())
			.def("ClassName", &level::ClassName)
			.def(tostring(self))
			.def("GetType", LUA_CMEMFN(HPLevelType, level, Type, void))
			.def("SetType", LUA_MEMFN(void, level, Type, HPLevelType))
			.def("GetValue", LUA_CMEMFN(double, level, Value, void))
			.def("SetValue", LUA_MEMFN(void, level, Value, double))
		,
		class_<raw_time>("raw_time")
			.def(constructor<const std::string&>())
			.def("ClassName", &raw_time::ClassName)
			.def("String", LUA_CMEMFN(std::string, raw_time, String, const std::string&))
			.def("Adjust", &raw_time::Adjust)
			.def("Empty", &raw_time::Empty)
		,
		class_<forecast_time>("forecast_time")
			.def(constructor<const raw_time&, const raw_time&>())
			.def("ClassName", &forecast_time::ClassName)
			.def("GetOriginDateTime", LUA_MEMFN(raw_time&, forecast_time, OriginDateTime, void), return_reference_to(_1))
			.def("GetValidDateTime", LUA_MEMFN(raw_time&, forecast_time, ValidDateTime, void))
			.def("SetOriginDateTime", LUA_MEMFN(void, forecast_time, OriginDateTime, const std::string&, const std::string&))
			.def("SetValidDateTime", LUA_MEMFN(void, forecast_time, ValidDateTime, const std::string&, const std::string&))
			.def("GetStep", LUA_CMEMFN(int, forecast_time, Step, void))
			.def("GetStepResolution", LUA_CMEMFN(HPTimeResolution, forecast_time, StepResolution, void))
			.def("SetStepResolution", LUA_MEMFN(void, forecast_time, StepResolution, HPTimeResolution))
		,
		class_<point>("point")
			.def(constructor<double, double>())
			.def("ClassName", &point::ClassName)
			.def("SetX", LUA_MEMFN(void, point, X, double))
			.def("SetY", LUA_MEMFN(void, point, Y, double))
			.def("GetX", LUA_CMEMFN(double, point, X, void))
			.def("GetY", LUA_CMEMFN(double, point, Y, void))
		,
		class_<producer>("producer")
			.def(constructor<>())
			.def("ClassName", &producer::ClassName)
			.def("SetName", LUA_MEMFN(void, producer, Name, const std::string&))
			.def("GetName", LUA_CMEMFN(std::string, producer, Name, void))
			.def("SetId", LUA_MEMFN(void, producer, Id, long))
			.def("GetId", LUA_CMEMFN(long, producer, Id, void))
			.def("SetProcess", LUA_MEMFN(void, producer, Process, long))
			.def("GetProcess", LUA_CMEMFN(long, producer, Process, void))
			.def("SetCentre", LUA_MEMFN(void, producer, Centre, long))
			.def("GetCentre", LUA_CMEMFN(long, producer, Centre, void))
			// TableVersion intentionally left out since in RADON it will be only
			// a parameter property
		,
		class_<logger>("logger")
			.def(constructor<>())
			.def("Trace", &logger::Trace)
			.def("Debug", &logger::Debug)
			.def("Info", &logger::Info)
			.def("Warning", &logger::Warning)
			.def("Error", &logger::Error)
			.def("Fatal", &logger::Fatal)
		,
		class_<aggregation>("aggregation")
			.def(constructor<HPAggregationType, HPTimeResolution, int>())
			.def("ClassName", &aggregation::ClassName)
			.def("GetType", LUA_CMEMFN(HPAggregationType, aggregation, Type, void))
			.def("SetType", LUA_MEMFN(void, aggregation, Type, HPAggregationType))
			.def("GetTimeResolution", LUA_CMEMFN(HPTimeResolution, aggregation, TimeResolution, void))
			.def("SetTimeResolution", LUA_MEMFN(void, aggregation, TimeResolution, HPTimeResolution))
			.def("GetTimeResolutionValue", LUA_CMEMFN(int, aggregation, TimeResolutionValue, void))
			.def("SetTimeResolutionValue", LUA_MEMFN(void, aggregation, TimeResolutionValue, int))
		,
		class_<configuration, std::shared_ptr<configuration>>("configuration")
			.def(constructor<>())
			.def("ClassName", &configuration::ClassName)
			.def("GetOutputFileType", LUA_CMEMFN(HPFileType, configuration, OutputFileType, void))
			.def("GetSourceProducer", LUA_CMEMFN(const producer&, configuration, SourceProducer, size_t))
			.def("GetTargetProducer", LUA_CMEMFN(const producer&, configuration, TargetProducer, void))
			.def("GetForecastStep", &configuration::ForecastStep)

		,
		class_<plugin_configuration, configuration, std::shared_ptr<plugin_configuration>>("plugin_configuration")
			.def(constructor<>())
			.def("ClassName", &plugin_configuration::ClassName)
			.def("GetValue", &plugin_configuration::GetValue)
			.def("GetValue", &plugin_configuration::GetValueList)
			.def("Exists", &plugin_configuration::Exists)
		,
		class_<lcl_t>("lcl_t")
			.def(constructor<>())
			.def_readwrite("T", &lcl_t::T)
			.def_readwrite("P", &lcl_t::P)
			.def_readwrite("Q", &lcl_t::Q)
		,
		// util namespace
		def("Filter2D", &util::Filter2D)
		,
		// metutil namespace
		def("LCL_", &metutil::LCL_)
		,
		def("Es_", &metutil::Es_)
		,
		def("Gammas_", &metutil::Gammas_)
		,
		def("Gammaw_", &metutil::Gammaw_)
		,
		def("MixingRatio_", &metutil::MixingRatio_)
		,
		def("MoistLift_", &metutil::MoistLift_)
		,
		def("DryLift_", &metutil::DryLift_)
	];

	// STL
	
	module(L)
	[
		class_<std::vector<double>>("vector")
			.def(constructor<>())
			.def("size", &std::vector<double>::size)
	];

}

void BindPlugins(lua_State* L)
{
	module(L) [
		class_<compiled_plugin_base>("compiled_plugin_base")
			.def(constructor<>())
			.def("WriteToFile", LUA_CMEMFN(void, compiled_plugin_base, WriteToFile, const info&))
		,
		class_<luatool, compiled_plugin_base>("luatool")
			.def(constructor<>())
			.def("ClassName", &luatool::ClassName)
			.def("Fetch", LUA_CMEMFN(std::shared_ptr<info>, luatool, Fetch, const forecast_time&, const level&, const param&))
		,
		class_<hitool, std::shared_ptr<hitool>>("hitool")
			.def(constructor<>())
			.def("ClassName", &hitool::ClassName)
			// Local functions to luatool
			.def("VerticalMaximumMultiParamGrid", &hitool_wrapper::VerticalMaximumMultiParamGrid)
			.def("VerticalMaximumMultiParam", &hitool_wrapper::VerticalMaximumMultiParam)
			.def("VerticalMaximumGrid", &hitool_wrapper::VerticalMaximumGrid)
			.def("VerticalMaximum", &hitool_wrapper::VerticalMaximum)
			.def("VerticalMinimumMultiParamGrid", &hitool_wrapper::VerticalMinimumMultiParamGrid)
			.def("VerticalMinimumMultiParam", &hitool_wrapper::VerticalMinimumMultiParam)
			.def("VerticalMinimumGrid", &hitool_wrapper::VerticalMinimumGrid)
			.def("VerticalMinimum", &hitool_wrapper::VerticalMinimum)
			.def("VerticalSumMultiParamGrid", &hitool_wrapper::VerticalSumMultiParamGrid)
			.def("VerticalSumMultiParam", &hitool_wrapper::VerticalSumMultiParam)
			.def("VerticalSumGrid", &hitool_wrapper::VerticalSumGrid)
			.def("VerticalSum", &hitool_wrapper::VerticalSum)
			.def("VerticalAverageMultiParamGrid", &hitool_wrapper::VerticalAverageMultiParamGrid)
			.def("VerticalAverageMultiParam", &hitool_wrapper::VerticalAverageMultiParam)
			.def("VerticalAverageGrid", &hitool_wrapper::VerticalAverageGrid)
			.def("VerticalAverage", &hitool_wrapper::VerticalAverage)
			.def("VerticalCountMultiParamGrid", &hitool_wrapper::VerticalCountMultiParamGrid)
			.def("VerticalCountMultiParam", &hitool_wrapper::VerticalCountMultiParam)
			.def("VerticalCountGrid", &hitool_wrapper::VerticalCountGrid)
			.def("VerticalCount", &hitool_wrapper::VerticalCount)
			.def("VerticalHeightMultiParamGrid", &hitool_wrapper::VerticalHeightMultiParamGrid)
			.def("VerticalHeightMultiParam", &hitool_wrapper::VerticalHeightMultiParam)
			.def("VerticalHeightGrid", &hitool_wrapper::VerticalHeightGrid)
			.def("VerticalHeight", &hitool_wrapper::VerticalHeight)
			.def("VerticalValueMultiParamGrid", &hitool_wrapper::VerticalValueMultiParamGrid)
			.def("VerticalValueMultiParam", &hitool_wrapper::VerticalValueMultiParam)
			.def("VerticalValueGrid", &hitool_wrapper::VerticalValueGrid)
			.def("VerticalValue", &hitool_wrapper::VerticalValue)
			.def("VerticalPlusMinusAreaMultiParamGrid", &hitool_wrapper::VerticalPlusMinusAreaMultiParamGrid)
			.def("VerticalPlusMinusAreaMultiParam", &hitool_wrapper::VerticalPlusMinusAreaMultiParam)
			.def("VerticalPlusMinusAreaGrid", &hitool_wrapper::VerticalPlusMinusAreaGrid)
			.def("VerticalPlusMinusArea", &hitool_wrapper::VerticalPlusMinusArea)
	];
}

void luatool::Run(info_t myTargetInfo, unsigned short threadIndex)
{
	while (AdjustLeadingDimension(myTargetInfo))
	{
		ResetNonLeadingDimension(myTargetInfo);

		while (AdjustNonLeadingDimension(myTargetInfo))
		{
			Calculate(myTargetInfo, threadIndex);

			if (itsConfiguration->StatisticsEnabled())
			{
				itsConfiguration->Statistics()->AddToMissingCount(myTargetInfo->Data().MissingCount());
				itsConfiguration->Statistics()->AddToValueCount(myTargetInfo->Data().Size());
			}
		}
	}
}

void luatool::Finish() const
{

	if (itsConfiguration->StatisticsEnabled())
	{
		itsTimer->Stop();
		itsConfiguration->Statistics()->AddToProcessingTime(itsTimer->GetTime());
	}
}

std::shared_ptr<info> luatool::Fetch(const forecast_time& theTime, const level& theLevel, const param& theParam) const
{
	return compiled_plugin_base::Fetch(theTime,theLevel,theParam,false);
}
