/*
 * writer.cpp
 *
 *  Created on: Nov 26, 2012
 *	  Author: partio
 */

#include "writer.h"
#include "plugin_factory.h"
#include "logger_factory.h"
#include <fstream>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include "util.h"
#include "timer_factory.h"

#define HIMAN_AUXILIARY_INCLUDE

#include "grib.h"
#include "querydata.h"
#include "neons.h"
#include "cache.h"

#undef HIMAN_AUXILIARY_INCLUDE

using namespace himan::plugin;

writer::writer()
{
	itsLogger = std::unique_ptr<logger> (logger_factory::Instance()->GetLog("writer"));
}

bool writer::ToFile(std::shared_ptr<info> theInfo,
					std::shared_ptr<const plugin_configuration> conf,
					const std::string& theOutputFile)
{

	std::unique_ptr<himan::timer> t = std::unique_ptr<himan::timer> (timer_factory::Instance()->GetTimer());

	if (conf->StatisticsEnabled())
	{
		t->Start();
	}

	bool activeOnly = (conf->FileWriteOption() == kSingleFile) ? false : true;

	if (conf->UseCache())
	{
		std::shared_ptr<cache> c = std::dynamic_pointer_cast<plugin::cache> (plugin_factory::Instance()->Plugin("cache"));

		c->Insert(theInfo, activeOnly);
	}

	namespace fs = boost::filesystem;

	bool ret = false;

	std::string correctFileName = theOutputFile;
	
	HPFileWriteOption fileWriteOption = conf->FileWriteOption();
	HPFileType fileType = conf->OutputFileType();

	if ((fileWriteOption == kNeons || fileWriteOption == kMultipleFiles) || correctFileName.empty())
	{
		correctFileName = util::MakeFileName(fileWriteOption, theInfo);
	}

	fs::path pathname(correctFileName);

	if (!pathname.parent_path().empty() && !fs::is_directory(pathname.parent_path()))
	{
		fs::create_directories(pathname.parent_path());
	}

	switch (fileType)
	{

		case kGRIB:
		case kGRIB1:
		case kGRIB2:
		{

			std::shared_ptr<grib> theGribWriter = std::dynamic_pointer_cast<grib> (plugin_factory::Instance()->Plugin("grib"));

			correctFileName += ".grib";

			if (fileType == kGRIB2)
			{
				correctFileName += "2";
			}

			ret = theGribWriter->ToFile(theInfo, correctFileName, fileType, fileWriteOption);

			break;
		}
		case kQueryData:
		{
			std::shared_ptr<querydata> theWriter = std::dynamic_pointer_cast<querydata> (plugin_factory::Instance()->Plugin("querydata"));

			correctFileName += ".fqd";

			ret = theWriter->ToFile(theInfo, correctFileName, fileWriteOption);

			break;
		}
		case kNetCDF:
			break;

			// Must have this or compiler complains
		default:
			throw std::runtime_error(ClassName() + ": Invalid file type: " + HPFileTypeToString.at(fileType));
			break;

	}

	if (ret && fileWriteOption == kNeons)
	{
		std::shared_ptr<neons> n = std::dynamic_pointer_cast<neons> (plugin_factory::Instance()->Plugin("neons"));

		// Save file information to neons

		ret = n->Save(theInfo, correctFileName);

		if (!ret)
		{
			itsLogger->Warning("Saving file information to neons failed");
		   // unlink(correctFileName.c_str());
		}

	}

	if (!conf->UseCache())
	{
		theInfo->Grid()->Data()->Clear();
	}

	if (conf->StatisticsEnabled())
	{
		t->Stop();

		conf->Statistics()->AddToWritingTime(t->GetTime());
	}

	return ret;
}
