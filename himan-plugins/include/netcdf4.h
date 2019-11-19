/**
 * @file netcdf4.h
 *
 * @brief Class to implement netcdf writing and reading. Actual netcdf opening and reading is done by fminc4 library.
 */

#ifndef NETCDF4_H
#define NETCDF4_H

#include "auxiliary_plugin.h"
#include "info.h"

namespace himan
{
namespace plugin
{

class netcdf4 : public io_plugin
{
   public:
	netcdf4();

	//virtual ~netcdf4();

	netcdf4(const netcdf4& other) = delete;
	netcdf4& operator=(const netcdf4& other) = delete;

	virtual std::string ClassName() const
	{
		return "himan::plugin::netcdf4";
	};
	virtual HPPluginClass PluginClass() const
	{
		return kAuxiliary;
	};

	template <typename T>
	bool ToFile(info<T>& anInfo, std::string& outputFile, bool appendToFile = false);
	bool ToFile(info<double>& anInfo, std::string& outputFile, bool appendToFile = false);

	template <typename T>
	bool InitFile(const info<T>& baseInfo, const std::string& outputFile);

	template <typename T>
	std::vector<std::shared_ptr<himan::info<T>>> FromFile(const std::string& theInputFile, const search_options& options) const;
};

#ifndef HIMAN_AUXILIARY_INCLUDE

// the class factory

extern "C" std::shared_ptr<himan_plugin> create()
{
	return std::make_shared<netcdf4>();
}
#define HIMAN_AUXILIARY_INCLUDE
#endif /* HIMAN_AUXILIARY_INCLUDE */

}  // namespace plugin
}  // namespace himan

#endif /* NETCDF4_H */
