/*
 * playground.h
 *
 */

#ifndef TRANSFORMER_H
#define TRANSFORMER_H

#include "compiled_plugin.h"
#include "compiled_plugin_base.h"
#include "playground.cuh"
#include <boost/property_tree/ptree.hpp>
#include <vector>

namespace himan
{
namespace plugin
{

class playground : public compiled_plugin, private compiled_plugin_base
{
public:
	playground();

	inline virtual ~playground() {}
	playground(const playground& other) = delete;
	playground& operator=(const playground& other) = delete;

	virtual void Process(std::shared_ptr<const plugin_configuration> conf);

	virtual std::string ClassName() const { return "himan::plugin::playground"; }
	virtual HPPluginClass PluginClass() const { return kCompiled; }
	virtual HPVersionNumber Version() const { return HPVersionNumber(1, 2); }
private:
	virtual void Calculate(std::shared_ptr<info> theTargetInfo, unsigned short theThreadIndex);

	// Check and write json parameters needed for playground plug-in to local variables.
	void SetAdditionalParameters();
	std::vector<level> LevelsFromString(const std::string& levelType, const std::string& levelValues) const;

#ifdef HAVE_CUDA
	std::unique_ptr<playground_cuda::options> CudaPrepare(std::shared_ptr<info> myTargetInfo,
	                                                       std::shared_ptr<info> sourceInfo);
#endif

	double itsBase;
	double itsScale;
	std::string itsSourceParam;
	std::string itsTargetParam;
	int itsTargetUnivID;
	std::vector<level> itsSourceLevels;
	bool itsApplyLandSeaMask;
	double itsLandSeaMaskThreshold;
	HPInterpolationMethod itsInterpolationMethod;
};

// the class factory

extern "C" std::shared_ptr<himan_plugin> create() { return std::shared_ptr<playground>(new playground()); }
}  // namespace plugin
}  // namespace himan

#endif /* TRANSFORMER */
