/**
 * @file plugin_configuration.h
 *
 * @date Feb 11, 2013
 * @author partio
 */

#ifndef PLUGIN_CONFIGURATION_H
#define PLUGIN_CONFIGURATION_H

#include "himan_common.h"
#include <map>
#include "statistics.h"
#include "configuration.h"

namespace himan
{

class plugin_configuration : public configuration
{

public:

	plugin_configuration();
	plugin_configuration(const plugin_configuration& other);
	plugin_configuration& operator=(const plugin_configuration& other) = delete;

	plugin_configuration(const configuration& theConfiguration);
	plugin_configuration(const std::string& theName, const std::map<std::string,std::string>& theOptions);

	~plugin_configuration() = default;

    /**
     * @return Class name
     */

    std::string ClassName() const
    {
        return "himan::plugin_configuration";
    }

    std::ostream& Write(std::ostream& file) const;

	/**
	 * @brief Set plugin name
     * @param theName New name
     */

    void Name(const std::string& theName);

	/**
	 *
     * @return Plugin name
     */
	std::string Name() const;

	/**
	 * @brief Return a map of options defined for a plugin.
	 *
	 * Options are extra specifiers, for example with windvector one can specify
	 * 'for_ice' = true, these options are stored here.
	 *
	 * @return const reference to options map
	 */

    const std::map<std::string,std::string>& Options() const;
    void Options(const std::map<std::string,std::string>& theOptions);

	/**
	 * @brief Add new element to options map
     * @param key Key name
     * @param value Value
     */
	
    void AddOption(const std::string& key, const std::string& value);

	/**
	 * @brief Check if a key exists in options map
	 */

    bool Exists(const std::string & key) const;

	/**
	 * @brief Get value for a given key
     * @param key Key name
     * @return Value as string
     */

	std::string GetValue(const std::string & key) const;
    
	void Info(std::shared_ptr<info> theInfo);
	std::shared_ptr<info> Info() const;

	std::shared_ptr<statistics> Statistics() const;

	bool StatisticsEnabled() const;
	void StartStatistics();
	void WriteStatistics();
	
private:

	std::string itsName;
	std::map<std::string,std::string> itsOptions;
	std::shared_ptr<info> itsInfo;
	std::shared_ptr<statistics> itsStatistics;
};

inline
std::ostream& operator<<(std::ostream& file, const plugin_configuration& ob)
{
    return ob.Write(file);
}

} // namespace himan

#endif /* PLUGIN_CONFIGURATION_H */
