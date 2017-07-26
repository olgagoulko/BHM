/** @file iniparser_frontend.hpp

    C++ front-end to iniparser library
*/

#ifndef INIPARSER_FRONTEND_HPP_34257acc45e34cfdb5efe2b72931f618
#define INIPARSER_FRONTEND_HPP_34257acc45e34cfdb5efe2b72931f618

#include <iniparser.h>
#include <string>
#include <stdexcept>

namespace iniparser {
    class param {
      private:
        std::string inifile_;
        dictionary* dict_;
      public:
        /// Load and parse a parameter file
        explicit param(const std::string& inifile): inifile_(inifile), dict_(0)
        {
            dict_=iniparser_load(inifile_.c_str());
            if (!dict_) throw std::runtime_error("Error loading file: " + inifile_);
            // iniparser_dump(dict_, stderr);
        }

        ~param()
        {
            iniparser_freedict(dict_);
        }
        
        int get(const std::string& key, int deflt)
        {
            return iniparser_getint(dict_, key.c_str(), deflt);
        }

        double get(const std::string& key, double deflt)
        {
            return iniparser_getdouble(dict_, key.c_str(), deflt);
        }

        bool get(const std::string& key, bool deflt)
        {
            return iniparser_getboolean(dict_, key.c_str(), deflt);
        }

        std::string get(const std::string& key, const std::string& deflt)
        {
            return iniparser_getstring(dict_, key.c_str(), deflt.c_str());
        }

        std::string get(const std::string& key, const char* deflt)
        {
            return iniparser_getstring(dict_, key.c_str(), deflt);
        }

    };
}



#endif /* INIPARSER_FRONTEND_HPP_34257acc45e34cfdb5efe2b72931f618 */
