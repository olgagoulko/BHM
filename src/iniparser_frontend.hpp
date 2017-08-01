/** @file iniparser_frontend.hpp

    C++ front-end to iniparser library
*/

#ifndef INIPARSER_FRONTEND_HPP_34257acc45e34cfdb5efe2b72931f618
#define INIPARSER_FRONTEND_HPP_34257acc45e34cfdb5efe2b72931f618

#include <iniparser.h>
#include <string>
#include <stdexcept>
#include <algorithm> // for swap()

namespace iniparser {
    /// An (indefinite) error from the iniparser
    class Error : public std::runtime_error {
      public:
        Error(const std::string& msg) : std::runtime_error(msg) {}
    };

    class MissingFileError : public Error {
      public:
        std::string filename;
        MissingFileError(const std::string& msg, const std::string& fname)
            : Error(msg), filename(fname)
        {}
        ~MissingFileError() throw() {}
    };

    class param {
      private:
        std::string inifile_;
        dictionary* dict_;

        /// Assignment is not allowed
        param& operator=(const param&) { throw std::logic_error("Assignment of parameter objects is not allowed"); }
        
        /// Copy-constructing is not allowed
        param(const param&): dict_(0) { throw std::logic_error("Copy-constructing parameter objects is not allowed"); }
        
      public:
        /// Default constructor: no files loaded/parsed
        explicit param(): inifile_(), dict_(0)
        { }
        
        /// Load and parse a parameter file
        explicit param(const std::string& inifile): inifile_(inifile), dict_(0)
        {
            dict_=iniparser_load(inifile_.c_str());
            if (!dict_) throw MissingFileError("Error loading file: " + inifile_, inifile_);
            // iniparser_dump(dict_, stderr);
        }

        ~param()
        {
            if (dict_) iniparser_freedict(dict_);
        }

        void swap(param& other)
        {
            using std::swap;
            swap(inifile_, other.inifile_);
            swap(dict_, other.dict_);
        }
        
        void load(const std::string& inifile) {
            param other(inifile);
            swap(other);
        }
        
        int get(const std::string& key, int deflt)
        {
            if (!dict_) return deflt;
            return iniparser_getint(dict_, key.c_str(), deflt);
        }

        double get(const std::string& key, double deflt)
        {
            if (!dict_) return deflt;
            return iniparser_getdouble(dict_, key.c_str(), deflt);
        }

        bool get(const std::string& key, bool deflt)
        {
            if (!dict_) return deflt;
            return iniparser_getboolean(dict_, key.c_str(), deflt);
        }

        std::string get(const std::string& key, const std::string& deflt)
        {
            if (!dict_) return deflt;
            return iniparser_getstring(dict_, key.c_str(), deflt.c_str());
        }

        std::string get(const std::string& key, const char* deflt)
        {
            if (!dict_) return deflt;
            return iniparser_getstring(dict_, key.c_str(), deflt);
        }

    };
}



#endif /* INIPARSER_FRONTEND_HPP_34257acc45e34cfdb5efe2b72931f618 */
