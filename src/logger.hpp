/**
   @file logger.hpp defines a simplistic logger,
   as described at https://stackoverflow.com/a/8337882
*/
#ifndef BHM_SRC_LOGGER_HPP_6379cafa8eb74c1aaaca0a8fe9caad93
#define BHM_SRC_LOGGER_HPP_6379cafa8eb74c1aaaca0a8fe9caad93

#include <iostream>
#include <sstream>

namespace logger {
    class LogLine {
      private:
        std::stringstream stream_;
        std::ostream& out_;

        // this trick ensures that static `verbosity` variable is always initialized
        static int* verbosity_ptr_() {
            static int verbosity=0;
            return &verbosity;
        }

      public:
        LogLine(std::ostream& out = std::cout) : out_(out) {}

        ~LogLine() { // destructor does printing
            if (*verbosity_ptr_()==0) return;
            stream_ << "\n";
            out_ << stream_.rdbuf(); // copy the string stream to the output stream
            out_.flush();
        }

        static void set_verbosity(int verb) { *verbosity_ptr_()=verb; }
        
        template <class T>
        LogLine& operator<<(const T& val) { stream_ << val; return *this; }
    };
}
#define LOGGER logger::LogLine()
#define LOGGER_VERBOSITY(v) logger::LogLine::set_verbosity(v)

#endif /* BHM_SRC_LOGGER_HPP_6379cafa8eb74c1aaaca0a8fe9caad93 */
