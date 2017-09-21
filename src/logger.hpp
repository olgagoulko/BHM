/*** LICENCE: ***
Bin histogram method for restoration of smooth functions from noisy integrals. Copyright (C) 2017 Olga Goulko

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or (at
your option) any later version.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
02110-1301 USA.

*** END OF LICENCE ***/
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
