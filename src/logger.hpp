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
        std::ostream* out_ptr_;

        // this trick ensures that static `verbosity` variable is always initialized
        static int* verbosity_ptr_() {
            static int verbosity=0;
            return &verbosity;
        }

        // this trick ensures that static `stream` variable is always initialized
        static std::ostream*& default_ostream_ptr_() {
            static std::ostream* stream_ptr=0;
            return stream_ptr;
        }


      public:
        LogLine() : out_ptr_(default_ostream_ptr_()) {}
        LogLine(std::ostream& strm) : out_ptr_(&strm) {}

        ~LogLine() { // destructor does printing
            if (*verbosity_ptr_()==0) return;
            if (!out_ptr_) return;

            stream_ << "\n";
            *out_ptr_ << stream_.rdbuf(); // copy the string stream to the output stream
            out_ptr_->flush();
        }

        static void set_verbosity(int verb) { *verbosity_ptr_()=verb; }
        static void set_output_stream(std::ostream& ostrm) { default_ostream_ptr_()=&ostrm; }

        template <class T>
        LogLine& operator<<(const T& val) { stream_ << val; return *this; }
    };
}
#define LOGGER logger::LogLine()
#define LOGGER_VERBOSITY(v) logger::LogLine::set_verbosity(v)
#define LOGGER_OUTPUT(s) logger::LogLine::set_output_stream(s)

#endif /* BHM_SRC_LOGGER_HPP_6379cafa8eb74c1aaaca0a8fe9caad93 */
