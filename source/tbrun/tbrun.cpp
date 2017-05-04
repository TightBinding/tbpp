/* ===========================================================================
 * Copyright (c) 2016-2017 Giacomo Resta
 *
 * This file is part of TightBinding++.
 *
 * TightBinding++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TightBinding++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * ===========================================================================
 */

#include <iostream>
#include <iomanip>
#include <tbpp/context.h>

#include <chrono>
#include <fstream>
#include <sstream>

using namespace std;
using namespace tbpp;
using time_point = std::chrono::time_point<std::chrono::steady_clock>;

void pse(time_point start, time_point end, const string& desc) {
    cout << desc << " dt="
         << chrono::duration_cast<chrono::milliseconds>(end-start).count()/1000.0
         << " sec" << endl;
}

void print_usage() {
    cout << "Usage: tbrun [FILE]...\n"
         << "Solves all TBPP nodes in FILEs which are currently unsolved\n\n";
}

bool file_exists(const string& filename) {
    ifstream file(filename);
    return file.good();
}

void hline() {
    for(size_t i=0; i<78; i++)
        cout << '=';
    cout << '\n';
}

int main(int argc, char *argv[])
{
    hline();
    cout << "TightBinding++ v" << TBPP_VERSION_MAJOR << "." << TBPP_VERSION_MINOR << "\n";
    hline();

    if(argc < 2) {
        print_usage();
        return 0;
    }

    for(size_t i=1; i<argc; i++) {
        string infile(argv[i]);

        // Determine whether extension is (.tbpp, .h5, .hdf5)
        // in which case use format: fileprefix_run_XXX.extension
        // otherwise format is: filename._run_XXX
        string out_pre = infile;
        string out_ext = "";
        size_t dot_index =  infile.rfind('.');
        if( dot_index != string::npos) {
            string ext = infile.substr(dot_index+1);
            if((ext == "tbpp") or (ext == "h5") or (ext == "hdf5")) {
                out_pre = infile.substr(0, dot_index);
                out_ext = "." + ext;
            }
        }

        // Determine output filename
        string outfile = out_pre + ".run_000" + out_ext;
        for(size_t j=0; file_exists(outfile); j++) {
            stringstream stream;
            stream << out_pre << ".run_" << setfill('0') << setw(3) << j << out_ext;
            outfile = stream.str();
        }

        cout << "Running " << infile << " -> " << outfile << '\n';

        Context ct;
        ct.load(infile);
        ct.make_ready();
        ct.save(outfile);
    }

    return 0;
}
