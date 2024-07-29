/*
	This program takes the light curve as input and resample a new lightcurve to
	get inversion results. This program is used for uncertainty estimation.

	Syntax:
		yorp lc input_parameter out_area out_par out_lcs yorp_chi2 nbootstrap nth

*/

/*
Copyright (C) 2024  Shuai Feng

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/

#include <iostream>
#include <iomanip>
#include <string>
#include <future>
#include <thread>
#include <sstream>
#include <vector>
#include <fstream>

using namespace std;

vector<string> yorp(string lcFileName, string inputFileName, string areaFileName, string parFileName, string outlcFileName, int nth, int n)
{
	vector<string> vs;
	for (int i = 0; i < n; i++)
	{
		stringstream sscmd;
		stringstream bslcFileName, area, par, outlc;
		bslcFileName << lcFileName << "_bs_" << nth << "_" << i << ".txt";
		area << areaFileName << "_bs_" << nth << "_" << i << ".txt";
		par << parFileName << "_bs_" << nth << "_" << i << ".txt";
		outlc << outlcFileName << "_bs_" << nth << "_" << i << ".txt";
		sscmd << "./bootstrap " << lcFileName << " > " << bslcFileName.str();
		cout << sscmd.str() << endl;
		system(sscmd.str().c_str());
		sscmd = stringstream("");
		sscmd << "cat " << bslcFileName.str() << " | ./convexinv -s -o " << area.str() << " -p " << par.str() << " " << inputFileName << " " << outlc.str();
		cout << sscmd.str() << endl;
		system(sscmd.str().c_str());
		fstream fp(par.str(), ios::in);
		string line;
		getline(fp, line);
		getline(fp, line);
		vs.push_back(line);
	}
	return vs;
}

int main(int argc, char** argv)
{
	if (argc < 6)
	{
		cout << "yorp lc input_parameter out_area out_par out_lcs yorp_chi2 nbootstrap nth" << endl;
		cout << "yorp_chi2: parameters and correlated chi2 file" << endl;
		cout << "nbootstrap: bootstrap time, 8000+ would be great" << endl;
		cout << "nth: count of thread to do this calculation" << endl;

		exit(-1);
	}
	string lcFileName(argv[1]);
	string inputFileName(argv[2]);
	string areaFileName(argv[3]);
	string parFileName(argv[4]);
	string outlcFileName(argv[5]);
	string yorp_chi2FileName(argv[6]);
	vector<future<vector<string>>> vfvs;
	vector<vector<string>> vvs;
	int bsn = stoi(argv[7]);
	int nth = stoi(argv[8]);
	int n = bsn / nth;
	for (int i = 0; i < nth; i++)
	{
		vfvs.push_back(async(launch::async, [=]() {
			return yorp(lcFileName, inputFileName, areaFileName, parFileName, outlcFileName, i, n);
			}));
	}
	for (int i = 0; i < nth; i++)
	{
		vvs.push_back(vfvs[i].get());
	}
	fstream yorpf(yorp_chi2FileName, ios::out);
	yorpf << "# lambda(deg) beta(deg) p(hour) yorp(rad/d^2) rchisq jd0(JD) phi0(deg) a d k b c" << endl;
	for (auto i : vvs)
	{
		for (auto j : i)
		{
			yorpf << j << endl;
		}
	}
	yorpf.close();

	return 0;
}