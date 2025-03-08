/*
	This program takes the light curve as input and resample a new lightcurve to
	get inversion results. This program is used for uncertainty estimation.

	Syntax:
		yorp_error_t lc input_parameter out_area out_par out_lcs yorp_chi2 nbootstrap nth

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
#include "lData.h"

using namespace std;

vector<string> yorp(string lcFileName, string inputFileName, string areaFileName, string parFileName, string outlcFileName,
	int ns, int ne)
{
	vector<string> vs;
	for (int i = ns; i < ne; i++)
	{
		stringstream sscmd;
		stringstream bslcFileName, area, par, outlc;
		bslcFileName << lcFileName << "_fn_" << i << ".txt";
		area << areaFileName << "_fn_" << i << ".txt";
		par << parFileName << "_fn_" << i << ".txt";
		outlc << outlcFileName << "_fn_" << i << ".txt";
		sscmd << "./firstnlc " << lcFileName << " " << i << " > " << bslcFileName.str();
		cout << sscmd.str() << endl;
		system(sscmd.str().c_str());
		lData ld(bslcFileName.str());
		double jdx = ld.jdx();
		sscmd = stringstream("");
		sscmd << "cat " << bslcFileName.str() << " | ./convexinv -s -o " << area.str() << " -p " << par.str() << " "
			<< inputFileName << " " << outlc.str();
		cout << sscmd.str() << endl;
		system(sscmd.str().c_str());
		lData ldo(outlc.str());
		double rms = ld.relative_rms(ldo);
		fstream fp(par.str(), ios::in);
		string line;
		getline(fp, line);
		getline(fp, line);
		line += " " + to_string(jdx) + " " + to_string(rms);
		vs.push_back(line);
	}
	return vs;
}

int main(int argc, char** argv)
{
	if (argc < 8)
	{
		cout << "yorp_error_t lc input_parameter out_area out_par out_lcs yorp_chi2 nth" << endl;
		cout << "yorp_chi2: parameters and correlated chi2 file" << endl;

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
	
	int nth = stoi(argv[7]);
	lData ld(lcFileName);
	int min_nlc = ld.min_nlc();
	int bsn = ld.get_n() - min_nlc + 1;
	int n = bsn / nth;
	int nmod = bsn % nth;
	int m = 0;
	for (int i = 0; i < nth; i++)
	{
		if (i < nmod)
		{
			vfvs.push_back(async(launch::async, [=]() {
				return yorp(lcFileName, inputFileName, areaFileName, parFileName, outlcFileName, min_nlc + m, min_nlc + m + n + 1);
				}));
			m = m + n + 1;
		}
		else
		{
			vfvs.push_back(async(launch::async, [=]() {
				return yorp(lcFileName, inputFileName, areaFileName, parFileName, outlcFileName, min_nlc + m, min_nlc + m + n);
				}));
			m = m + n;
		}
	}
	for (int i = 0; i < nth; i++)
	{
		vvs.push_back(vfvs[i].get());
	}
	fstream yorpf(yorp_chi2FileName, ios::out);
	yorpf << "# lambda(deg) beta(deg) p(hour) yorp(rad/d^2) rchisq jd0(JD) phi0(deg) a d k b c jdx rms" << endl;
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
