/* 
	This program takes the light curve path and initial parameters as the input and 
	output the results of out areas, lightcurves, parameters and related chi2. 
	This program works with 40 processes by default, and each process calculate 250
	different initial random parameters in the range of the intervals set in lines 42-50 
	int this yorp.cpp file, while the numbers 40 and 250 can be changed in lines 131 and 133.

	Syntax:
		yorp lc input_parameter out_area out_par out_lcs yorp_chi2

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
#include <thread>
#include <fstream>
#include <future>
#include <sstream>
#include <vector>
#include <iomanip>
#include <random>
#include <ctime>

using namespace std;

const double PI = 3.14159265358979323846;

unsigned seed = chrono::system_clock::now().time_since_epoch().count();
default_random_engine gen(seed);

uniform_real_distribution<double> dislambda(17.23122838 - 0.001, 17.23122838 + 0.001);
uniform_real_distribution<double> disbeta(11.35821331 - 0.001, 11.35821331 + 0.001);
uniform_real_distribution<double> disp(360.0 / (1639.38928 + 0.00001) * 24.0, 360.0 / (1639.38928 - 0.00001) * 24.0);
uniform_real_distribution<double> disnu(-1.25e-9, 0.25e-9);
uniform_real_distribution<double> disfa(0.039, 0.041);
uniform_real_distribution<double> disfd(8.4 * PI / 180.0, 8.6 * PI / 180.0);
uniform_real_distribution<double> disfk(-4.8e-4 * 180.0 / PI, -4.6e-4 * 180.0 / PI);
uniform_real_distribution<double> disfc(0.11, 0.15);
uniform_real_distribution<double> disfb(0.049, 0.051);

vector<string> yorp(string lcFileName, string inputFileName, string areaFileName, string parFileName, string outlcFileName,
	vector<double> vlambda, vector<double> vbeta, vector<double> vP, vector<double> vnu)
{
	vector<string> vs;
	for (int i = 0; i < vlambda.size(); i++)
	{
		double lambda = vlambda[i];
		double beta = vbeta[i];
		double P = vP[i];
		double nu = vnu[i];
		char ch;
		string line;
		stringstream ssInputFileName("");
		ssInputFileName << setprecision(10) << inputFileName << "_" << lambda << "_" << beta << "_" << P << "_" << nu << ".txt";
		fstream testfile(ssInputFileName.str(), ios::in);
		testfile >> ch;
		if ((!testfile.is_open()) || testfile.eof())
		{
			fstream outInFile(ssInputFileName.str(), ios::out);
			outInFile << setprecision(15);
			outInFile << lambda << "\t0\tinital lambda [deg] (0/1 - fixed/free)" << endl;
			outInFile << beta << "\t0\tinitial beta [deg] (0/1 - fixed/free)" << endl;
			outInFile << P << "\t0\tinital period [hours] (0/1 - fixed/free)" << endl;
			outInFile << nu << "\t1\tinital YORP strength [rad/day^2] (0/1 - fixed/free)" << endl;
			outInFile << "2451545.0\tzero time [JD]" << endl;
			outInFile << "32.64\tinitial rotation angle [deg]" << endl;
			outInFile << "0.1\tconvexity regularization" << endl;
			outInFile << "6 6\tdegree and order of spherical harmonics expansion" << endl;
			outInFile << "8\tnumber of rows" << endl;
			outInFile << disfa(gen) << "\t0\tphase funct. param. 'a' (0/1 - fixed/free) " << endl;
			outInFile << disfd(gen) << "\t0\tphase funct.param. 'd' (0/1 - fixed/free)" << endl;
			outInFile << disfk(gen) << "\t0\tphase funct. param. 'k' (0/1 - fixed/free) " << endl;
			outInFile << disfb(gen) << "\t0\tphase funct. param. 'b' (0/1 - fixed/free) " << endl;
			outInFile << disfc(gen) << "\t0\tLambert coefficient 'c' (0/1 - fixed/free)" << endl;
			outInFile << "100\titeration stop condition" << endl;
			outInFile.close();
		}
		stringstream ssAreaFileName("");
		ssAreaFileName << setprecision(10) << areaFileName << "_" << lambda << "_" << beta << "_" << P << "_" << nu << ".txt";
		stringstream ssParFileName("");
		ssParFileName << setprecision(10) << parFileName << "_" << lambda << "_" << beta << "_" << P << "_" << nu << ".txt";
		stringstream ssOutLcFileName("");
		ssOutLcFileName << setprecision(10) << outlcFileName << "_" << lambda << "_" << beta << "_" << P << "_" << nu << ".txt";
		fstream areaFile(ssAreaFileName.str(), ios::in); areaFile >> ch;
		fstream parFile(ssParFileName.str(), ios::in); parFile >> ch;
		fstream outlcFile(ssOutLcFileName.str(), ios::in);	outlcFile >> ch;
		if ((!areaFile.is_open()) || (!parFile.is_open()) || (!outlcFile.is_open()) || areaFile.eof() || parFile.eof() || outlcFile.eof())
		{
			stringstream cmd("");
			cmd << "cat " << lcFileName << " | ./convexinv -s -o " << ssAreaFileName.str() << " -p " << ssParFileName.str() << " " << ssInputFileName.str() << " " << ssOutLcFileName.str();
			//cout << cmd.str() << endl;
			system(cmd.str().c_str());
		}
		parFile = fstream(ssParFileName.str(), ios::in);
		getline(parFile, line);
		getline(parFile, line);
		vs.push_back(line);
	}
	return vs;
	//return cmd.str();
}

int main(int argc, char** argv)
{
	if (argc < 2)
	{
		cout << "yorp lc input_parameter out_area out_par out_lcs yorp_chi2" << endl;
		exit(-1);
	}
	string lcFileName(argv[1]);
	string inputFileName(argv[2]);
	string areaFileName(argv[3]);
	string parFileName(argv[4]);
	string outlcFileName(argv[5]);
	string yorp_chi2FileName(argv[6]);
	fstream yorpFile(yorp_chi2FileName, ios::out);
	yorpFile << "# lambda(deg) beta(deg) p(hour) yorp(rad/d^2) rchisq jd0(JD) phi0(deg) a d k b c" << endl;
	vector<future<vector<string>>> vfvs;
	vector<double> vlambda, vbeta, vP, vnu;
	for (int i = 0; i < 40; i++)
	{
		for (int th = 0; th < 250; th++)
		{
			vlambda.push_back(dislambda(gen));
			vbeta.push_back(disbeta(gen));
			vP.push_back(disp(gen));
			vnu.push_back(disnu(gen));
		}
		vfvs.push_back(async(yorp, lcFileName, inputFileName, areaFileName, parFileName, outlcFileName, vlambda, vbeta, vP, vnu));
		vlambda.clear();
		vbeta.clear();
		vP.clear();
		vnu.clear();
	}
	for (int i = 0; i < vfvs.size(); i++)
	{
		vector<string> vs = vfvs[i].get();
		for (int j = 0; j < vs.size(); j++)
			yorpFile << vs[j] << endl;
	}
	yorpFile.close();
	stringstream cmd("");
	cmd << "sort -k 5n " << yorp_chi2FileName << " > " << yorp_chi2FileName << ".sort.txt";
	system(cmd.str().c_str());

	return 0;
}