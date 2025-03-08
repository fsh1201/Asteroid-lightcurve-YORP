/* 
	This program takes the light curve and the range of initial parameters as the input and 
	output the results of out areas, lightcurves, parameters and related chi2. 
	This program works with 40 processes by default, and calculate with totally 1000
	random initial parameters.
	
	The range of parameters can be pass through input_parameter file.
	If the input_parameter file does not exsit or cannot be open, then this program will adopt
	the default parameters already set in this yorp.cpp file.

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

/*default parameters*/
uniform_real_distribution<double> dislambda(0., 360.);
uniform_real_distribution<double> disbeta(-90., 90.);
uniform_real_distribution<double> disp(2., 10.);
uniform_real_distribution<double> disnu(-1e-7, 1e-7);
uniform_real_distribution<double> disfa(0.009, 0.048);
uniform_real_distribution<double> disfd(3.3 * PI / 180.0, 9.45 * PI / 180.0);
uniform_real_distribution<double> disfk(-7.7e-4 * 180.0 / PI, -2.2e-4 * 180.0 / PI);
uniform_real_distribution<double> disfc(0.05, 0.5);
uniform_real_distribution<double> disfb(0.0142, 0.087);
double JD0 = 0;
double phi0 = 0;
double con_reg = 0.1;
int lm_l = 6, lm_m = 6;
int rows = 8;
double stopcondition = 100;
int fflambda = 1, ffbeta = 1, ffp = 1, ffnu = 1, ffa = 0, ffd = 0, ffk = 0, ffb = 0, ffc = 0;	//fixed (0) or free (1)
int ntotal = 1000, nth = 40;
int nbootstrap = 10000, nbsth = 40;

vector<string> yorp(string lcFileName, string inputFileName, string areaFileName, string parFileName, string outlcFileName,
	vector<double> vlambda, vector<double> vbeta, vector<double> vP, vector<double> vnu, int nt)
{
	vector<string> vs;
	int n = vlambda.size();
	for (int i = 0; i < n; i++)
	{
		double lambda = vlambda[i];
		double beta = vbeta[i];
		double P = vP[i];
		double nu = vnu[i];
		char ch;
		string line;
		stringstream ssInputFileName("");
		ssInputFileName << setprecision(10) << inputFileName << "_" << nt << "_" << i << ".txt";
		fstream testfile(ssInputFileName.str(), ios::in);
		testfile >> ch;
		if ((!testfile.is_open()) || testfile.eof())
		{
			fstream outInFile(ssInputFileName.str(), ios::out);
			outInFile << setprecision(15);
			outInFile << lambda << "\t" << fflambda << "\tinital lambda[deg](0 / 1 - fixed / free)" << endl;
			outInFile << beta << "\t" << ffbeta << "\tinitial beta [deg] (0/1 - fixed/free)" << endl;
			outInFile << P << "\t" << ffp << "\tinital period [hours] (0/1 - fixed/free)" << endl;
			outInFile << nu << "\t" << ffnu << "\tinital YORP strength [rad/day^2] (0/1 - fixed/free)" << endl;
			outInFile << JD0 << "\tzero time [JD]" << endl;
			outInFile << phi0 << "\tinitial rotation angle [deg]" << endl;
			outInFile << con_reg << "\tconvexity regularization" << endl;
			outInFile << lm_l << " " << lm_m << "\tdegree and order of spherical harmonics expansion" << endl;
			outInFile << rows << "\tnumber of rows" << endl;
			outInFile << disfa(gen) << "\t" << ffa << "\tphase funct. param. 'a' (0/1 - fixed/free) " << endl;
			outInFile << disfd(gen) << "\t" << ffd << "\tphase funct.param. 'd' (0/1 - fixed/free)" << endl;
			outInFile << disfk(gen) << "\t" << ffk << "\tphase funct. param. 'k' (0/1 - fixed/free) " << endl;
			outInFile << disfb(gen) << "\t" << ffb << "\tphase funct. param. 'b' (0/1 - fixed/free) " << endl;
			outInFile << disfc(gen) << "\t" << ffc << "\tLambert coefficient 'c' (0/1 - fixed/free)" << endl;
			outInFile << "100\titeration stop condition" << endl;
			outInFile.close();
		}
		stringstream ssAreaFileName("");
		ssAreaFileName << setprecision(10) << areaFileName << "_" << nt << "_" << i << ".txt";
		stringstream ssParFileName("");
		ssParFileName << setprecision(10) << parFileName << "_" << nt << "_" << i << ".txt";
		stringstream ssOutLcFileName("");
		ssOutLcFileName << setprecision(10) << outlcFileName << "_" << nt << "_" << i << ".txt";
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
		line += " " + to_string(nt) + " " + to_string(i);
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

	ifstream infile(inputFileName, ios::in);
	if (infile)
	{
		double low, high;
		string line;
		infile >> low >> high >> fflambda;	getline(infile, line);	//lambda
		dislambda = uniform_real_distribution<double>(low, high);
		infile >> low >> high >> ffbeta;	getline(infile, line);	//beta
		disbeta = uniform_real_distribution<double>(low, high);
		infile >> low >> high >> ffp;		getline(infile, line);	//p
		disp = uniform_real_distribution<double>(low, high);
		infile >> low >> high >> ffnu;		getline(infile, line);	//yorp
		disnu = uniform_real_distribution<double>(low, high);
		infile >> JD0;						getline(infile, line);	//jd0
		infile >> phi0;						getline(infile, line);	//phi0
		infile >> con_reg;					getline(infile, line);	//convexity regularization
		infile >> lm_l >> lm_m;				getline(infile, line);	//degree and order of spherical harmonics expansion
		infile >> rows;						getline(infile, line);	//number of rows
		infile >> low >> high >> ffa;		getline(infile, line);	//phase funct. param. 'a'
		disfa = uniform_real_distribution<double>(low, high);
		infile >> low >> high >> ffd;		getline(infile, line);	//phase funct. param. 'd'
		disfd = uniform_real_distribution<double>(low, high);
		infile >> low >> high >> ffk;		getline(infile, line);	//phase funct. param. 'k'
		disfk = uniform_real_distribution<double>(low, high);
		infile >> low >> high >> ffb;		getline(infile, line);	//phase funct. param. 'b'
		disfb = uniform_real_distribution<double>(low, high);
		infile >> low >> high >> ffc;		getline(infile, line);	//phase funct. param. 'c'
		disfc = uniform_real_distribution<double>(low, high);
		infile >> stopcondition;			getline(infile, line);	//iteration stop condition
		infile >> ntotal;					getline(infile, line);	//number of total random parameter sets
		infile >> nth;						getline(infile, line);	//number of threads
	}
	else
	{
		cout << "failed to open the config file" << endl;
		cout << "using the default parameters" << endl;
	}

	fstream yorpFile(yorp_chi2FileName, ios::out);
	yorpFile << "# lambda(deg) beta(deg) p(hour) yorp(rad/d^2) rchisq jd0(JD) phi0(deg) a d k b c nt i" << endl;

	int n = ntotal / nth;

	vector<future<vector<string>>> vfvs;
	vector<double> vlambda, vbeta, vP, vnu;
	for (int i = 0; i < nth; i++)
	{
		for (int th = 0; th < n; th++)
		{
			vlambda.push_back(dislambda(gen));
			vbeta.push_back(disbeta(gen));
			vP.push_back(disp(gen));
			vnu.push_back(disnu(gen));
		}
		vfvs.push_back(async(yorp, lcFileName, inputFileName, areaFileName, parFileName, outlcFileName, vlambda, vbeta, vP, vnu, i));
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