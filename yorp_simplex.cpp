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
#include <algorithm>
#include <vector>
#include <iomanip>
#include <random>
#include <ctime>
#include <numeric>
#include <cmath>

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
int ntotal, nth;
int fflambda = 1, ffbeta = 1, ffp = 1, ffnu = 1, ffa = 0, ffd = 0, ffk = 0, ffb = 0, ffc = 0;	//fixed (0) or free (1)

struct par
{
	double lambda;
	double beta;
	double P;
	double nu;
	double fa;
	double fd;
	double fk;
	double fb;
	double fc;
	par(double lambda = 0, double beta = 0, double P = 0, double nu = 0, double fa = 0, double fd = 0, double fk = 0, double fb = 0, double fc = 0) :
		lambda(lambda), beta(beta), P(P), nu(nu), fa(fa), fd(fd), fk(fk), fb(fb), fc(fc) {};
	par(const par& p) :lambda(p.lambda), beta(p.beta), P(p.P), nu(p.nu), fa(p.fa), fd(p.fd), fk(p.fk), fb(p.fb), fc(p.fc) {};
	par centroid(vector<double> vlambda, vector<double> vbeta, vector<double> vP, vector<double> vnu, vector<double> vfa, vector<double> vfd,
		vector<double> vfk, vector<double> vfb, vector<double> vfc)
	{
		int n = vlambda.size();
		double lambda = accumulate(vlambda.begin(), vlambda.end(), 0.) / n;
		double beta = accumulate(vbeta.begin(), vbeta.end(), 0.) / n;
		double P = accumulate(vP.begin(), vP.end(), 0.) / n;
		double nu = accumulate(vnu.begin(), vnu.end(), 0.) / n;
		double fa = accumulate(vfa.begin(), vfa.end(), 0.) / n;
		double fd = accumulate(vfd.begin(), vfd.end(), 0.) / n;
		double fk = accumulate(vfk.begin(), vfk.end(), 0.) / n;
		double fb = accumulate(vfb.begin(), vfb.end(), 0.) / n;
		double fc = accumulate(vfc.begin(), vfc.end(), 0.) / n;
		return par(lambda, beta, P, nu, fa, fd, fk, fb, fc);
	}
	// 🔹 重载 + 运算符（par + par）
	par operator+(const par& other) const {
		return par(lambda + other.lambda, beta + other.beta, P + other.P, nu + other.nu,
			fa + other.fa, fd + other.fd, fk + other.fk, fb + other.fb, fc + other.fc);
	}

	// 🔹 重载 - 运算符（par - par）
	par operator-(const par& other) const {
		return par(lambda - other.lambda, beta - other.beta, P - other.P, nu - other.nu,
			fa - other.fa, fd - other.fd, fk - other.fk, fb - other.fb, fc - other.fc);
	}

	// 🔹 重载 * 运算符（par * par）
	par operator*(const par& other) const {
		return par(lambda * other.lambda, beta * other.beta, P * other.P, nu * other.nu,
			fa * other.fa, fd * other.fd, fk * other.fk, fb * other.fb, fc * other.fc);
	}

	// 🔹 重载 / 运算符（par / par）
	par operator/(const par& other) const {
		return par(lambda / other.lambda, beta / other.beta, P / other.P, nu / other.nu,
			fa / other.fa, fd / other.fd, fk / other.fk, fb / other.fb, fc / other.fc);
	}

	// 🔹 重载 + 运算符（par + double）
	par operator+(double value) const {
		return par(lambda + value, beta + value, P + value, nu + value,
			fa + value, fd + value, fk + value, fb + value, fc + value);
	}

	// 🔹 重载 - 运算符（par - double）
	par operator-(double value) const {
		return par(lambda - value, beta - value, P - value, nu - value,
			fa - value, fd - value, fk - value, fb - value, fc - value);
	}

	// 🔹 重载 * 运算符（par * double）
	par operator*(double value) const {
		return par(lambda * value, beta * value, P * value, nu * value,
			fa * value, fd * value, fk * value, fb * value, fc * value);
	}

	// 🔹 重载 / 运算符（par / double）
	par operator/(double value) const {
		return par(lambda / value, beta / value, P / value, nu / value,
			fa / value, fd / value, fk / value, fb / value, fc / value);
	}

	// 🔹 重载 << 运算符（输出流）
	friend ostream& operator<<(ostream& os, const par& p) {
		os << "par(lambda: " << p.lambda << ", beta: " << p.beta << ", P: " << p.P
			<< ", nu: " << p.nu << ", fa: " << p.fa << ", fd: " << p.fd
			<< ", fk: " << p.fk << ", fb: " << p.fb << ", fc: " << p.fc << ")";
		return os;
	}
};



string yorp(string lcFileName, string inputFileName, string areaFileName, string parFileName, string outlcFileName,
	par& parin, int nt)
{
	vector<string> vs;

	char ch;
	string line;
	stringstream ssInputFileName("");
	ssInputFileName << setprecision(10) << inputFileName << "_" << nt << ".txt";
	fstream testfile(ssInputFileName.str(), ios::in);
	testfile >> ch;
	if ((!testfile.is_open()) || testfile.eof())
	{
		fstream outInFile(ssInputFileName.str(), ios::out);
		outInFile << setprecision(15);
		outInFile << parin.lambda << "\t" << fflambda << "\tinital lambda[deg](0 / 1 - fixed / free)" << endl;
		outInFile << parin.beta << "\t" << ffbeta << "\tinitial beta [deg] (0/1 - fixed/free)" << endl;
		outInFile << parin.P << "\t" << ffp << "\tinital period [hours] (0/1 - fixed/free)" << endl;
		outInFile << parin.nu << "\t" << ffnu << "\tinital YORP strength [rad/day^2] (0/1 - fixed/free)" << endl;
		outInFile << JD0 << "\tzero time [JD]" << endl;
		outInFile << phi0 << "\tinitial rotation angle [deg]" << endl;
		outInFile << con_reg << "\tconvexity regularization" << endl;
		outInFile << lm_l << " " << lm_m << "\tdegree and order of spherical harmonics expansion" << endl;
		outInFile << rows << "\tnumber of rows" << endl;
		outInFile << parin.fa << "\t" << ffa << "\tphase funct. param. 'a' (0/1 - fixed/free) " << endl;
		outInFile << parin.fd << "\t" << ffd << "\tphase funct.param. 'd' (0/1 - fixed/free)" << endl;
		outInFile << parin.fk << "\t" << ffk << "\tphase funct. param. 'k' (0/1 - fixed/free) " << endl;
		outInFile << parin.fb << "\t" << ffb << "\tphase funct. param. 'b' (0/1 - fixed/free) " << endl;
		outInFile << parin.fc << "\t" << ffc << "\tLambert coefficient 'c' (0/1 - fixed/free)" << endl;
		outInFile << "100\titeration stop condition" << endl;
		outInFile.close();
	}
	stringstream ssAreaFileName("");
	ssAreaFileName << setprecision(10) << areaFileName << "_" << nt << ".txt";
	stringstream ssParFileName("");
	ssParFileName << setprecision(10) << parFileName << "_" << nt << ".txt";
	stringstream ssOutLcFileName("");
	ssOutLcFileName << setprecision(10) << outlcFileName << "_" << nt << ".txt";
	fstream areaFile(ssAreaFileName.str(), ios::in); areaFile >> ch;
	fstream parFile(ssParFileName.str(), ios::in); parFile >> ch;
	fstream outlcFile(ssOutLcFileName.str(), ios::in);	outlcFile >> ch;
	if ((!areaFile.is_open()) || (!parFile.is_open()) || (!outlcFile.is_open()) || areaFile.eof() || parFile.eof() || outlcFile.eof())
	{
		stringstream cmd("");
		cmd << "cat " << lcFileName << " | ./convexinv -s -o " << ssAreaFileName.str() << " -p " << ssParFileName.str() 
			<< " " << ssInputFileName.str() << " " << ssOutLcFileName.str();
		//cout << cmd.str() << endl;
		system(cmd.str().c_str());
	}
	parFile = fstream(ssParFileName.str(), ios::in);
	getline(parFile, line);
	getline(parFile, line);

	return line;
}

pair<par, double> nelder_mead(string lcFileName, string inputFileName, string areaFileName, string parFileName, string outlcFileName, fstream& yorpFile,
	vector<par>& simplex, vector<double>& chi2, double tol = 1e-6, int max_iter = 100)
{
	int n = simplex.size();
	double alpha = 1.0;   // 反射系数
	double gamma = 2.0;   // 扩展系数
	double rho = 0.5;     // 收缩系数
	double sigma = 0.5;   // 缩小系数

	int m = n;
	string line;
	stringstream ssline;

	// 1️⃣ **排序**：按 chi2 排序
	vector<int> indices(n);
	iota(indices.begin(), indices.end(), 0);
	sort(indices.begin(), indices.end(), [&](int i, int j) { return chi2[i] < chi2[j]; });

	for (int iter = 0; iter < max_iter; ++iter)
	{

		

		// 2️⃣ **计算质心**（去掉最差的点）
		par centroid;
		for (int i = 0; i < n - 1; ++i)
		{
			centroid = centroid + simplex[indices[i]];
		}
		centroid = centroid * (1.0 / (n - 1));

		// 3️⃣ **反射**
		par xr = centroid + (centroid - simplex[indices[n - 1]]) * alpha;
		xr.lambda = fmod(xr.lambda + 720, 360);
		xr.beta = asin(sin(xr.beta * PI / 180)) * 180 / PI;
		line = yorp(lcFileName, inputFileName, areaFileName, parFileName, outlcFileName, xr, m);
		m += 1;
		double jd0r, phi0r;
		double fr;
		ssline = stringstream(line);
		for (int i = 0; i < 5; i++)
		{
			ssline >> fr;
		}
		ssline >> jd0r >> phi0r;

		if (fr < chi2[indices[n - 2]])
		{
			// 替换最差点
			simplex[indices[n - 1]] = xr;
			chi2[indices[n - 1]] = fr;
		}
		else
		{
			// 4️⃣ **收缩**
			par xc = centroid + (simplex[indices[n - 1]] - centroid) * rho;
			xc.lambda = fmod(xc.lambda + 720, 360);
			xc.beta = asin(sin(xc.beta * PI / 180)) * 180 / PI;
			line = yorp(lcFileName, inputFileName, areaFileName, parFileName, outlcFileName, xc, m);
			m += 1;
			ssline = stringstream(line);
			double fc;
			for (int i = 0; i < 5; i++)
			{
				ssline >> fc;
			}
			ssline >> jd0r >> phi0r;

			if (fc < chi2[indices[n - 1]])
			{
				simplex[indices[n - 1]] = xc;
				chi2[indices[n - 1]] = fc;
			}
			else
			{
				// 5️⃣ **缩小**
				for (int i = 1; i < n; ++i)
				{
					simplex[indices[i]] = simplex[indices[0]] + (simplex[indices[i]] - simplex[indices[0]]) * sigma;
					line = yorp(lcFileName, inputFileName, areaFileName, parFileName, outlcFileName, simplex[indices[i]], m);
					m += 1;
					ssline = stringstream(line);
					for (int j = 0; j < 5; j++)
					{
						ssline >> chi2[indices[i]];
					}
					ssline >> jd0r >> phi0r;
				}
			}
		}

		iota(indices.begin(), indices.end(), 0);
		sort(indices.begin(), indices.end(), [&](int i, int j) { return chi2[i] < chi2[j]; });

		yorpFile << setprecision(15) << simplex[indices[0]].lambda << " " << simplex[indices[0]].beta << " " << simplex[indices[0]].P << " " <<
			simplex[indices[0]].nu << " " << chi2[indices[0]] << " " << jd0r << " " << phi0r << " " << simplex[indices[0]].fa << " " <<
			simplex[indices[0]].fd << " " << simplex[indices[0]].fk << " " << simplex[indices[0]].fb << " " << simplex[indices[0]].fc << endl;

		// 6️⃣ **检查收敛性**
		double max_diff = 0;
		for (int i = 1; i < n; ++i)
		{
			max_diff = max(max_diff, fabs(chi2[indices[i]] - chi2[indices[0]]));
		}
		if (max_diff < tol)
		{
			cout << "Converged after " << iter << " iterations." << endl;
			break;
		}
	}
	m += 1;

	iota(indices.begin(), indices.end(), 0);
	sort(indices.begin(), indices.end(), [&](int i, int j) { return chi2[i] < chi2[j]; });

	
	return make_pair(simplex[indices[0]], chi2[indices[0]]); // 返回最佳参数
}

pair<par, double> simplex(int argc, char** argv, int num)
{
	string lcFileName(argv[1]);
	string inputFileName(argv[2]);
	string areaFileName(argv[3]);
	string parFileName(argv[4]);
	string outlcFileName(argv[5]);
	string yorp_chi2FileName(argv[6]);
	string yorp_chi2_sorted;

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

	int n_par = fflambda + ffbeta + ffp + ffnu + ffa + ffd + ffk + ffb + ffc;	//number of parameters
	int n_init = n_par + 1;	//number of sets of initial parameters


	fstream yorpFile(yorp_chi2FileName + "_" + to_string(num) + ".txt", ios::out);
	yorpFile << "# lambda(deg) beta(deg) p(hour) yorp(rad/d^2) rchisq jd0(JD) phi0(deg) a d k b c" << endl;

	vector<future<string>> vfs;
	vector<par> vpar;
	par tpar;
	vector<double> vchi2;
	double chi2;

	for (int i = 0; i < n_init; i++)
	{
		vpar.push_back(par(dislambda(gen), disbeta(gen), disp(gen), disnu(gen), disfa(gen), disfd(gen), disfk(gen), disfb(gen), disfc(gen)));
	}
	for (int i = 0; i < n_init; i++)
	{
		vfs.push_back(async(yorp, lcFileName, inputFileName + "_" + to_string(num), areaFileName + "_" + to_string(num),
			parFileName + "_" + to_string(num), outlcFileName + "_" + to_string(num), ref(vpar[i]), i));
	}
	for (int i = 0; i < vfs.size(); i++)
	{
		string vs = vfs[i].get();
		yorpFile << vs << endl;
		stringstream tline(vs);
		for (int k = 0; k < 5; k++)
		{
			tline >> chi2;
		}
		vchi2.push_back(chi2);
	}


	pair<par, double> bestpar_chi2 = nelder_mead(lcFileName, inputFileName + "_" + to_string(num), areaFileName + "_" + to_string(num),
		parFileName + "_" + to_string(num), outlcFileName + "_" + to_string(num), yorpFile, vpar, vchi2, 1e-6, 10);
	yorpFile.close();
	stringstream cmd("");
	cmd << "sort -k 5n " << yorp_chi2FileName + "_" + to_string(num) + ".txt" << " > " 
		<< yorp_chi2FileName + "_" + to_string(num) + ".txt" << ".sort.txt";
	system(cmd.str().c_str());

	return bestpar_chi2;
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

	int n_par = fflambda + ffbeta + ffp + ffnu + ffa + ffd + ffk + ffb + ffc;	//number of parameters
	int n_init = n_par + 1;	//number of sets of initial parameters
	n_init = 10;
	
	vector<future<pair<par, double>>> vfp;
	pair<par, double> par_chi2;
	vector<par> vpar;
	vector<double> vchi2;
	for (int i = 0; i < n_init; i++)
	{
		vfp.push_back(async(simplex, argc, argv, i));
	}
	for (int i = 0; i < n_init; i++)
	{
		par_chi2 = vfp[i].get();
		vpar.push_back(par_chi2.first);
		vchi2.push_back(par_chi2.second);
	}
	fstream yorpFile(yorp_chi2FileName + "_" + to_string(n_init) + ".txt", ios::out);
	vector<int> indices(n_init);
	iota(indices.begin(), indices.end(), 0);
	sort(indices.begin(), indices.end(), [&](int i, int j) { return vchi2[i] < vchi2[j]; });
	vector<par> sortpar;
	vector<double> sortchi2;
	for (int i = 0; i < n_par + 1; i++)
	{
		sortpar.push_back(vpar[indices[i]]);
		sortchi2.push_back(vchi2[indices[i]]);
	}
	yorpFile << "# lambda(deg) beta(deg) p(hour) yorp(rad/d^2) rchisq jd0(JD) phi0(deg) a d k b c" << endl;
	par_chi2 = nelder_mead(lcFileName, inputFileName + "_" + to_string(n_init), areaFileName + "_" + to_string(n_init),
		parFileName + "_" + to_string(n_init), outlcFileName + "_" + to_string(n_init), yorpFile, sortpar, sortchi2, 1e-6, 10);
	
	//yorpFile << setprecision(15) << par_chi2.first.lambda << " " << par_chi2.first.beta << " " << par_chi2.first.P << " " << par_chi2.first.nu << " "
	//	<< par_chi2.second << " " << par_chi2.first.fa << " " << par_chi2.first.fd << " " << par_chi2.first.fk << " " << par_chi2.first.fb << " "
	//	<< par_chi2.first.fc << endl;
	yorpFile.close();
	stringstream cmd("");
	cmd << "sort -k 5n " << yorp_chi2FileName + "_" + to_string(n_init) + ".txt" << " > "
		<< yorp_chi2FileName + "_" + to_string(n_init) + ".txt" << ".sort.txt";
	system(cmd.str().c_str());

	return 0;
}