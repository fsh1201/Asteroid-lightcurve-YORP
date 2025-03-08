/* lData library, a lData consists numbers of lightcurves */

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

#include "lData.h"
#include <fstream>
#include <algorithm>
#include <numeric>
#include <future>
#include <utility>
#include <random>
#include <chrono>

using namespace std;

lData::lData(std::string lcFileName)
{
	std::ifstream lcFile(lcFileName, std::ios::in);
	if (!lcFile.is_open())
	{
		std::cout << "lc File open failed" << std::endl;
		exit(-1);
	}
	n = 0;
	lcFile >> n;
	vlc.reserve(n);
	for (int i = 0; i < n; i++)
	{
		int m, mt;
		lcFile >> m >> mt;
		//lcFile >> m;
		//std::cout << m << " " << mt << std::endl;
		double temp[8];
		lightCurve lc;
		for (int j = 0; j < m; j++)
		{
			for (int k = 0; k < 8; k++)
			{
				lcFile >> temp[k];
			}
			string stemp;
			getline(lcFile, stemp);
			lc.vdp.push_back(dataPoint(temp[0], temp[1], point(temp[2], temp[3], temp[4]), point(temp[5], temp[6], temp[7])));
			lc.length = lc.vdp.size();
		}
		vlc.push_back(lc);
	}
}

lData::lData(double*** lc, int nlc, int** nlp)
{
	n = nlc;
	vlc.reserve(nlc);
	for (int i = 0; i < n; i++)
	{
		vlc.push_back(lightCurve(lc[i], nlp[i][0]));
	}
}

lData::lData()
{
	n = 0;
	vlc = std::vector<lightCurve>();
}

lData::lData(const lData& ld)
{
	n = ld.n;
	vlc = ld.vlc;
}

void lData::normalize(int method)
{
	if (method == normal::mean)
	{
		for (int i = 0; i < n; i++)
		{
			vlc[i].normalize(normal::mean);
		}
	}
	if (method == normal::max_min)
	{
		for (int i = 0; i < n; i++)
		{
			vlc[i].normalize(normal::max_min);
		}
	}
}

lData& lData::operator=(const lData& ld)
{
	this->n = ld.n;
	this->vlc = ld.vlc;
	return *this;
}

double lData::period(int th, double s, double e, double step, int m)
{
	lData ld(*this);
	ld.normalize(normal::mean);
	std::vector<std::vector<double>> lcn;
	for (int i = 0; i < ld.n; i++)
	{
		for (int j = 0; j < ld.vlc[i].length; j++)
		{
			std::vector<double> vtemp;
			vtemp.push_back(ld.vlc[i].vdp[j].time);
			vtemp.push_back(ld.vlc[i].vdp[j].brightness);
			lcn.push_back(vtemp);
		}
	}
	double s1 = s;
	double dse = (e - s) / th;
	std::vector<std::future<std::pair<double, double>>> vfdt;
	vfdt.resize(th);
	for (int i = 0; i < th; i++)
	{
		s1 = s + i * dse;
		vfdt[i] = std::async(std::launch::async, period_t, lcn, s1, s1 + dse, step, m);
	}
	std::vector<std::pair<double, double>> vdt;
	for (int i = 0; i < th; i++)
	{
		vdt.push_back(vfdt[i].get());
	}
	return (*std::max_element(vdt.begin(), vdt.end(), my_greater())).first;
}

std::pair<double, double> period_t(std::vector<std::vector<double>> vvlc, double s, double e, double step, int m)
{

	double temp = 1e26, T = 0;
	for (double t = s; t < e; t += step)
	{
		std::vector<std::vector<double>> lcn;
		lcn.resize(m);
		std::vector<double> mean;
		mean.resize(m);
		std::vector<double> v2;
		v2.resize(m);
		for (int i = 0; i < vvlc.size(); i++)
		{
			int k = (int)(fmod(vvlc[i][0], t) / (t / m));
			//std::cout << k << std::endl;
			if (k == m)
				k--;
			//std::cout << lcn[k].size() << std::endl;
			lcn[k].push_back(vvlc[i][1]);
		}
		for (int i = 0; i < m; i++)
		{
			mean[i] = std::accumulate(lcn[i].begin(), lcn[i].end(), 0.0) / lcn[i].size();
			v2[i] = 0;
			for (int j = 0; j < lcn[i].size(); j++)
			{
				v2[i] += pow(lcn[i][j], 2);
			}
			//std::cout << v2[i] << std::endl;
			v2[i] = v2[i] - lcn[i].size() * pow(mean[i], 2);

		}
		double sum_v2 = std::accumulate(v2.begin(), v2.end(), 0.0);
		//std::cout << sum_v2 << std::endl;
		if (temp > sum_v2)
		{
			temp = sum_v2;
			T = t;
		}
	}
	//std::cout << T << " " << temp << std::endl;
	return std::pair<double, double>(T, temp);
}

std::vector<std::vector<double>> lData::relative_residual(lData ld)
{
	lData ldn1(*this);
	lData ldn2(ld);
	ldn1.normalize(normal::mean);
	ldn2.normalize(normal::mean);
	std::vector<std::vector<double>> vvd;
	for (int i = 0; i < n; i++)
	{
		std::vector<double> vdre;
		for (int j = 0; j < vlc[i].length; j++)
		{
			vdre.push_back(ldn1.vlc[i].vdp[j].brightness - ldn2.vlc[i].vdp[j].brightness);
		}
		vvd.push_back(vdre);
	}
	return vvd;
}

double lData::relative_rms(lData ld)
{
	std::vector<std::vector<double>> vvd = this->relative_residual(ld);
	double rms = 0;
	int nnn = 0;
	for (int i = 0; i < vvd.size(); i++)
	{
		for (int j = 0; j < vvd[i].size(); j++)
		{
			rms += pow(vvd[i][j], 2);
		}
		nnn += vvd[i].size();
	}
	return sqrt(rms / nnn);
}

void lData::mergelc(string outFileName)
{
	ofstream outFile(outFileName, ios::out);
	outFile << n << endl;
	random_device seed;
	ranlux48 engine(seed());
	uniform_int_distribution<int> distribute(2, 5);
	for (int i = 0; i < n; i++)
	{
		int k = 0;

		for (int j = 0; j < vlc[i].length;)
		{
			int randint = distribute(engine);

		}
	}
}

lData lData::separate_day(const lData& ld)
{
	lData ld1;
	double t0 = ld.vlc[0].vdp[0].time;
	double t1;
	vector<dataPoint> vdp1;
	for (int i = 0; i < ld.n; i++)
	{
		for (int j = 0; j < ld.vlc[i].length; j++)
		{
			t1 = ld.vlc[i].vdp[j].time;
			if (t1 - t0 < 0.2)
				vdp1.push_back(ld.vlc[i].vdp[j]);
			else
			{
				ld1.vlc.push_back(lightCurve(vdp1));
				vdp1.clear();
				ld1.n++;
				vdp1.push_back(ld.vlc[i].vdp[j]);
			}
			t0 = t1;
		}
	}
	if (!vdp1.empty())
	{
		ld1.vlc.push_back(lightCurve(vdp1));
		ld1.n++;
	}
	return ld1;
}

lData lData::all_data()
{
	lData ld;
	ld.n = 1;
	ld.vlc.resize(1);
	ld.vlc[0].vdp.resize(0);
	for (auto i : vlc)
	{
		for (auto j : i.vdp)
		{
			ld.vlc[0].vdp.push_back(j);
		}
	}
	return ld;
}

lData lData::bootstrap(int num)
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine gen(seed);
	vector<int> vi;	//the position of each datapoint (in which light curve)
	for (int i = 0; i < vlc.size(); i++)
	{
		for (int j = 0; j < vlc[i].length; j++)
		{
			vi.push_back(i);
		}
	}
	lData ldq = all_data();
	uniform_int_distribution<int> disn(0, ldq.vlc[0].vdp.size() - 1);
	int nn = ldq.vlc[0].vdp.size();
	if (num != 0)
		nn = num;
	lData ld;
	ld.vlc.resize(n);
	for (int i = 0; i < nn; i++)
	{
		int x = disn(gen);
		ld.vlc[vi[x]].vdp.push_back(ldq.vlc[0].vdp[x]);
	}
	lData ldnew;
	for (auto i : ld.vlc)
	{
		if (i.vdp.empty())
			continue;
		sort(i.vdp.begin(), i.vdp.end(), timecmp());
		i.length = i.vdp.size();
		ldnew.vlc.push_back(i);
	}
	ldnew.n = ldnew.vlc.size();
	ldnew.normalize(normal::mean);
	return ldnew;
}

lData lData::bootstraplc(int num)
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine gen(seed);
	uniform_int_distribution<int> disn(0, n - 1);
	lData ldnew;
	int nn = n;
	if (num > 0)
		nn = num;
	ldnew.n = nn;
	for (int i = 0; i < nn; i++)
	{
		int x = disn(gen);
		ldnew.vlc.push_back(vlc[x]);
	}

	return ldnew;
}

lData lData::firstn(int n)
{
	lData nldata(*this);
	nldata.n = n;
	nldata.vlc = vector<lightCurve>(vlc.begin(), vlc.begin() + n);

	return nldata;
}

int lData::min_nlc(int n)
{
	int m = 0;
	for (int i = 0; i < n; i++)
	{
		m += vlc[i].length;
		if (m >= n)
			return i + 1;
	}
}

double lData::jd0()
{
	double jd = 1e26;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < vlc[i].length; j++)
		{
			jd = jd < vlc[i][j].time ? jd : vlc[i][j].time;
		}
	}

	return jd;
}

double lData::jdn()
{
	double jd = -1e26;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < vlc[i].length; j++)
		{
			jd = jd > vlc[i][j].time ? jd : vlc[i][j].time;
		}
	}

	return jd;
}

double lData::jdx()
{
	return jdn() - jd0();
}
