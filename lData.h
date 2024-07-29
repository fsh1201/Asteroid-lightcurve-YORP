#ifndef LDATAH
#define LDATAH

#include <iostream>
#include <string>
#include "lightCurve.h"

class lData
{
private:
	int n;	//length of the light curve
	std::vector<lightCurve> vlc;
public:
	int get_n() { return n; }
	std::vector<lightCurve> get_vlc() { return vlc; }
	lData(std::string lcFileName);
	lData(double*** lc, int nlc, int** nlp);
	lData();
	lData(const lData& ld);
	/*period [d]*/
	double period(int th = 8, double s = 0.02, double e = 1.6, double step = 0.0000001, int m = 8);
	void normalize(int method);
	lData& operator=(const lData& ld);
	std::vector<std::vector<double>> relative_residual(lData ld);
	double relative_rms(lData ld);
	friend std::ostream& operator<<(std::ostream& out, lData ld)
	{
		out << ld.n << std::endl;
		for (int i = 0; i < ld.n; i++)
		{
			out << ld.vlc[i];
		}
		return out;
	}
	void mergelc(std::string outFileName);
	lData separate_day(const lData& ld);
	lData all_data();
	lData bootstrap(int num = 0);
};

std::pair<double, double> period_t(std::vector<std::vector<double>> vvlc, double s = 0.02, double e = 1.6, double step = 1e-6, int m = 8);

class my_greater
{
public:
	bool operator()(const std::pair<double, double>& p1, const std::pair<double, double>& p2)
	{
		return p1.second > p2.second;
	}
};

class timecmp
{
public:
	bool operator()(const dataPoint& dp1, const dataPoint& dp2)
	{
		return dp1.time < dp2.time;
	}
};

#endif //LDATAH