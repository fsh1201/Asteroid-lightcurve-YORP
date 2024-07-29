#ifndef DATAPOINTH
#define DATAPOINTH

#include "dataPoint.h"
#include <cmath>
#include <vector>

namespace normal
{
	constexpr int max_min = 0;
	constexpr int mean = 1;
}

class lightCurve
{
private:
	int length;
	std::vector<dataPoint> vdp;
public:
	friend class lData;
	int get_length() { return length; }
	std::vector<dataPoint> get_vdp() { return vdp; }
	lightCurve(std::vector<dataPoint> lc);
	lightCurve(double** lc, int n);
	lightCurve(const lightCurve& lc);
	lightCurve();
	dataPoint& operator[](int k);
	const dataPoint& operator[](int k) const;
	void normalize(int method);
	lightCurve& operator=(const lightCurve& lc);
	friend std::ostream& operator<<(std::ostream& out, lightCurve lc)
	{
		out << lc.length << " 0" << std::endl;
		for (int i = 0; i < lc.length; i++)
		{
			out << lc.vdp[i] << std::endl;
		}
		return out;
	}
	double relativeChi2(lightCurve lc);
};

#endif //DATAPOINTH