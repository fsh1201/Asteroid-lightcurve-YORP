#ifndef POINTH
#define POINTF

#include "point.h"
#include <vector>

class lightCurve;

class dataPoint
{
private:
	double time;
	double brightness;
	point as;
	point ae;
	int n;
	std::vector<double> vtime;
	std::vector<double> vbrightness;
	std::vector<point> vas;
	std::vector<point> vae;
public:
	friend class lightCurve;
	friend class lData;
	friend class timecmp;
	double get_time() { return time; }
	double get_brightness() { return brightness; }
	point get_as() { return as; }
	point get_ae() { return ae; }
	int get_n() { return n; }
	std::vector<double> get_vtime() { return vtime; }
	std::vector<double> get_vbrightness() { return vbrightness; }
	std::vector<point> get_vas() { return vas; }
	std::vector<point> get_vae() { return vae; }
	dataPoint(double time = 0, double brightness = 0, point as = point(0, 0, 0), point ae = point(0, 0, 0));
	dataPoint(std::vector<double> vtime, std::vector<double> vbrightness, std::vector<point> vas, std::vector<point> vae);
	dataPoint(const dataPoint& dp);
	dataPoint& operator=(const dataPoint& dp);
	friend std::ostream& operator<<(std::ostream& out, const dataPoint& dp)
	{
		out << dp.time << " " << dp.brightness << " " << dp.as << " " << dp.ae;
		if (dp.n > 1)
		{
			out << dp.n << std::endl;
			for (int i = 0; i < dp.n; i++)
			{
				out << dp.vtime[i] << " " << dp.vbrightness[i] << " " << dp.vas[i] << " " << dp.vae[i] << std::endl;
			}
		}
		return out;
	}
	point PAB();
	double phaseAngle();
};
#endif //POINTH