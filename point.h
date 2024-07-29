#ifndef POINTH
#define POINTH

#include <iostream>
class point
{
private:
	double x;
	double y;
	double z;
public:
	point(double x = 0, double y = 0, double z = 0);
	point(const point& p);
	double get_x() { return x; }
	double get_y() { return y; }
	double get_z() { return z; }
	//到原点距离
	double r();
	//到点的距离
	double r(const point& p);
	//单位向量
	point rn();
	point operator+(const point& p);
	point operator-(const point& p);
	point operator*(const point& p);
	point operator/(const point& p);
	point& operator=(const point& p);
	point operator+(const double& p);
	point operator-(const double& p);
	point operator*(const double& p);
	point operator/(const double& p);
	friend std::ostream& operator<<(std::ostream& out, const point& p)
	{
		out << p.x << " " << p.y << " " << p.z;
		return out;
	}
	//内积
	double innerProduct(const point& p);
	//叉积
	point crossProduct(const point& p);
	//角度余弦
	double cosAngleP(point p);
	//角度
	double angleP(point p);
	//绕x旋转
	point rotateX(double angle);
	//绕y旋转
	point rotateY(double angle);
	//绕z旋转
	point rotateZ(double angle);
	//经度
	double longitude();
	//纬度
	double latitude();
};

#endif //POINTH