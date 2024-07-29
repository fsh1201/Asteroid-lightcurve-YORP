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
	//��ԭ�����
	double r();
	//����ľ���
	double r(const point& p);
	//��λ����
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
	//�ڻ�
	double innerProduct(const point& p);
	//���
	point crossProduct(const point& p);
	//�Ƕ�����
	double cosAngleP(point p);
	//�Ƕ�
	double angleP(point p);
	//��x��ת
	point rotateX(double angle);
	//��y��ת
	point rotateY(double angle);
	//��z��ת
	point rotateZ(double angle);
	//����
	double longitude();
	//γ��
	double latitude();
};

#endif //POINTH