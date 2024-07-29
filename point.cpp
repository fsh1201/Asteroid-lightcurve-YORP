/* point (vec3d) library */

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

#include "point.h"
#include <cmath>

point::point(double x, double y, double z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

point::point(const point& p)
{
	this->x = p.x;
	this->y = p.y;
	this->z = p.z;
}

double point::r()
{
	return sqrt(pow(this->x, 2) + pow(this->y, 2) + pow(this->z, 2));
}

double point::r(const point& p)
{
	return (*this - p).r();
}

point point::rn()
{
	double r = this->r();
	return point(x / r, y / r, z / r);
}

point point::operator+(const point& p)
{
	return point(x + p.x, y + p.y, z + p.z);
}

point point::operator-(const point& p)
{
	return point(x - p.x, y - p.y, z - p.z);
}

point point::operator*(const point& p)
{
	return point(x * p.x, y * p.y, z * p.z);
}

point point::operator/(const point& p)
{
	return point(x / p.x, y / p.y, z / p.z);
}

point& point::operator=(const point& p)
{
	this->x = p.x;
	this->y = p.y;
	this->z = p.z;
	return *this;
}

point point::operator+(const double& p)
{
	return point(x + p, y + p, z + p);
}
point point::operator-(const double& p)
{
	return point(x - p, y - p, z - p);
}
point point::operator*(const double& p)
{
	return point(x * p, y * p, z * p);
}
point point::operator/(const double& p)
{
	return point(x / p, y / p, z / p);
}

double point::innerProduct(const point& p)
{
	return this->x * p.x + this->y * p.y + this->z * p.z;
}

point point::crossProduct(const point& p)
{
	double nx = this->y * p.z - this->z * p.y;
	double ny = this->z * p.x - this->x * p.z;
	double nz = this->x * p.y - this->y * p.x;
	return point(nx, ny, nz);
}

double point::cosAngleP(point p)
{
	return this->rn().innerProduct(p.rn());
}

double point::angleP(point p)
{
	return acos(this->rn().innerProduct(p.rn()));
}

point point::rotateX(double angle)
{
	double nx = this->x;
	double ny = this->y * cos(angle) + this->z * sin(angle);
	double nz = -this->y * sin(angle) + this->z * cos(angle);
	return point(nx, ny, nz);
}

point point::rotateY(double angle)
{
	double nx = this->x * cos(angle) - this->z * sin(angle);
	double ny = this->y;
	double nz = this->x * sin(angle) + this->z * cos(angle);
	return point(nx, ny, nz);
}

point point::rotateZ(double angle)
{
	double nx = this->x * cos(angle) + this->y * sin(angle);
	double ny = -this->x * sin(angle) + this->y * cos(angle);
	double nz = this->z;
	return point(nx, ny, nz);
}

double point::longitude()
{
	return atan2(this->y, this->x);
}

double point::latitude()
{
	return atan2(this->z, sqrt(x * x + y * y));
}