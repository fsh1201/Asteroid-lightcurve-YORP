/* 
	datapoint library, a datapoint consists a time, a brightness, and the asteroidcentric 
	ecliptic coordinates of the sun and the earth.
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

#include "dataPoint.h"

dataPoint::dataPoint(double time, double brightness, point as, point ae)
{
	this->time = time;
	this->brightness = brightness;
	this->as = as;
	this->ae = ae;

	this->vtime.push_back(this->time);
	this->vbrightness.push_back(this->brightness);
	this->vas.push_back(this->as);
	this->vae.push_back(this->ae);

	this->n = 1;
}

dataPoint::dataPoint(std::vector<double> vtime, std::vector<double> vbrightness, std::vector<point> vas, std::vector<point> vae)
{
	this->n = vtime.size();

	this->vtime = vtime;
	this->vbrightness = vbrightness;
	this->vas = vas;
	this->vae = vae;

	this->time = 0;
	this->brightness = 0;
	this->as = point();
	this->ae = point();

	for (int i = 0; i < this->n; i++)
	{
		this->time += vtime[i];
		this->brightness += vbrightness[i];
		this->as = this->as + vas[i];
		this->ae = this->ae + vae[i];
	}
	this->time /= n;
	this->brightness /= n;
	this->as = this->as / n;
	this->ae = this->ae / n;
}

dataPoint::dataPoint(const dataPoint& dp)
{
	this->time = dp.time;
	this->brightness = dp.brightness;
	this->as = dp.as;
	this->ae = dp.ae;

	this->vtime = dp.vtime;
	this->vbrightness = dp.vbrightness;
	this->vas = dp.vas;
	this->vae = dp.vae;

	this->n = dp.n;
}

dataPoint& dataPoint::operator=(const dataPoint& dp)
{
	this->time = dp.time;
	this->brightness = dp.brightness;
	this->as = dp.as;
	this->ae = dp.ae;

	this->n = dp.n;

	this->vtime = dp.vtime;
	this->vbrightness = dp.vbrightness;
	this->vas = dp.vas;
	this->vae = dp.vae;

	return *this;
}

point dataPoint::PAB()
{
	return this->as.rn() + this->ae.rn();
}

double dataPoint::phaseAngle()
{
	return this->as.angleP(this->ae);
}