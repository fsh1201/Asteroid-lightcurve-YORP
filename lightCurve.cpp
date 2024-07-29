/* lightcurve library, a lightcurve consists a number of datapoints */

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

#include "lightCurve.h"


lightCurve::lightCurve(std::vector<dataPoint> lc)
{
	this->vdp = lc;
	this->length = lc.size();
}

lightCurve::lightCurve(double** lc, int n)
{
	for (int i = 0; i < n; i++)
	{
		this->vdp.push_back(dataPoint(lc[i][0], lc[i][1], point(lc[i][2], lc[i][3], lc[i][4]), point(lc[i][5], lc[i][6], lc[i][7])));
	}
	this->length = n;
}

lightCurve::lightCurve(const lightCurve& lc)
{
	this->length = lc.length;
	this->vdp = lc.vdp;
}

lightCurve::lightCurve()
{
	length = 0;
	this->vdp = std::vector<dataPoint>();
}

dataPoint& lightCurve::operator[](int k)
{
	return this->vdp[k];
}

const dataPoint& lightCurve::operator[](int k) const
{
	return this->vdp[k];
}

void lightCurve::normalize(int method)
{
	double scale = 0;
	double lsum = 0;
	double lmax = -1e26, lmin = 1e26;
	for (int i = 0; i < this->length; i++)
	{
		lsum += this->vdp[i].brightness;
		lmax = fmax(lmax, vdp[i].brightness);
		lmin = fmin(lmin, vdp[i].brightness);
	}
	if (method == normal::mean)
	{
		scale = lsum / length;
		for (int i = 0; i < this->length; i++)
		{
			this->vdp[i].brightness /= scale;
		}
	}
	if (method == normal::max_min)
	{
		scale = lmax - lmin;
		for (int i = 0; i < this->length; i++)
		{
			this->vdp[i].brightness = (this->vdp[i].brightness - lmin) / scale;
		}
	}
}

lightCurve& lightCurve::operator=(const lightCurve& lc)
{
	this->length = lc.length;
	this->vdp = lc.vdp;
	return *this;
}

double lightCurve::relativeChi2(lightCurve lc)
{
	normalize(normal::mean);
	lc.normalize(normal::mean);
	double chi2 = 0;
	for (int i = 0; i < vdp.size(); i++)
	{
		chi2 += pow(vdp[i].brightness - lc.vdp[i].brightness, 2);
	}
	return chi2;
}