/* this computes brightness and its derivatives w.r.t. parameters */

/*
Copyright (C) 2006  Mikko Kaasalainen, Josef Durech
Modified by Shuai Feng in 2024

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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"
#include "constants.h"

double bright(double* ee, double* ee0, double t, double* cg,
	double* dyda, int ncoef, int Nphpar, int Numfac,
	double* Scale, double** Blmat, double*** Dblm,
	double** Nor, double* Area, double* Darea, double** Dg, double Phi_0)
{
	int ncoef0, i, j, k,
		* incl;

	double cos_alpha, br, cl, cls, alpha, sum, sum0, dnom,
		e[4], e0[4],
		* php, * dphp,
		* mu, * mu0, * s, * dbr,
		* dsmu, * dsmu0,
		de[4][5], de0[4][5], tmat[4][4],
		dtm[5][4][4];

	incl = vector_int(MAX_N_FAC + 1);

	php = vector_double(N_PHOT_PAR + 1);
	dphp = vector_double(N_PHOT_PAR + 1);
	mu = vector_double(MAX_N_FAC + 1);
	mu0 = vector_double(MAX_N_FAC + 1 + 1);
	s = vector_double(MAX_N_FAC + 1);
	dbr = vector_double(MAX_N_FAC + 1);
	dsmu = vector_double(MAX_N_FAC + 1);
	dsmu0 = vector_double(MAX_N_FAC + 1);

	ncoef0 = ncoef - 2 - Nphpar;
	cl = exp(cg[ncoef - 1]); /* Lambert */
	cls = cg[ncoef];       /* Lommel-Seeliger */
	cos_alpha = dot_product(ee, ee0);
	alpha = acos(cos_alpha);
	for (i = 1; i <= Nphpar; i++)
		php[i] = cg[ncoef0 + i];

	phasec(dphp, alpha, php, Scale); /* computes also Scale */
	//printf("%.7f %.7f\n", t, alpha);
	matrix(cg[ncoef0 - 1], cg[ncoef0], t, tmat, dtm, Blmat, Dblm, Phi_0);

	br = 0;
	/* Directions (and ders.) in the rotating system */
	for (i = 1; i <= 3; i++)
	{
		e[i] = 0;
		e0[i] = 0;
		for (j = 1; j <= 3; j++)
		{
			e[i] = e[i] + tmat[i][j] * ee[j];
			e0[i] = e0[i] + tmat[i][j] * ee0[j];
		}
		for (j = 1; j <= 4; j++)
		{
			de[i][j] = 0;
			de0[i][j] = 0;
			for (k = 1; k <= 3; k++)
			{
				de[i][j] = de[i][j] + dtm[j][i][k] * ee[k];
				de0[i][j] = de0[i][j] + dtm[j][i][k] * ee0[k];
			}
		}
	}

	/*Integrated brightness (phase coeff. used later) */
	for (i = 1; i <= Numfac; i++)
	{
		incl[i] = 0;
		mu[i] = e[1] * Nor[i][1] + e[2] * Nor[i][2] + e[3] * Nor[i][3];
		mu0[i] = e0[1] * Nor[i][1] + e0[2] * Nor[i][2] + e0[3] * Nor[i][3];
		if ((mu[i] > TINY) && (mu0[i] > TINY))
		{
			incl[i] = 1;
			dnom = mu[i] + mu0[i];
			s[i] = mu[i] * mu0[i] * (cl + cls / dnom);
			br = br + Area[i] * s[i];
			dsmu[i] = cls * pow(mu0[i] / dnom, 2) + cl * mu0[i];
			dsmu0[i] = cls * pow(mu[i] / dnom, 2) + cl * mu[i];
			dbr[i] = Darea[i] * s[i];
		}
	}

	/* Derivatives of brightness w.r.t. g-coeffs */
	for (i = 1; i <= ncoef0 - 4; i++)
	{
		dyda[i] = 0;
		for (j = 1; j <= Numfac; j++)
			if (incl[j] == 1)
				dyda[i] = dyda[i] + dbr[j] * Dg[j][i];
		dyda[i] = *Scale * dyda[i];
	}
	/*   printf("%f \n", dyda[1]);   */
	   /* Ders. of brightness w.r.t. rotation parameters */
	for (k = 1; k <= 4; k++)
	{
		dyda[ncoef0 - 4 + k] = 0;
		for (i = 1; i <= Numfac; i++)
			if (incl[i] == 1)
			{
				sum = 0;
				sum0 = 0;
				for (j = 1; j <= 3; j++)
				{
					sum = sum + Nor[i][j] * de[j][k];
					sum0 = sum0 + Nor[i][j] * de0[j][k];
				}
				dyda[ncoef0 - 4 + k] = dyda[ncoef0 - 4 + k] + Area[i] * (dsmu[i] * sum + dsmu0[i] * sum0);
			}
		dyda[ncoef0 - 4 + k] = *Scale * dyda[ncoef0 - 4 + k];
	}


	/* Ders. of br. w.r.t. phase function params. */
	for (i = 1; i <= Nphpar; i++)
		dyda[ncoef0 + i] = br * dphp[i];

	/* Ders. of br. w.r.t. cl, cls */
	dyda[ncoef - 1] = 0;
	dyda[ncoef] = 0;
	for (i = 1; i <= Numfac; i++)
		if (incl[i] == 1)
		{
			dyda[ncoef - 1] = dyda[ncoef - 1] + mu[i] * mu0[i] * Area[i];
			dyda[ncoef] = dyda[ncoef] + Area[i] * mu[i] * mu0[i] / (mu[i] + mu0[i]);
		}
	dyda[ncoef - 1] = *Scale * dyda[ncoef - 1] * cl;
	dyda[ncoef] = *Scale * dyda[ncoef];
	//printf("%f\n", dyda[ncoef0]);

	/* Scaled brightness */
	br = br * *Scale;

	deallocate_vector((void*)incl);
	deallocate_vector((void*)php);
	deallocate_vector((void*)dphp);
	deallocate_vector((void*)mu);
	deallocate_vector((void*)mu0);
	deallocate_vector((void*)s);
	deallocate_vector((void*)dbr);
	deallocate_vector((void*)dsmu);
	deallocate_vector((void*)dsmu0);

	return(br);
}
