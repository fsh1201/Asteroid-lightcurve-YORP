/* this is a modified version of the routine
   from Press, Teukolsky, Vetterling, and Flannery - Numerical Recipes in C, CUP 1992 */

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


#include <stdio.h>
#include "constants.h"
#include "declarations.h"

void mrqmin(double** x1, double** x2, double* x3, double* y,
	double* sig, double* a, int* ia, int ma,
	double** covar, double** alpha, double (*funcs)(), int Nphpar, int Deallocate, int Lastcall, int Numfac, int Mmax, int Lmax,
	int Lcurves, int* Inrel, int* Lpoints,
	double* Alamda, double* Scale, double** Blmat, double*** Dblm,
	double** Nor, double* Area, double* Darea, double** Dg, double* Chisq, double* Ochisq,
	double** Fc, double** Fs, double*** Pleg, double** Dsph, double* Yout, double* Sclnw,
	int Ncoef, double Phi_0)
{

	int j, k, l;
	static int mfit; /* it is set in the firs call */

	static double* atry, * beta, * da;

	double temp;

	/* dealocates memory when usen in period_search */
	if (Deallocate == 1)
	{

		deallocate_vector((void*)atry);
		deallocate_vector((void*)beta);
		deallocate_vector((void*)da);
		return;
	}

	if (Lastcall != 1)
	{
		if (*Alamda < 0)
		{
			atry = vector_double(ma);
			beta = vector_double(ma);
			da = vector_double(ma);

			/* number of fitted parameters */
			mfit = 0;
			for (j = 1; j <= ma; j++)
				if (ia[j]) mfit++;

			*Alamda = ALAMDA_START; /* initial alambda */

			temp = mrqcof(x1, x2, x3, y, sig, a, ia, ma, alpha, beta, funcs, Nphpar, Numfac, Mmax, Lmax,
				Lcurves, Inrel, Lpoints, Lastcall,
				Alamda, Scale, Blmat, Dblm,
				Nor, Area, Darea, Dg,
				Fc, Fs, Pleg, Dsph, Yout, Sclnw,
				Ncoef, Phi_0);    //chi2, calculate alpha[],beta[]

			*Ochisq = temp; //chi20
			for (j = 1; j <= ma; j++)
				atry[j] = a[j];
		}
		for (j = 1; j <= mfit; j++)
		{
			for (k = 1; k <= mfit; k++)
				covar[j][k] = alpha[j][k];
			covar[j][j] = alpha[j][j] * (1 + *Alamda);
			da[j] = beta[j];
		}

		gauss(covar, mfit, da);
		//printf("%g\n", da[Ncoef + 4]);

		if (*Alamda == 0)
		{
			covsrt(covar, ma, ia, mfit);
			return;
		}

		j = 0;
		for (l = 1; l <= ma; l++)
			if (ia[l])
			{
				j++;
				atry[l] = a[l] + da[j];
			}
	} /* Lastcall != 1 */

	if (Lastcall == 1)
		for (l = 1; l <= ma; l++)
			atry[l] = a[l];

	temp = mrqcof(x1, x2, x3, y, sig, atry, ia, ma, covar, da, funcs, Nphpar, Numfac, Mmax, Lmax,
		Lcurves, Inrel, Lpoints, Lastcall,
		Alamda, Scale, Blmat, Dblm,
		Nor, Area, Darea, Dg,
		Fc, Fs, Pleg, Dsph, Yout, Sclnw,
		Ncoef, Phi_0);
	*Chisq = temp;    //chi2


	if (Lastcall == 1)
	{
		deallocate_vector((void*)atry);
		deallocate_vector((void*)beta);
		deallocate_vector((void*)da);
		return;
	}

	if (temp < *Ochisq)
	{
		*Alamda = *Alamda / ALAMDA_COEFF;
		for (j = 1; j <= mfit; j++)
		{
			for (k = 1; k <= mfit; k++)
				alpha[j][k] = covar[j][k];
			beta[j] = da[j];
		}
		for (l = 1; l <= ma; l++)
			a[l] = atry[l];
	}
	else
	{
		*Alamda = ALAMDA_COEFF * *Alamda;
		*Chisq = *Ochisq;
	}

}

