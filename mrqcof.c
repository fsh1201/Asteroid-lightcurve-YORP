/* this is a modified veresion of the routine from
   Press, Teukolsky, Vetterling, and Flannery - Numerical Recipes in C, CUP 1992 */

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
#include <stdlib.h>
#include "declarations.h"
#include "constants.h"

double mrqcof(double** x1, double** x2, double* x3, double* y,
	double* sig, double* a, int* ia, int ma,
	double** alpha, double* beta, double (*funcs)(), int Nphpar, int Numfac, int Mmax, int Lmax,
	int Lcurves, int* Inrel, int* Lpoints, int Lastcall,
	double* Alamda, double* Scale, double** Blmat, double*** Dblm,
	double** Nor, double* Area, double* Darea, double** Dg,
	double** Fc, double** Fs, double*** Pleg, double** Dsph, double* Yout, double* Sclnw,
	int Ncoef, double Phi_0)
{

	int mfit, i, j, k, l, m, np, np1, np2, jp, ic;

	double xx1[4], xx2[4], dy, sig2i, wt, * dyda, ymod,
		* ytemp, ** dytemp,
		* dave,
		coef, ave = 0, trial_chisq;

	dyda = vector_double(MAX_N_PAR + 1);
	ytemp = vector_double(POINTS_MAX + 1);
	dytemp = matrix_double(POINTS_MAX + 1, MAX_N_PAR + 1);
	dave = vector_double(MAX_N_PAR + 1);


	/* N.B. curv and blmatrix called outside bright
	   because output same for all points */
	curv(a, Numfac, Mmax, Lmax,
		Fc, Fs, Pleg, Area, Darea, Dsph, Dg); //calculate g[],Area[],Dg[][]

	blmatrix(a[ma - 5 - Nphpar], a[ma - 4 - Nphpar], Blmat, Dblm); //calculate Blmat[][],Dblm[][][]

	mfit = 0;
	for (j = 1; j <= ma; j++)
		if (ia[j])
			mfit++;
	for (j = 1; j <= mfit; j++)
	{
		for (k = 1; k <= j; k++)
			alpha[j][k] = 0;
		beta[j] = 0;
	}
	trial_chisq = 0;
	np = 0;
	np1 = 0;
	np2 = 0;

	for (i = 1; i <= Lcurves; i++)
	{
		if (Inrel[i] == 1) /* is the LC relative? */
		{
			ave = 0;
			for (l = 1; l <= ma; l++)
				dave[l] = 0;
		}
		for (jp = 1; jp <= Lpoints[i]; jp++)
		{
			np++;
			for (ic = 1; ic <= 3; ic++) /* position vectors */
			{
				xx1[ic] = x1[np][ic];
				xx2[ic] = x2[np][ic];
			}

			if (i < Lcurves)
				ymod = funcs(xx1, xx2, x3[np], a, dyda, ma, Nphpar, Numfac,
					Scale, Blmat, Dblm,
					Nor, Area, Darea, Dg, Phi_0);	//calculate dyda[], brightness
			else
				ymod = conv(jp, dyda, ma, Numfac, ma,
					Area, Nor, Darea, Dg);

			ytemp[jp] = ymod;

			if (Inrel[i] == 1)
				ave = ave + ymod;

			for (l = 1; l <= ma; l++)
			{
				dytemp[jp][l] = dyda[l];
				if (Inrel[i] == 1)
					dave[l] = dave[l] + dyda[l];
			}
			/* save lightcurves */

			if (Lastcall == 1)
			{
				Yout[np] = ymod;
				/*FILE* scaleFile = fopen("/mnt/d/wuli/damit/files/asteroid_101/scale.txt", "a+");
				fprintf(scaleFile, "%f\n", *Scale);
				fclose(scaleFile);*/
			}
		} /* jp, lpoints */

		if (Lastcall != 1)
		{
			for (jp = 1; jp <= Lpoints[i]; jp++)
			{
				np1++;
				if (Inrel[i] == 1)
				{
					coef = sig[np1] * Lpoints[i] / ave;
					for (l = 1; l <= ma; l++)
						dytemp[jp][l] = coef * (dytemp[jp][l] - ytemp[jp] * dave[l] / ave);
					ytemp[jp] = coef * ytemp[jp];
					/* Set the size scale coeff. deriv. explicitly zero for relative lcurves */
					dytemp[jp][1] = 0;
				}
			}

			for (jp = 1; jp <= Lpoints[i]; jp++)
			{
				ymod = ytemp[jp];
				for (l = 1; l <= ma; l++)
					dyda[l] = dytemp[jp][l];
				np2++;
				sig2i = 1 / (sig[np2] * sig[np2]);
				dy = y[np2] - ymod;
				j = 0;
				for (l = 1; l <= ma; l++)
				{
					if (ia[l])
					{
						j++;
						wt = dyda[l] * sig2i;
						k = 0;
						for (m = 1; m <= l; m++)
						{
							if (ia[m])
							{
								k++;
								alpha[j][k] = alpha[j][k] + wt * dyda[m];	//alpha=JtJ
							}
						} /* m */
						beta[j] = beta[j] + dy * wt;
					}
				} /* l */
				trial_chisq = trial_chisq + dy * dy * sig2i;
			} /* jp */
		} /* Lastcall != 1 */

		if ((Lastcall == 1) && (Inrel[i] == 1))
			Sclnw[i] = *Scale * Lpoints[i] * sig[np] / ave;

	} /* i,  lcurves */

	for (j = 2; j <= mfit; j++)
		for (k = 1; k <= j - 1; k++)
			alpha[k][j] = alpha[j][k];


	deallocate_vector((void*)dyda);
	deallocate_vector((void*)ytemp);
	deallocate_vector((void*)dave);
	deallocate_matrix((void**)dytemp, POINTS_MAX + 1);

	return trial_chisq;

}

