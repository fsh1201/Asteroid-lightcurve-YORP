
void trifac(int nrows, int** ifp);
void areanorm(double* t, double* f, int ndir, int nfac, int** ifp,
	double* at, double* af,
	double** Nor, double* Darea);
void sphfunc(int ndir, double* at, double* af,
	int Mmax, int Lmax,
	double** Ts, double** Tc, double** Fs, double** Fc, double*** Pleg, double** Dsph);
void ellfit(double* cg, double a, double b, double c,
	int ndir, int ncoef, double* at, double* af,
	int Mmax, int Lmax,
	double*** Pleg);
void lubksb(double** a, int n, int indx[], double b[]);
void ludcmp(double** a, int n, int indx[], double d[]);
void mrqmin(double** x1, double** x2, double* x3, double* y,
	double* sig, double* a, int* ia, int ma,
	double** covar, double** alpha, double (*funcs)(), int Nphpar, int Deallocate, int Lastcall, int Numfac, int Mmax, int Lmax,
	int Lcurves, int* Inrel, int* Lpoints,
	double* Alamda, double* Scale, double** Blmat, double*** Dblm,
	double** Nor, double* Area, double* Darea, double** Dg, double* Chisq, double* Ochisq,
	double** Fc, double** Fs, double*** Pleg, double** Dsph, double* Yout, double* Sclnw,
	int Ncoef, double Phi_0);
double mrqcof(double** x1, double** x2, double* x3, double* y,
	double* sig, double* a, int* ia, int ma,
	double** alpha, double* beta, double (*funcs)(), int Nphpar, int Numfac, int Mmax, int Lmax,
	int Lcurves, int* Inrel, int* Lpoints, int Lastcall,
	double* Alamda, double* Scale, double** Blmat, double*** Dblm,
	double** Nor, double* Area, double* Darea, double** Dg,
	double** Fc, double** Fs, double*** Pleg, double** Dsph, double* Yout, double* Sclnw,
	int Ncoef, double Phi_0);
void curv(double cg[], int Numfac, int Mmax, int Lmax,
	double** Fc, double** Fs, double*** Pleg, double* Area, double* Darea, double** Dsph, double** Dg);
void blmatrix(double bet, double lam, double** Blmat, double*** Dblm);
double conv(int nc, double dres[], int ma, int Numfac, int Ncoef,
	double* Area, double** Nor, double* Darea, double** Dg);
void gauss(double** aa, int n, double b[]);
void covsrt(double** covar, int ma, int ia[], int mfit);
void phasec(double dcdp[], double alpha, double p[], double* Scale);
void matrix(double omg, double nu, double t, double tmat[][4], double dtm[][4][4],
	double** Blmat, double*** Dblm, double Phi_0);
double bright(double* ee, double* ee0, double t, double* cg,
	double* dyda, int ncoef, int Nphpar, int Numfac,
	double* Scale, double** Blmat, double*** Dblm,
	double** Nor, double* Area, double* Darea, double** Dg, double Phi_0);

double* vector_double(int length);
int* vector_int(int length);
double** matrix_double(int rows, int columns);
int** matrix_int(int rows, int columns);
double*** matrix_3_double(int n_1, int n_2, int n_3);
void deallocate_vector(void* p_x);
void deallocate_matrix(void** p_x, int rows);
void deallocate_matrix_3(void*** p_x, int n_1, int n_2);

double dot_product(double a[], double b[]);
