#pragma once
#include "CommonDef.h"
//  mgmres.hpp
//
//  Thanks to Philip Sakievich for correcting mismatches between this HPP
//  file and the corresponding CPP file, 23 March 2016!
//

class GMRESLib
{
public:
	GMRESLib();
	~GMRESLib();

	static void atx_cr(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[]);
	static void atx_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[]);
	static void ax_cr(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[]);
	static void ax_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[]);
	static void diagonal_pointer_cr(int n, int nz_num, int ia[], int ja[], int ua[]);
	static void ilu_cr(int n, int nz_num, int ia[], int ja[], double a[], int ua[], double l[]);
	static void lus_cr(int n, int nz_num, int ia[], int ja[], double l[], int ua[], double r[], double z[]);
	static void mgmres_st(int n, int nz_num, int ia[], int ja[], double a[], double x[],
		double rhs[], int itr_max, int mr, double tol_abs, double tol_rel);
	static void mult_givens(double c, double s, int k, double g[]);
	static void pmgmres_ilu_cr(int n, int nz_num, int ia[], int ja[], double a[],
		double x[], double rhs[], int itr_max, int mr, double tol_abs,
		double tol_rel);
	static double r8vec_dot(int n, double a1[], double a2[]);
	static double *r8vec_uniform_01(int n, int &seed);
	static void rearrange_cr(int n, int nz_num, int ia[], int ja[], double a[]);

	static void test01();
	static void test02();
	static void test03();
	static void test04();

private:

};

GMRESLib::GMRESLib()
{
}

GMRESLib::~GMRESLib()
{
}

//****************************************************************************80

inline void GMRESLib::atx_cr(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[])

//****************************************************************************80
//
//  Purpose:
//
//    ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
//
//  Discussion:
//
//    The Sparse Compressed Row storage format is used.
//
//    The matrix A is assumed to be sparse.  To save on storage, only
//    the nonzero entries of A are stored.  The vector JA stores the
//    column index of the nonzero value.  The nonzero values are sorted
//    by row, and the compressed row vector IA then has the property that
//    the entries in A and JA that correspond to row I occur in indices
//    IA[I] through IA[I+1]-1.
//
//    For this version of MGMRES, the row and column indices are assumed
//    to use the C/C++ convention, in which indexing begins at 0.
//
//    If your index vectors IA and JA are set up so that indexing is based 
//    at 1, then each use of those vectors should be shifted down by 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2007
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
//    Charles Romine, Henk van der Vorst,
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Methods,
//    SIAM, 1994,
//    ISBN: 0898714710,
//    LC: QA297.8.T45.
//
//    Tim Kelley,
//    Iterative Methods for Linear and Nonlinear Equations,
//    SIAM, 2004,
//    ISBN: 0898713528,
//    LC: QA297.8.K45.
//
//    Yousef Saad,
//    Iterative Methods for Sparse Linear Systems,
//    Second Edition,
//    SIAM, 2003,
//    ISBN: 0898715342,
//    LC: QA188.S17.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, int NZ_NUM, the number of nonzeros.
//
//    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
//    of the matrix values.  The row vector has been compressed.
//
//    Input, double A[NZ_NUM], the matrix values.
//
//    Input, double X[N], the vector to be multiplied by A'.
//
//    Output, double W[N], the value of A'*X.
//
{
	int i;
	int k;
	int k1;
	int k2;

	for (i = 0; i < n; i++)
	{
		w[i] = 0.0;
		k1 = ia[i];
		k2 = ia[i + 1];
		for (k = k1; k < k2; k++)
		{
			w[ja[k]] = w[ja[k]] + a[k] * x[i];
		}
	}
	return;
}
//****************************************************************************80

inline void GMRESLib::atx_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[])

//****************************************************************************80
//
//  Purpose:
//
//    ATX_ST computes A'*x for a matrix stored in sparse triplet form.
//
//  Discussion:
//
//    The matrix A is assumed to be sparse.  To save on storage, only
//    the nonzero entries of A are stored.  For instance, the K-th nonzero
//    entry in the matrix is stored by:
//
//      A(K) = value of entry,
//      IA(K) = row of entry,
//      JA(K) = column of entry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2007
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
//    Charles Romine, Henk van der Vorst,
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Methods,
//    SIAM, 1994,
//    ISBN: 0898714710,
//    LC: QA297.8.T45.
//
//    Tim Kelley,
//    Iterative Methods for Linear and Nonlinear Equations,
//    SIAM, 2004,
//    ISBN: 0898713528,
//    LC: QA297.8.K45.
//
//    Yousef Saad,
//    Iterative Methods for Sparse Linear Systems,
//    Second Edition,
//    SIAM, 2003,
//    ISBN: 0898715342,
//    LC: QA188.S17.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, int NZ_NUM, the number of nonzeros.
//
//    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
//    of the matrix values.
//
//    Input, double A[NZ_NUM], the matrix values.
//
//    Input, double X[N], the vector to be multiplied by A'.
//
//    Output, double W[N], the value of A'*X.
//
{
	//int i;
	//int k;
#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		w[i] = 0.0;
	}

	int i;
	int j;
#pragma omp parallel for private(i,j)
	for (int k = 0; k < nz_num; k++)
	{
		i = ia[k];
		j = ja[k];
		w[j] = w[j] + a[k] * x[i];
	}
	return;
}
//****************************************************************************80

inline void GMRESLib::ax_cr(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[])

//****************************************************************************80
//
//  Purpose:
//
//    AX_CR computes A*x for a matrix stored in sparse compressed row form.
//
//  Discussion:
//
//    The Sparse Compressed Row storage format is used.
//
//    The matrix A is assumed to be sparse.  To save on storage, only
//    the nonzero entries of A are stored.  The vector JA stores the
//    column index of the nonzero value.  The nonzero values are sorted
//    by row, and the compressed row vector IA then has the property that
//    the entries in A and JA that correspond to row I occur in indices
//    IA[I] through IA[I+1]-1.
//
//    For this version of MGMRES, the row and column indices are assumed
//    to use the C/C++ convention, in which indexing begins at 0.
//
//    If your index vectors IA and JA are set up so that indexing is based 
//    at 1, then each use of those vectors should be shifted down by 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2007
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
//    Charles Romine, Henk van der Vorst,
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Methods,
//    SIAM, 1994,
//    ISBN: 0898714710,
//    LC: QA297.8.T45.
//
//    Tim Kelley,
//    Iterative Methods for Linear and Nonlinear Equations,
//    SIAM, 2004,
//    ISBN: 0898713528,
//    LC: QA297.8.K45.
//
//    Yousef Saad,
//    Iterative Methods for Sparse Linear Systems,
//    Second Edition,
//    SIAM, 2003,
//    ISBN: 0898715342,
//    LC: QA188.S17.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, int NZ_NUM, the number of nonzeros.
//
//    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
//    of the matrix values.  The row vector has been compressed.
//
//    Input, double A[NZ_NUM], the matrix values.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double W[N], the value of A*X.
//
{
	//int i;
	//int k;
	int k1;
	int k2;
#pragma omp parallel for private(k1,k2)
	for (int i = 0; i < n; i++)
	{
		w[i] = 0.0;
		k1 = ia[i];
		k2 = ia[i + 1];
		for (int k = k1; k < k2; k++)
		{
			w[i] = w[i] + a[k] * x[ja[k]];
		}
	}
	return;
}
//****************************************************************************80

inline void GMRESLib::ax_st(int n, int nz_num, int ia[], int ja[], double a[], double x[], double w[])

//****************************************************************************80
//
//  Purpose:
//
//    AX_ST computes A*x for a matrix stored in sparse triplet form.
//
//  Discussion:
//
//    The matrix A is assumed to be sparse.  To save on storage, only
//    the nonzero entries of A are stored.  For instance, the K-th nonzero
//    entry in the matrix is stored by:
//
//      A(K) = value of entry,
//      IA(K) = row of entry,
//      JA(K) = column of entry.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2007
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
//    Charles Romine, Henk van der Vorst,
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Methods,
//    SIAM, 1994,
//    ISBN: 0898714710,
//    LC: QA297.8.T45.
//
//    Tim Kelley,
//    Iterative Methods for Linear and Nonlinear Equations,
//    SIAM, 2004,
//    ISBN: 0898713528,
//    LC: QA297.8.K45.
//
//    Yousef Saad,
//    Iterative Methods for Sparse Linear Systems,
//    Second Edition,
//    SIAM, 2003,
//    ISBN: 0898715342,
//    LC: QA188.S17.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, int NZ_NUM, the number of nonzeros.
//
//    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
//    of the matrix values.
//
//    Input, double A[NZ_NUM], the matrix values.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double W[N], the value of A*X.
//
{

#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		w[i] = 0.0;
	}

	int i;
	int j;
#pragma omp parallel for private(i,j)
	for (int k = 0; k < nz_num; k++)
	{
		i = ia[k];
		j = ja[k];
		w[i] = w[i] + a[k] * x[j];
	}
	return;
}
//****************************************************************************80

inline void GMRESLib::diagonal_pointer_cr(int n, int nz_num, int ia[], int ja[], int ua[])

//****************************************************************************80
//
//  Purpose:
//
//    DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
//
//  Discussion:
//
//    The matrix A is assumed to be stored in compressed row format.  Only
//    the nonzero entries of A are stored.  The vector JA stores the
//    column index of the nonzero value.  The nonzero values are sorted
//    by row, and the compressed row vector IA then has the property that
//    the entries in A and JA that correspond to row I occur in indices
//    IA[I] through IA[I+1]-1.
//
//    The array UA can be used to locate the diagonal elements of the matrix.
//
//    It is assumed that every row of the matrix includes a diagonal element,
//    and that the elements of each row have been ascending sorted.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2007
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, int NZ_NUM, the number of nonzeros.
//
//    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
//    of the matrix values.  The row vector has been compressed.  On output,
//    the order of the entries of JA may have changed because of the sorting.
//
//    Output, int UA[N], the index of the diagonal element of each row.
//
{
	//int i;
	//int j;
	int j1;
	int j2;
	int k;

	for (int i = 0; i < n; i++)
	{
		ua[i] = -1;
		j1 = ia[i];
		j2 = ia[i + 1];

		for (int j = j1; j < j2; j++)
		{
			if (ja[j] == i)
			{
				ua[i] = j;
			}
		}

	}
	return;
}
//****************************************************************************80

inline void GMRESLib::ilu_cr(int n, int nz_num, int ia[], int ja[], double a[], int ua[], double l[])

//****************************************************************************80
//
//  Purpose:
//
//    ILU_CR computes the incomplete LU factorization of a matrix.
//
//  Discussion:
//
//    The matrix A is assumed to be stored in compressed row format.  Only
//    the nonzero entries of A are stored.  The vector JA stores the
//    column index of the nonzero value.  The nonzero values are sorted
//    by row, and the compressed row vector IA then has the property that
//    the entries in A and JA that correspond to row I occur in indices
//    IA[I] through IA[I+1]-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 July 2007
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, int NZ_NUM, the number of nonzeros.
//
//    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
//    of the matrix values.  The row vector has been compressed.
//
//    Input, double A[NZ_NUM], the matrix values.
//
//    Input, int UA[N], the index of the diagonal element of each row.
//
//    Output, double L[NZ_NUM], the ILU factorization of A.
//
{
	int *iw;
	//int i;
	//int j;
	//int jj;
	int jrow;
	int jw;
	//int k;
	double tl;

	iw = new int[n];
	//
	//  Copy A.
	//
#pragma omp parallel for
	for (int k = 0; k < nz_num; k++)
	{
		l[k] = a[k];
	}

	for (int i = 0; i < n; i++)
	{
		//
		//  IW points to the nonzero entries in row I.
		//
#pragma omp parallel for
		for (int j = 0; j < n; j++)
		{
			iw[j] = -1;
		}
#pragma omp parallel for
		for (int k = ia[i]; k <= ia[i + 1] - 1; k++)
		{
			iw[ja[k]] = k;
		}

		int j = ia[i];
		do
		{
			jrow = ja[j];
			if (i <= jrow)
			{
				break;
			}
			tl = l[j] * l[ua[jrow]];
			l[j] = tl;
			for (int jj = ua[jrow] + 1; jj <= ia[jrow + 1] - 1; jj++)
			{
				jw = iw[ja[jj]];
				if (jw != -1)
				{
					l[jw] = l[jw] - tl * l[jj];
				}
			}
			j = j + 1;
		} while (j <= ia[i + 1] - 1);

		ua[i] = j;

		if (jrow != i)
		{
			cerr << "\n";
			cerr << "ILU_CR - Fatal error!\n";
			cerr << "  JROW != I\n";
			cerr << "  JROW = " << jrow << "\n";
			cerr << "  I    = " << i << "\n";
			exit(1);
		}

		if (l[j] == 0.0)
		{
			cerr << "\n";
			cerr << "ILU_CR - Fatal error!\n";
			cerr << "  Zero pivot on step I = " << i << "\n";
			cerr << "  L[" << j << "] = 0.0\n";
			exit(1);
		}

		l[j] = 1.0 / l[j];
	}
#pragma omp parallel for
	for (int k = 0; k < n; k++)
	{
		l[ua[k]] = 1.0 / l[ua[k]];
	}
	//
	//  Free memory.
	//
	delete[] iw;

	return;
}
//****************************************************************************80

inline void GMRESLib::lus_cr(int n, int nz_num, int ia[], int ja[], double l[], int ua[], double r[], double z[])

//****************************************************************************80
//
//  Purpose:
//
//    LUS_CR applies the incomplete LU preconditioner.
//
//  Discussion:
//
//    The linear system M * Z = R is solved for Z.  M is the incomplete
//    LU preconditioner matrix, and R is a vector supplied by the user.
//    So essentially, we're solving L * U * Z = R.
//
//    The matrix A is assumed to be stored in compressed row format.  Only
//    the nonzero entries of A are stored.  The vector JA stores the
//    column index of the nonzero value.  The nonzero values are sorted
//    by row, and the compressed row vector IA then has the property that
//    the entries in A and JA that correspond to row I occur in indices
//    IA[I] through IA[I+1]-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2007
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, int NZ_NUM, the number of nonzeros.
//
//    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
//    of the matrix values.  The row vector has been compressed.
//
//    Input, double L[NZ_NUM], the matrix values.
//
//    Input, int UA[N], the index of the diagonal element of each row.
//
//    Input, double R[N], the right hand side.
//
//    Output, double Z[N], the solution of the system M * Z = R.
//
{
	//int i;
	//int j;
	double *w;

	w = new double[n];
	//
	//  Copy R in.
	//
#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		w[i] = r[i];
	}
	//
	//  Solve L * w = w where L is unit lower triangular.
	//
	for (int i = 1; i < n; i++)
	{
		for (int j = ia[i]; j < ua[i]; j++)
		{
			w[i] = w[i] - l[j] * w[ja[j]];
		}
	}
	//
	//  Solve U * w = w, where U is upper triangular.
	//
	for (int i = n - 1; 0 <= i; i--)
	{
		for (int j = ua[i] + 1; j < ia[i + 1]; j++)
		{
			w[i] = w[i] - l[j] * w[ja[j]];
		}
		w[i] = w[i] / l[ua[i]];
	}
	//
	//  Copy Z out.
	//
#pragma omp parallel for
	for (int i = 0; i < n; i++)
	{
		z[i] = w[i];
	}
	//
	//  Free memory.
	//
	delete[] w;

	return;
}
//****************************************************************************80

inline void GMRESLib::mgmres_st(int n, int nz_num, int ia[], int ja[], double a[], double x[],
	double rhs[], int itr_max, int mr, double tol_abs, double tol_rel)

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    MGMRES_ST applies restarted GMRES to a matrix in sparse triplet form.
	//
	//  Discussion:
	//
	//    The linear system A*X=B is solved iteratively.
	//
	//    The matrix A is assumed to be stored in sparse triplet form.  Only
	//    the nonzero entries of A are stored.  For instance, the K-th nonzero
	//    entry in the matrix is stored by:
	//
	//      A(K) = value of entry,
	//      IA(K) = row of entry,
	//      JA(K) = column of entry.
	//
	//    The "matrices" H and V are treated as one-dimensional vectors
	//    which store the matrix data in row major form.
	//
	//    This requires that references to H[I][J] be replaced by references
	//    to H[I+J*(MR+1)] and references to V[I][J] by V[I+J*N].
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license. 
	//
	//  Modified:
	//
	//    18 July 2007
	//
	//  Author:
	//
	//    Original C version by Lili Ju.
	//    C++ version by John Burkardt.
	//
	//  Reference:
	//
	//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
	//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
	//    Charles Romine, Henk van der Vorst,
	//    Templates for the Solution of Linear Systems:
	//    Building Blocks for Iterative Methods,
	//    SIAM, 1994,
	//    ISBN: 0898714710,
	//    LC: QA297.8.T45.
	//
	//    Tim Kelley,
	//    Iterative Methods for Linear and Nonlinear Equations,
	//    SIAM, 2004,
	//    ISBN: 0898713528,
	//    LC: QA297.8.K45.
	//
	//    Yousef Saad,
	//    Iterative Methods for Sparse Linear Systems,
	//    Second Edition,
	//    SIAM, 2003,
	//    ISBN: 0898715342,
	//    LC: QA188.S17.
	//
	//  Parameters:
	//
	//    Input, int N, the order of the linear system.
	//
	//    Input, int NZ_NUM, the number of nonzero matrix values.
	//
	//    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
	//    of the matrix values.
	//
	//    Input, double A[NZ_NUM], the matrix values.
	//
	//    Input/output, double X[N]; on input, an approximation to
	//    the solution.  On output, an improved approximation.
	//
	//    Input, double RHS[N], the right hand side of the linear system.
	//
	//    Input, int ITR_MAX, the maximum number of (outer) iterations to take.
	//
	//    Input, int MR, the maximum number of (inner) iterations to take.
	//    MR must be less than N.
	//
	//    Input, double TOL_ABS, an absolute tolerance applied to the
	//    current residual.
	//
	//    Input, double TOL_REL, a relative tolerance comparing the
	//    current residual to the initial residual.
	//
{
	double av;
	double *c;
	double delta = 1.0e-03;
	double *g;
	double *h;
	double htmp;
	//int i;
	int itr;
	int itr_used;
	//int j;
	int k;
	int k_copy;
	double mu;
	double *r;
	double rho;
	double rho_tol;
	double *s;
	double *v;
	bool verbose = true;
	double *y;

	c = new double[mr];
	g = new double[mr + 1];
	h = new double[(mr + 1)*mr];
	r = new double[n];
	s = new double[mr];
	v = new double[n*(mr + 1)];
	y = new double[mr + 1];

	itr_used = 0;

	if (n < mr)
	{
		cerr << "\n";
		cerr << "MGMRES_ST - Fatal error!\n";
		cerr << "  N < MR.\n";
		cerr << "  N = " << n << "\n";
		cerr << "  MR = " << mr << "\n";
		exit(1);
	}

	for (itr = 1; itr <= itr_max; itr++)
	{
		ax_st(n, nz_num, ia, ja, a, x, r);

#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			r[i] = rhs[i] - r[i];
		}

		rho = sqrt(r8vec_dot(n, r, r));

		if (verbose)
		{
			cout << "  ITR = " << itr << "  Residual = " << rho << "\n";
		}

		if (itr == 1)
		{
			rho_tol = rho * tol_rel;
		}
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			v[i + 0 * n] = r[i] / rho;
		}

		g[0] = rho;
#pragma omp parallel for
		for (int i = 1; i <= mr; i++)
		{
			g[i] = 0.0;
		}
#pragma omp parallel for
		for (int i = 0; i < mr + 1; i++)
		{
			for (int j = 0; j < mr; j++)
			{
				h[i + j*(mr + 1)] = 0.0;
			}
		}

		for (k = 1; k <= mr; k++)
		{
			k_copy = k;

			ax_st(n, nz_num, ia, ja, a, v + (k - 1)*n, v + k*n);

			av = sqrt(r8vec_dot(n, v + k*n, v + k*n));

			for (int j = 1; j <= k; j++)
			{
				h[(j - 1) + (k - 1)*(mr + 1)] = r8vec_dot(n, v + k*n, v + (j - 1)*n);
				for (int i = 0; i < n; i++)
				{
					v[i + k*n] = v[i + k*n] - h[(j - 1) + (k - 1)*(mr + 1)] * v[i + (j - 1)*n];
				}
			}

			h[k + (k - 1)*(mr + 1)] = sqrt(r8vec_dot(n, v + k*n, v + k*n));

			if ((av + delta * h[k + (k - 1)*(mr + 1)]) == av)
			{
				for (int j = 1; j <= k; j++)
				{
					htmp = r8vec_dot(n, v + k*n, v + (j - 1)*n);
					h[(j - 1) + (k - 1)*(mr + 1)] = h[(j - 1) + (k - 1)*(mr + 1)] + htmp;
					for (int i = 0; i < n; i++)
					{
						v[i + k*n] = v[i + k*n] - htmp * v[i + (j - 1)*n];
					}
				}
				h[k + (k - 1)*(mr + 1)] = sqrt(r8vec_dot(n, v + k*n, v + k*n));
			}

			if (h[k + (k - 1)*(mr + 1)] != 0.0)
			{
#pragma omp parallel for
				for (int i = 0; i < n; i++)
				{
					v[i + k*n] = v[i + k*n] / h[k + (k - 1)*(mr + 1)];
				}
			}

			if (1 < k)
			{
#pragma omp parallel for
				for (int i = 1; i <= k + 1; i++)
				{
					y[i - 1] = h[(i - 1) + (k - 1)*(mr + 1)];
				}
				for (int j = 1; j <= k - 1; j++)
				{
					mult_givens(c[j - 1], s[j - 1], j - 1, y);
				}
#pragma omp parallel for
				for (int i = 1; i <= k + 1; i++)
				{
					h[i - 1 + (k - 1)*(mr + 1)] = y[i - 1];
				}
			}
			mu = sqrt(pow(h[(k - 1) + (k - 1)*(mr + 1)], 2)
				+ pow(h[k + (k - 1)*(mr + 1)], 2));
			c[k - 1] = h[(k - 1) + (k - 1)*(mr + 1)] / mu;
			s[k - 1] = -h[k + (k - 1)*(mr + 1)] / mu;
			h[(k - 1) + (k - 1)*(mr + 1)] = c[k - 1] * h[(k - 1) + (k - 1)*(mr + 1)]
				- s[k - 1] * h[k + (k - 1)*(mr + 1)];
			h[k + (k - 1)*(mr + 1)] = 0;
			mult_givens(c[k - 1], s[k - 1], k - 1, g);

			rho = fabs(g[k]);

			itr_used = itr_used + 1;

			if (verbose && k == mr)
			{
				cout << "  K =   " << k << "  Residual = " << rho << "\n";
			}

			if (rho <= rho_tol && rho <= tol_abs)
			{
				break;
			}
		}

		k = k_copy - 1;
		y[k] = g[k] / h[k + k*(mr + 1)];

		for (int i = k; 1 <= i; i--)
		{
			y[i - 1] = g[i - 1];
			for (int j = i + 1; j <= k + 1; j++)
			{
				y[i - 1] = y[i - 1] - h[(i - 1) + (j - 1)*(mr + 1)] * y[j - 1];
			}
			y[i - 1] = y[i - 1] / h[(i - 1) + (i - 1)*(mr + 1)];
		}
#pragma omp parallel for
		for (int i = 1; i <= n; i++)
		{
			for (int j = 1; j <= k + 1; j++)
			{
				x[i - 1] = x[i - 1] + v[(i - 1) + (j - 1)*n] * y[j - 1];
			}
		}

		if (rho <= rho_tol && rho <= tol_abs)
		{
			break;
		}
	}

	if (verbose)
	{
		cout << "\n";
		cout << "MGMRES_ST\n";
		cout << "  Number of iterations = " << itr_used << "\n";
		cout << "  Final residual = " << rho << "\n";
	}
	//
	//  Free memory.
	//
	delete[] c;
	delete[] g;
	delete[] h;
	delete[] r;
	delete[] s;
	delete[] v;
	delete[] y;

	return;
}
//****************************************************************************80

inline void GMRESLib::mult_givens(double c, double s, int k, double g[])

//****************************************************************************80
//
//  Purpose:
//
//    MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 August 2006
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
//    Charles Romine, Henk van der Vorst,
//    Templates for the Solution of Linear Systems:
//    Building Blocks for Iterative Methods,
//    SIAM, 1994,
//    ISBN: 0898714710,
//    LC: QA297.8.T45.
//
//    Tim Kelley,
//    Iterative Methods for Linear and Nonlinear Equations,
//    SIAM, 2004,
//    ISBN: 0898713528,
//    LC: QA297.8.K45.
//
//    Yousef Saad,
//    Iterative Methods for Sparse Linear Systems,
//    Second Edition,
//    SIAM, 2003,
//    ISBN: 0898715342,
//    LC: QA188.S17.
//
//  Parameters:
//
//    Input, double C, S, the cosine and sine of a Givens
//    rotation.
//
//    Input, int K, indicates the location of the first vector entry.
//
//    Input/output, double G[K+2], the vector to be modified.  On output,
//    the Givens rotation has been applied to entries G(K) and G(K+1).
//
{
	double g1;
	double g2;

	g1 = c * g[k] - s * g[k + 1];
	g2 = s * g[k] + c * g[k + 1];

	g[k] = g1;
	g[k + 1] = g2;

	return;
}
//****************************************************************************80

inline void GMRESLib::pmgmres_ilu_cr(int n, int nz_num, int ia[], int ja[], double a[],
	double x[], double rhs[], int itr_max, int mr, double tol_abs,
	double tol_rel)

	//****************************************************************************80
	//
	//  Purpose:
	//
	//    PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
	//
	//  Discussion:
	//
	//    The matrix A is assumed to be stored in compressed row format.  Only
	//    the nonzero entries of A are stored.  The vector JA stores the
	//    column index of the nonzero value.  The nonzero values are sorted
	//    by row, and the compressed row vector IA then has the property that
	//    the entries in A and JA that correspond to row I occur in indices
	//    IA[I] through IA[I+1]-1.
	//
	//    This routine uses the incomplete LU decomposition for the
	//    preconditioning.  This preconditioner requires that the sparse
	//    matrix data structure supplies a storage position for each diagonal
	//    element of the matrix A, and that each diagonal element of the
	//    matrix A is not zero.
	//
	//    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
	//    corrections to the code on 31 May 2007.
	//
	//    This implementation of the code stores the doubly-dimensioned arrays
	//    H and V as vectors.  However, it follows the C convention of storing
	//    them by rows, rather than my own preference for storing them by
	//    columns.   I may come back and change this some time.
	//
	//  Licensing:
	//
	//    This code is distributed under the GNU LGPL license. 
	//
	//  Modified:
	//
	//    26 July 2007
	//
	//  Author:
	//
	//    Original C version by Lili Ju.
	//    C++ version by John Burkardt.
	//
	//  Reference:
	//
	//    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
	//    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
	//    Charles Romine, Henk van der Vorst,
	//    Templates for the Solution of Linear Systems:
	//    Building Blocks for Iterative Methods,
	//    SIAM, 1994.
	//    ISBN: 0898714710,
	//    LC: QA297.8.T45.
	//
	//    Tim Kelley,
	//    Iterative Methods for Linear and Nonlinear Equations,
	//    SIAM, 2004,
	//    ISBN: 0898713528,
	//    LC: QA297.8.K45.
	//
	//    Yousef Saad,
	//    Iterative Methods for Sparse Linear Systems,
	//    Second Edition,
	//    SIAM, 2003,
	//    ISBN: 0898715342,
	//    LC: QA188.S17.
	//
	//  Parameters:
	//
	//    Input, int N, the order of the linear system.
	//
	//    Input, int NZ_NUM, the number of nonzero matrix values.
	//
	//    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
	//    of the matrix values.  The row vector has been compressed.
	//
	//    Input, double A[NZ_NUM], the matrix values.
	//
	//    Input/output, double X[N]; on input, an approximation to
	//    the solution.  On output, an improved approximation.
	//
	//    Input, double RHS[N], the right hand side of the linear system.
	//
	//    Input, int ITR_MAX, the maximum number of (outer) iterations to take.
	//
	//    Input, int MR, the maximum number of (inner) iterations to take.
	//    MR must be less than N.
	//
	//    Input, double TOL_ABS, an absolute tolerance applied to the
	//    current residual.
	//
	//    Input, double TOL_REL, a relative tolerance comparing the
	//    current residual to the initial residual.
	//
{
	double av;
	double *c;
	double delta = 1.0e-03;
	double *g;
	double *h;
	double htmp;
	//int i;
	int itr;
	int itr_used;
	//int j;
	int k;
	int k_copy;
	double *l;
	double mu;
	double *r;
	double rho;
	double rho_tol;
	double *s;
	int *ua;
	double *v;
	int verbose = 1;
	double *y;

	itr_used = 0;

	c = new double[mr + 1];
	g = new double[mr + 1];
	h = new double[(mr + 1)*mr];
	l = new double[ia[n] + 1];
	r = new double[n];
	s = new double[mr + 1];
	ua = new int[n];
	v = new double[(mr + 1)*n];
	y = new double[mr + 1];

	rearrange_cr(n, nz_num, ia, ja, a);

	diagonal_pointer_cr(n, nz_num, ia, ja, ua);

	ilu_cr(n, nz_num, ia, ja, a, ua, l);

	if (verbose)
	{
		cout << "\n";
		cout << "PMGMRES_ILU_CR\n";
		cout << "  Number of unknowns = " << n << "\n";
	}

	for (itr = 0; itr < itr_max; itr++)
	{
		ax_cr(n, nz_num, ia, ja, a, x, r);

#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			r[i] = rhs[i] - r[i];
		}

		lus_cr(n, nz_num, ia, ja, l, ua, r, r);

		rho = sqrt(r8vec_dot(n, r, r));

		if (verbose)
		{
			cout << "  ITR = " << itr << "  Residual = " << rho << "\n";
		}

		if (itr == 0)
		{
			rho_tol = rho * tol_rel;
		}
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			v[0 * n + i] = r[i] / rho;
		}

		g[0] = rho;
#pragma omp parallel for
		for (int i = 1; i < mr + 1; i++)
		{
			g[i] = 0.0;
		}
#pragma omp parallel for
		for (int i = 0; i < mr + 1; i++)
		{
			for (int j = 0; j < mr; j++)
			{
				h[i*(mr)+j] = 0.0;
			}
		}

		for (k = 0; k < mr; k++)
		{
			k_copy = k;

			ax_cr(n, nz_num, ia, ja, a, v + k*n, v + (k + 1)*n);

			lus_cr(n, nz_num, ia, ja, l, ua, v + (k + 1)*n, v + (k + 1)*n);

			av = sqrt(r8vec_dot(n, v + (k + 1)*n, v + (k + 1)*n));

			for (int j = 0; j <= k; j++)
			{
				h[j*mr + k] = r8vec_dot(n, v + (k + 1)*n, v + j*n);
#pragma omp parallel for
				for (int i = 0; i < n; i++)
				{
					v[(k + 1)*n + i] = v[(k + 1)*n + i] - h[j*mr + k] * v[j*n + i];
				}
			}
			h[(k + 1)*mr + k] = sqrt(r8vec_dot(n, v + (k + 1)*n, v + (k + 1)*n));

			if ((av + delta * h[(k + 1)*mr + k]) == av)
			{
#pragma omp parallel for private(htmp)
				for (int j = 0; j < k + 1; j++)
				{
					htmp = r8vec_dot(n, v + (k + 1)*n, v + j*n);
					h[j*mr + k] = h[j*mr + k] + htmp;
					for (int i = 0; i < n; i++)
					{
						v[(k + 1)*n + i] = v[(k + 1)*n + i] - htmp * v[j*n + i];
					}
				}
				h[(k + 1)*mr + k] = sqrt(r8vec_dot(n, v + (k + 1)*n, v + (k + 1)*n));
			}

			if (h[(k + 1)*mr + k] != 0.0)
			{
#pragma omp parallel for
				for (int i = 0; i < n; i++)
				{
					v[(k + 1)*n + i] = v[(k + 1)*n + i] / h[(k + 1)*mr + k];
				}
			}

			if (0 < k)
			{
				for (int i = 0; i < k + 2; i++)
				{
					y[i] = h[i*mr + k];
				}
				for (int j = 0; j < k; j++)
				{
					mult_givens(c[j], s[j], j, y);
				}
				for (int i = 0; i < k + 2; i++)
				{
					h[i*mr + k] = y[i];
				}
			}
			mu = sqrt(h[k*mr + k] * h[k*mr + k] + h[(k + 1)*mr + k] * h[(k + 1)*mr + k]);
			c[k] = h[k*mr + k] / mu;
			s[k] = -h[(k + 1)*mr + k] / mu;
			h[k*mr + k] = c[k] * h[k*mr + k] - s[k] * h[(k + 1)*mr + k];
			h[(k + 1)*mr + k] = 0.0;
			mult_givens(c[k], s[k], k, g);

			rho = fabs(g[k + 1]);

			itr_used = itr_used + 1;

			if (verbose)
			{
				cout << "  K   = " << k << "  Residual = " << rho << "\n";
			}

			if (rho <= rho_tol && rho <= tol_abs)
			{
				break;
			}
		}

		k = k_copy;

		y[k] = g[k] / h[k*mr + k];
		for (int i = k - 1; 0 <= i; i--)
		{
			y[i] = g[i];
			for (int j = i + 1; j < k + 1; j++)
			{
				y[i] = y[i] - h[i*mr + j] * y[j];
			}
			y[i] = y[i] / h[i*mr + i];
		}
#pragma omp parallel for
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < k + 1; j++)
			{
				x[i] = x[i] + v[j*n + i] * y[j];
			}
		}

		if (rho <= rho_tol && rho <= tol_abs)
		{
			break;
		}
	}

	if (verbose)
	{
		cout << "\n";;
		cout << "PMGMRES_ILU_CR:\n";
		cout << "  Iterations = " << itr_used << "\n";
		cout << "  Final residual = " << rho << "\n";
	}
	//
	//  Free memory.
	//
	delete[] c;
	delete[] g;
	delete[] h;
	delete[] l;
	delete[] r;
	delete[] s;
	delete[] ua;
	delete[] v;
	delete[] y;

	return;
}
//****************************************************************************80

inline double GMRESLib::r8vec_dot(int n, double a1[], double a2[])

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_DOT computes the dot product of a pair of R8VEC's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT, the dot product of the vectors.
//
{
	int i;
	double value;

	value = 0.0;
#pragma omp parallel for reduction(+:value)
	for (i = 0; i < n; i++)
	{
		value += a1[i] * a2[i];
	}
	return value;
}
//****************************************************************************80

inline double* GMRESLib::r8vec_uniform_01(int n, int &seed)

//****************************************************************************80
//
//  Purpose:
//
//    R8VEC_UNIFORM_01 returns a unit pseudorandom R8VEC.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2^31 - 1 )
//      unif = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    19 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input/output, int &SEED, a seed for the random number generator.
//
//    Output, double R8VEC_UNIFORM_01[N], the vector of pseudorandom values.
//
{
	int i;
	int k;
	double *r;

	if (seed == 0)
	{
		cerr << "\n";
		cerr << "R8VEC_UNIFORM_01 - Fatal error!\n";
		cerr << "  Input value of SEED = 0.\n";
		exit(1);
	}

	r = new double[n];

	for (i = 0; i < n; i++)
	{
		k = seed / 127773;

		seed = 16807 * (seed - k * 127773) - k * 2836;

		if (seed < 0)
		{
			seed = seed + 2147483647;
		}

		r[i] = (double)(seed) * 4.656612875E-10;
	}

	return r;
}
//****************************************************************************80

inline void GMRESLib::rearrange_cr(int n, int nz_num, int ia[], int ja[], double a[])

//****************************************************************************80
//
//  Purpose:
//
//    REARRANGE_CR sorts a sparse compressed row matrix.
//
//  Discussion:
//
//    This routine guarantees that the entries in the CR matrix
//    are properly sorted.
//
//    After the sorting, the entries of the matrix are rearranged in such
//    a way that the entries of each column are listed in ascending order
//    of their column values.
//
//    The matrix A is assumed to be stored in compressed row format.  Only
//    the nonzero entries of A are stored.  The vector JA stores the
//    column index of the nonzero value.  The nonzero values are sorted
//    by row, and the compressed row vector IA then has the property that
//    the entries in A and JA that correspond to row I occur in indices
//    IA[I] through IA[I+1]-1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 July 2007
//
//  Author:
//
//    Original C version by Lili Ju.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int N, the order of the system.
//
//    Input, int NZ_NUM, the number of nonzeros.
//
//    Input, int IA[N+1], the compressed row index.
//
//    Input/output, int JA[NZ_NUM], the column indices.  On output,
//    the order of the entries of JA may have changed because of the sorting.
//
//    Input/output, double A[NZ_NUM], the matrix values.  On output, the
//    order of the entries may have changed because of the sorting.
//
{
	double dtemp;
	//int i;
	int is;
	int itemp;
	//int j;
	int j1;
	int j2;
	//int k;

	for (int i = 0; i < n; i++)
	{
		j1 = ia[i];
		j2 = ia[i + 1];
		is = j2 - j1;

		for (int k = 1; k < is; k++)
		{
			for (int j = j1; j < j2 - k; j++)
			{
				if (ja[j + 1] < ja[j])
				{
					itemp = ja[j + 1];
					ja[j + 1] = ja[j];
					ja[j] = itemp;

					dtemp = a[j + 1];
					a[j + 1] = a[j];
					a[j] = dtemp;
				}
			}
		}
	}
	return;
}
//****************************************************************************80



void GMRESLib::test01()

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 tests MGMRES_ST on the simple -1,2-1 matrix.
//
//  Discussion:
//
//    This is a very weak test, since the matrix has such a simple
//    structure, is diagonally dominant (though not strictly), 
//    and is symmetric.
//
//    To make the matrix bigger, simply increase the value of N.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    13 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 2000
# define NZ_NUM 3 * N - 2

	double a[NZ_NUM];
	int i;
	int ia[NZ_NUM];
	int itr_max;
	int j;
	int ja[NZ_NUM];
	int k;
	int mr;
	int n = N;
	int nz_num = NZ_NUM;
	double rhs[N];
	int test;
	double tol_abs;
	double tol_rel;
	double x_error;
	double x_estimate[N];
	double x_exact[N];

	cout << "\n";
	cout << "TEST01\n";
	cout << "  Test MGMRES_ST on the simple -1,2-1 matrix.\n";
	//
	//  Set the matrix.
	//  Note that we use zero based index values in IA and JA.
	//
	k = 0;

	for (i = 0; i < n; i++)
	{
		if (0 < i)
		{
			ia[k] = i;
			ja[k] = i - 1;
			a[k] = -1.0;
			k = k + 1;
		}

		ia[k] = i;
		ja[k] = i;
		a[k] = 2.0;
		k = k + 1;

		if (i < n - 1)
		{
			ia[k] = i;
			ja[k] = i + 1;
			a[k] = -1.0;
			k = k + 1;
		}

	}
	//
	//  Set the right hand side:
	//
	for (i = 0; i < n - 1; i++)
	{
		rhs[i] = 0.0;
	}
	rhs[N - 1] = (double)(n + 1);
	//
	//  Set the exact solution.
	//
	for (i = 0; i < n; i++)
	{
		x_exact[i] = (double)(i + 1);
	}

	for (test = 1; test <= 3; test++)
	{
		//
		//  Set the initial solution estimate.
		//
		for (i = 0; i < n; i++)
		{
			x_estimate[i] = 0.0;
		}

		x_error = 0.0;
		for (i = 0; i < n; i++)
		{
			x_error = x_error + pow(x_exact[i] - x_estimate[i], 2);
		}
		x_error = sqrt(x_error);

		if (test == 1)
		{
			itr_max = 1;
			mr = n / itr_max;
		}
		else if (test == 2)
		{
			itr_max = 2;
			mr = n / itr_max;
		}
		else if (test == 3)
		{
			itr_max = 5;
			mr = n / itr_max;
		}
		tol_abs = 1.0E-08;
		tol_rel = 1.0E-08;

		cout << "\n";
		cout << "  Test " << test << "\n";
		cout << "  Matrix order N = " << n << "\n";
		cout << "  Inner iteration limit = " << mr << "\n";
		cout << "  Outer iteration limit = " << itr_max << "\n";
		cout << "  Initial X_ERROR = " << x_error << "\n";

		mgmres_st(n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr,
			tol_abs, tol_rel);

		x_error = 0.0;
#pragma omp parallel for reduction(+:x_error)
		for (i = 0; i < n; i++)
		{
			x_error += pow(x_exact[i] - x_estimate[i], 2);
		}
		x_error = sqrt(x_error);

		cout << "  Final X_ERROR = " << x_error << "\n";
	}
	return;
# undef N
# undef NZ_NUM
}
//****************************************************************************80

void GMRESLib::test02()

//****************************************************************************80
//
//  Purpose:
//
//    TEST02 tests MGMRES_ST on a 9 by 9 matrix.
//
//  Discussion:
//
//    A = 
//      2  0  0 -1  0  0  0  0  0
//      0  2 -1  0  0  0  0  0  0
//      0 -1  2  0  0  0  0  0  0
//     -1  0  0  2 -1  0  0  0  0
//      0  0  0 -1  2 -1  0  0  0
//      0  0  0  0 -1  2 -1  0  0
//      0  0  0  0  0 -1  2 -1  0
//      0  0  0  0  0  0 -1  2 -1
//      0  0  0  0  0  0  0 -1  2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    18 October 2006
//
//  Author:
//
//    John Burkardt
//
{
# define N 9
# define NZ_NUM 23

	double a[NZ_NUM] = {
		2.0, -1.0,
		2.0, -1.0,
		-1.0, 2.0,
		-1.0, 2.0, -1.0,
		-1.0, 2.0, -1.0,
		-1.0, 2.0, -1.0,
		-1.0, 2.0, -1.0,
		-1.0, 2.0, -1.0,
		-1.0, 2.0 };
	int i;
	int ia[NZ_NUM] = {
		0, 0,
		1, 1,
		2, 2,
		3, 3, 3,
		4, 4, 4,
		5, 5, 5,
		6, 6, 6,
		7, 7, 7,
		8, 8 };
	int itr_max;
	int j;
	int ja[NZ_NUM] = {
		0, 3,
		1, 2,
		1, 2,
		0, 3, 4,
		3, 4, 5,
		4, 5, 6,
		5, 6, 7,
		6, 7, 8,
		7, 8 };
	int k;
	int mr;
	int n = N;
	int nz_num = NZ_NUM;
	double rhs[N] = {
		1.0,
		1.0,
		1.0,
		1.0,
		1.0,
		1.0,
		1.0,
		1.0,
		1.0 };
	int seed = 123456789;
	int test;
	double tol_abs;
	double tol_rel;
	double x_error;
	double *x_estimate;
	double x_exact[N] = {
		3.5,
		1.0,
		1.0,
		6.0,
		7.5,
		8.0,
		7.5,
		6.0,
		3.5 };

	cout << "\n";
	cout << "TEST02\n";
	cout << "  Test MGMRES_ST on matrix that is not quite the -1,2,-1 matrix,\n";
	cout << "  of order N = " << n << "\n";

	for (test = 1; test <= 2; test++)
	{
		if (test == 1)
		{
			cout << "\n";
			cout << "  First try, use zero initial vector:\n";
			cout << "\n";

			x_estimate = new double[n];
#pragma omp parallel for
			for (i = 0; i < n; i++)
			{
				x_estimate[i] = 0.0;
			}
		}
		else
		{
			cout << "\n";
			cout << "  Second try, use random initial vector:\n";
			cout << "\n";

			x_estimate = r8vec_uniform_01(n, seed);
		}
		//
		//  Set the initial solution estimate.
		//
		x_error = 0.0;
#pragma omp parallel for reduction(+:x_error)
		for (i = 0; i < n; i++)
		{
			x_error = x_error + pow(x_exact[i] - x_estimate[i], 2);
		}
		x_error = sqrt(x_error);

		cout << "  Before solving, X_ERROR = " << x_error << "\n";

		itr_max = 20;
		mr = n - 1;
		tol_abs = 1.0E-08;
		tol_rel = 1.0E-08;

		mgmres_st(n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr,
			tol_abs, tol_rel);

		x_error = 0.0;
#pragma omp parallel for reduction(+:x_error)
		for (i = 0; i < N; i++)
		{
			x_error = x_error + pow(x_exact[i] - x_estimate[i], 2);
		}
		x_error = sqrt(x_error);

		cout << "  After solving, X_ERROR = " << x_error << "\n";

		cout << "\n";
		cout << "  Final solution estimate:\n";
		cout << "\n";
		for (i = 0; i < n; i++)
		{
			cout << "  " << setw(8) << i
				<< "  " << setw(12) << x_estimate[i] << "\n";
		}

		delete[] x_estimate;
	}

	return;
# undef N
# undef NZ_NUM
}
//******************************************************************************

void GMRESLib::test03()

//******************************************************************************
//
//  Purpose:
//
//    TEST03 tests PMGMRES_ILU_CR on the simple -1,2-1 matrix.
//
//  Discussion:
//
//    This is a very weak test, since the matrix has such a simple
//    structure, is diagonally dominant (though not strictly), 
//    and is symmetric.
//
//    To make the matrix bigger, simply increase the value of N.
//
//    Note that PGMRES_ILU_CR expects the matrix to be stored using the
//    sparse compressed row format.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 20
# define NZ_NUM 3 * N - 2

	double a[NZ_NUM];
	int i;
	int ia[N + 1];
	int itr_max;
	int j;
	int ja[NZ_NUM];
	int k;
	int mr;
	int n = N;
	int nz_num = NZ_NUM;
	double rhs[N];
	int test;
	double tol_abs;
	double tol_rel;
	double x_error;
	double x_estimate[N];
	double x_exact[N];

	cout << "\n";
	cout << "TEST03\n";
	cout << "  Test PMGMRES_ILU_CR on the simple -1,2-1 matrix.\n";
	//
	//  Set the matrix.
	//  Note that we use zero based index valuesin IA and JA.
	//
	k = 0;
	ia[0] = 0;

	cout << "\n";
	cout << "  ia[" << 0 << "] = " << ia[0] << "\n";
	for (i = 0; i < n; i++)
	{
		ia[i + 1] = ia[i];
		if (0 < i)
		{
			ia[i + 1] = ia[i + 1] + 1;
			ja[k] = i - 1;
			a[k] = -1.0;
			k = k + 1;
		}

		ia[i + 1] = ia[i + 1] + 1;
		ja[k] = i;
		a[k] = 2.0;
		k = k + 1;

		if (i < N - 1)
		{
			ia[i + 1] = ia[i + 1] + 1;
			ja[k] = i + 1;
			a[k] = -1.0;
			k = k + 1;
		}
		cout << "  ia[" << i + 1 << "] = " << ia[i + 1] << "\n";
	}
	//
	//  Set the right hand side:
	//
#pragma omp parallel for
	for (i = 0; i < n - 1; i++)
	{
		rhs[i] = 0.0;
	}
	rhs[n - 1] = (double)(n + 1);
	//
	//  Set the exact solution.
	//
#pragma omp parallel for
	for (i = 0; i < n; i++)
	{
		x_exact[i] = (double)(i + 1);
	}

	for (test = 1; test <= 3; test++)
	{
		//
		//  Set the initial solution estimate.
		//
#pragma omp parallel for
		for (i = 0; i < n; i++)
		{
			x_estimate[i] = 0.0;
		}
		x_error = 0.0;
#pragma omp parallel for reduction(+:x_error)
		for (i = 0; i < n; i++)
		{
			x_error = x_error + pow(x_exact[i] - x_estimate[i], 2);
		}
		x_error = sqrt(x_error);

		if (test == 1)
		{
			itr_max = 1;
			mr = 20;
		}
		else if (test == 2)
		{
			itr_max = 2;
			mr = 10;
		}
		else if (test == 3)
		{
			itr_max = 5;
			mr = 4;
		}
		tol_abs = 1.0E-08;
		tol_rel = 1.0E-08;

		cout << "\n";
		cout << "  Test " << test << "\n";
		cout << "  Matrix order N = " << n << "\n";
		cout << "  Inner iteration limit = " << mr << "\n";
		cout << "  Outer iteration limit = " << itr_max << "\n";
		cout << "  Initial X_ERROR = " << x_error << "\n";

		pmgmres_ilu_cr(n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr,
			tol_abs, tol_rel);

		x_error = 0.0;
#pragma omp parallel for reduction(+:x_error)
		for (i = 0; i < n; i++)
		{
			x_error = x_error + pow(x_exact[i] - x_estimate[i], 2);
		}
		x_error = sqrt(x_error);

		cout << "  Final X_ERROR = " << x_error << "\n";
	}

	return;
# undef N
# undef NZ_NUM
}
//******************************************************************************

void GMRESLib::test04()

//******************************************************************************
//
//  Purpose:
//
//    TEST04 tests PMGMRES_ILU_CR on a simple 5 by 5 matrix.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    25 July 2007
//
//  Author:
//
//    John Burkardt
//
{
# define N 5
# define NZ_NUM 9

	double a[NZ_NUM] = {
		1.0, 2.0, 1.0,
		2.0,
		3.0, 3.0,
		4.0,
		1.0, 5.0 };
	int i;
	int ia[N + 1] = { 0, 3, 4, 6, 7, 9 };
	int itr_max;
	int j;
	int ja[NZ_NUM] = {
		0, 3, 4,
		1,
		0, 2,
		3,
		1, 4 };
	int k;
	int mr;
	int n = N;
	int nz_num = NZ_NUM;
	double rhs[N] = { 14.0, 4.0, 12.0, 16.0, 27.0 };
	int test;
	double tol_abs;
	double tol_rel;
	double x_error;
	double x_estimate[N];
	double x_exact[N] = { 1.0, 2.0, 3.0, 4.0, 5.0 };

	cout << "\n";
	cout << "TEST04\n";
	cout << "  Test PMGMRES_ILU_CR on a simple 5 x 5 matrix.\n";

	cout << "\n";
	for (i = 0; i <= n + 1; i++)
	{
		cout << "  ia[" << i << "] = " << ia[i] << "\n";
	}
	for (test = 1; test <= 3; test++)
	{
		//
		//  Set the initial solution estimate.
		//
#pragma omp parallel for
		for (i = 0; i < n; i++)
		{
			x_estimate[i] = 0.0;
		}
		x_error = 0.0;
#pragma omp parallel for reduction(+:x_error)
		for (i = 0; i < n; i++)
		{
			x_error += pow(x_exact[i] - x_estimate[i], 2);
		}
		x_error = sqrt(x_error);

		if (test == 1)
		{
			itr_max = 1;
			mr = 20;
		}
		else if (test == 2)
		{
			itr_max = 2;
			mr = 10;
		}
		else if (test == 3)
		{
			itr_max = 5;
			mr = 4;
		}
		tol_abs = 1.0E-08;
		tol_rel = 1.0E-08;

		cout << "\n";
		cout << "  Test " << test << "\n";
		cout << "  Matrix order N = " << n << "\n";
		cout << "  Inner iteration limit = " << mr << "\n";
		cout << "  Outer iteration limit = " << itr_max << "\n";
		cout << "  Initial X_ERROR = " << x_error << "\n";

		pmgmres_ilu_cr(n, nz_num, ia, ja, a, x_estimate, rhs, itr_max, mr,
			tol_abs, tol_rel);

		x_error = 0.0;
#pragma omp parallel for reduction(+:x_error)
		for (i = 0; i < n; i++)
		{
			x_error += pow(x_exact[i] - x_estimate[i], 2);
		}
		x_error = sqrt(x_error);

		cout << "  Final X_ERROR = " << x_error << "\n";
	}

	return;
# undef N
# undef NZ_NUM
}