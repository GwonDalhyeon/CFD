#pragma once
#include "CommonDef.h"
#include "VectorND.h"
#include "CSR.h"

class PCGSolver
{
public:
	PCGSolver();
	~PCGSolver();

	static void Solver(const CSR<double>& A, const CSR<double>& M, const VectorND<double>& b, VectorND<double>& x, const Array2D<int> & bc_input);

	static void IncompleteCholeskyDecomposition(const CSR<double>& A, CSR<double>& L, const Array2D<int>& bc_input);
	static void MultiplicationByMinverse(const CSR<double>& M, const VectorND<double>& b, VectorND<double>& x);
	static void MultiplicationByMinverse(const CSR<double>& M, const VectorND<double>& Mdiag, const VectorND<double>& b, VectorND<double>& x);
	static void MultiplicationByMDiagonalInverse(const CSR<double>& M, const VectorND<double>& b, VectorND<double>& x);
private:

};

PCGSolver::PCGSolver()
{
}

PCGSolver::~PCGSolver()
{
}

inline void PCGSolver::Solver(const CSR<double>& A, const CSR<double>& M, const VectorND<double>& b, VectorND<double>& x, const Array2D<int> & bc_input)
{
	clock_t before;
	double  result;
	before = clock();
	cout << "--------  Start : PCG  --------" << endl;

	x = 0;
	const int N = x.iLength;
	int	num_iteration = 0;
	double alpha, beta, res_old, res_new = 0, dot_result;

	VTN res(N);
	VTN Ap(N);
	A.ComputeResidual(x, b, res);

	VTN z(N);
	VTN Mdiag(N);

#pragma omp parallel for
	for (int i = 0; i < M.rowNum; i++)
	{
		Mdiag[i] = M(i, i);
	}
	cout << "Time: " << (double)(clock() - before) / CLOCKS_PER_SEC << "\n";

	clock_t bb = clock();
	MultiplicationByMinverse(M, Mdiag, res, z);
	cout << "Time: " << (double)(clock() - bb) / CLOCKS_PER_SEC << "\n";

	VTN p(z);
	
	res_old = DotProduct(res, z);
	while (num_iteration < 2 * A.rowNum)
	{
		A.Multiply(p, Ap);

		dot_result = DotProduct(p, Ap);
		if (num_iteration == 0 && abs(dot_result) < DBL_EPSILON)
		{
			cout << "First Time!!" << endl;
			cout << endl;
			res_new = 0;
			break;
		}
		alpha = res_old / dot_result;

#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{
			x[i] += alpha*p[i];
			res[i] -= alpha*Ap[i];
		}

		if (res_old < DBL_EPSILON)
		{
			cout << "Converge!!" << endl;
			break;
		}

		MultiplicationByMinverse(M, Mdiag, res, z);

		res_new = DotProduct(res, z);
		beta = res_new / res_old;
#pragma omp parallel for
		for (int i = 0; i < N; i++)
		{
			p[i] = z[i] + beta*p[i];
		}

		res_old = res_new;

		num_iteration++;
	}


	cout << "Iteration Number : " << num_iteration << endl;
	cout << "Time             : " << (double)(clock() - before) / CLOCKS_PER_SEC << "\n";
	cout << "Residual         : " << sqrt(res_new) << endl;
	cout << "--------  End : CG  -------- " << endl;
	cout << endl;
}

inline void PCGSolver::IncompleteCholeskyDecomposition(const CSR<double>& A, CSR<double>& L, const Array2D<int>& bc_input)
{
	const int N = A.rowNum;
	const int nz = A.valueNum;
	int iStart = bc_input.iStart, iEnd = bc_input.iEnd, jStart = bc_input.jStart, jEnd = bc_input.jEnd;
	int iRes = bc_input.iRes;
	int number = 0;
	bool tempBool = true;
	double sum, coef;
	for (int j = jStart; j <= jEnd; j++)
	{
		for (int i = iStart; i <= iEnd; i++)
		{
			sum = 0, coef = 0;

		
			if (j > jStart)
			{
				if (bc_input(i, j - 1) > -1)
				{
					if (bc_input(i, j) == iRes)
					{
						coef = 1 / L.values[L.indPrt[bc_input(i, j) - iRes]] * (A(bc_input(i, j - 1), bc_input(i, j)));
					}
					else if (((bc_input(i, j) > iRes) && (bc_input(i, j) < 2 * iRes)) || (bc_input(i, j) % iRes == 0))
					{
						coef = 1 / L.values[L.indPrt[bc_input(i, j) - iRes] + 1] * (A(bc_input(i, j - 1), bc_input(i, j)));
					}
					else
					{
						tempBool = true;
						if (i < iEnd)
						{
							if (bc_input(i + 1, j) == BC_PER)
							{
								coef = 1 / L.values[L.indPrt[bc_input(i, j) - iRes] + 3] * (A(bc_input(i, j - 1), bc_input(i, j)));
								tempBool = false;
							}
						}
						if (tempBool)
						{
							coef = 1 / L.values[L.indPrt[bc_input(i, j) - iRes] + 2] * (A(bc_input(i, j - 1), bc_input(i, j)));
						}
					}

					L.AssignValue(bc_input(i, j), bc_input(i, j - 1), coef);
					number += 1;
				}
			}
			
	
			if (i<iEnd)
			{
				if (bc_input(i + 1, j) == BC_PER)
				{
					coef = 1 / L.values[L.indPrt[bc_input(i, j) - 2] - 1] * (A(bc_input(i, j) - (iRes - 1), bc_input(i, j)));

					L.AssignValue(bc_input(i, j), bc_input(i, j) - (iRes - 1), coef);
					number += 1;
				}
			}
			

			if (i > iStart)
			{
				if (bc_input(i - 1, j) > -1)
				{
					tempBool = true;
					if (i<iEnd)
					{
						if (bc_input(i + 1, j) == BC_PER)
						{
							coef = 1 / L.values[L.indPrt[bc_input(i, j) - 1]] * (A(bc_input(i - 1, j), bc_input(i, j)));
							tempBool = false;
						}
					}
					if (tempBool)
					{
						if (bc_input(i, j) == 1)
						{
							coef = 1 / L.values[L.indPrt[bc_input(i, j) - 1]] * (A(bc_input(i - 1, j), bc_input(i, j)));
						}
						else if (bc_input(i, j) < iRes)
						{
							coef = 1 / L.values[L.indPrt[bc_input(i, j) - 1] + 1] * (A(bc_input(i - 1, j), bc_input(i, j)));
						}
						else
						{
							coef = 1 / L.values[L.indPrt[bc_input(i, j)] - 1] * (A(bc_input(i - 1, j), bc_input(i, j)));
						}
					}
					

					L.AssignValue(bc_input(i, j), bc_input(i - 1, j), coef);
					number += 1;
				}
			}
			

			if (bc_input(i, j) == 0)
			{
				coef = sqrt(A(bc_input(i, j), bc_input(i, j)));
			}
			else
			{
				sum = 0;

				for (int k = L.indPrt[bc_input(i, j)]; k < number; k++)
				{
					sum += L.values[k] * L.values[k];
				}
				coef = sqrt(A(bc_input(i, j), bc_input(i, j)) - sum);
			}
			L.AssignValue(bc_input(i, j), bc_input(i, j), coef);

			number += 1;
		}
	}
}

inline void PCGSolver::MultiplicationByMinverse(const CSR<double>& M, const VectorND<double>& b, VectorND<double>& x)
{
	double* y = new double[x.iLength];
	
	double one_over_M_start = 1 / M(0, 0);

	y[0] = b[0] * one_over_M_start;

	double sum;
	double one_over_Mii;
	int number(0), num_2(0);
	for (int i = 1; i < x.iLength; i++)
	{
		sum = 0;
		for (int k = M.indPrt[i]; k < (M.indPrt[i + 1] - 1); k++)
		{
			sum += M.values[k] * y[M.columns[k]];
		}

		one_over_Mii = 1 / M(i, i);
		y[i] = (b[i] - sum)*one_over_Mii;
	}

	double* summation = new double[x.iLength];
	for (int i = 0; i < x.iLength; i++)
	{
		summation[i] = 0;
	}

	// Matrix-Transpose-Vector Multiplication
	// Parallel Sparse Matrix-Vector and Matrix-Trnaspose-Vector Multiplication Using Compressed Sparse Blocks
	double one_over_M;
	double one_over_M_end = 1 / M(x.iLength - 1, x.iLength - 1);
	x[x.iLength - 1] = y[x.iLength - 1] * one_over_M_end;
	for (int i = x.iLength - 1; i >= 1; i--)
	{
		for (int k = M.indPrt[i + 1] - 1; k >= M.indPrt[i]; k--)
		{
			if (k != M.indPrt[i + 1] - 1)
			{
				summation[M.columns[k]] += M.values[k] * x[i];
			}
		}

		one_over_M = 1 / M(i - 1, i - 1);
		x[i - 1] = (y[i - 1] - summation[i - 1])*one_over_M;
	}

	delete[] y;
	delete[] summation;
}

inline void PCGSolver::MultiplicationByMinverse(const CSR<double>& M, const VectorND<double>& Mdiag, const VectorND<double>& b, VectorND<double>& x)
{
	//VTN y(x.iLength);
	double* y = new double[x.iLength];
	double one_over_M_start = 1 / Mdiag[0];

	y[0] = b[0] * one_over_M_start;

	double sum;
	double one_over_Mii;
	int number(0), num_2(0);
	for (int i = 1; i < x.iLength; i++)
	{
		sum = 0;
		for (int k = M.indPrt[i]; k < (M.indPrt[i + 1] - 1); k++)
		{
			sum += M.values[k] * y[M.columns[k]];
		}

		one_over_Mii = 1 / Mdiag(i);
		y[i] = (b[i] - sum)*one_over_Mii;
	}

	double* summation = new double[x.iLength];
#pragma omp parallel for
	for (int i = 0; i < x.iLength; i++)
	{
		summation[i] = 0;
	}
	//// Matrix-Transpose-Vector Multiplication
	//// Parallel Sparse Matrix-Vector and Matrix-Trnaspose-Vector Multiplication Using Compressed Sparse Blocks
	double one_over_M;
	double one_over_M_end = 1 / Mdiag(x.iLength - 1);
	x[x.iLength - 1] = y[x.iLength - 1] * one_over_M_end;
	for (int i = x.iLength - 1; i >= 1; i--)
	{
		for (int k = M.indPrt[i + 1] - 1; k >= M.indPrt[i]; k--)
		{
			if (k != M.indPrt[i + 1] - 1)
			{
				summation[M.columns[k]] += M.values[k] * x[i];
			}
		}

		one_over_M = 1 / Mdiag(i - 1);
		x[i - 1] = (y[i - 1] - summation[i - 1])*one_over_M;
	}

	delete[] y;
	delete[] summation;
}

inline void PCGSolver::MultiplicationByMDiagonalInverse(const CSR<double>& M, const VectorND<double>& b, VectorND<double>& x)
{
	double oneOverMii;
#pragma omp parallel for private(oneOverMii)
	for (int i = 0; i <x.iLength; i++)
	{
		oneOverMii = 1 / M(i, i);
		x[i] = oneOverMii*b[i];
	}
}
