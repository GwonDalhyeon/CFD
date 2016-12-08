#pragma once
#include "CommonDef.h"
#include "VectorND.h"
#include "CSR.h"

class PCGSolver
{
public:
	PCGSolver();
	~PCGSolver();

	static void Solver(const CSR<double>& A, const VectorND<double>& b, VectorND<double>& x, const Field2D<int> & bc_input);

	static void IncompleteCholeskyDecomposition(const int& i_res_input, const int& j_res_input, const CSR<double>& A, CSR<double>& L, const Field2D<int>& bc_input);
	static void MultiplicationByMinverse(const int& i_res_input, const int& j_res_input, const CSR<double>& M, VectorND<double>& x, const VectorND<double>& b);
private:

};

PCGSolver::PCGSolver()
{
}

PCGSolver::~PCGSolver()
{
}

inline void PCGSolver::Solver(const CSR<double>& A, const VectorND<double>& b, VectorND<double>& x, const Field2D<int> & bc_input)
{
	clock_t before;
	double  result;
	before = clock();
	cout << "--------  Start : PCG  --------" << endl;

	const int N = x.iLength;

	VTN res(N);
	VTN s(N);
	VTN p(N);
	VTN Ap(N);

	int	num_iteration = 0;

	//T *rval(res.values), *sval(s.values), *pval(p.values), *Apval(Ap.values), *xval(x.values);

	double alpha, beta, res_old, res_new, dot_result;

	A.ComputeResidual(x, b, res);

	const int nz = A.valueNum;

	const int num_of_nz = (nz - N) / 2 + N;

	CSR<double> M(N, num_of_nz);
	M.Initialize(N, num_of_nz);

	IncompleteCholeskyDecomposition(A.rowNum, A.rowNum, A, M, bc_input);

	/*ofstream fout;
	fout.open("value");
	for (int i = 0; i < M.nz; i++)
	{
	fout << M.values[i] << endl;
	}
	fout.close();

	fout.open("col");
	for (int i = 0; i < M.nz; i++)
	{
	fout << M.column_index[i] << endl;
	}
	fout.close();

	fout.open("row_ptr");
	for (int i = 0; i <= M.N; i++)
	{
	fout << M.row_ptr[i] << endl;
	}
	fout.close();*/

	//MultiplicationByMinverseAsDiagonal(A, p, res);
	MultiplicationByMinverse(A.rowNum, A.rowNum, M, p, res);

	res_new = DotProduct(res, p);
	int max_iteration = A.rowNum*A.rowNum;
	while (true)
	{
		A.Multiply(p, Ap);

		dot_result = DotProduct(p, Ap);
		if (abs(dot_result) < DBL_EPSILON)
		{
			cout << "First Time!!" << endl;
			cout << endl;
			res_new = 0;
			break;
		}

		alpha = res_new / dot_result;

		for (int i = 0; i < N; i++)
		{
			x[i] += alpha*p[i];
			res[i] -= alpha*Ap[i];
		}

		/*if (num_iteration % (int)sqrt(A.N) == 0)
		{
		A.ComputeResidual(x, b, res);
		}
		else
		{
		for (int i = 0; i < N; i++)
		{
		rval[i] -= alpha*Apval[i];
		}
		}*/

		//MultiplicationByMinverseAsDiagonal(A, s, res);
		MultiplicationByMinverse(A.rowNum, A.rowNum, M, s, res);

		res_old = res_new;

		res_new = DotProduct(res, s);

		beta = res_new / res_old;

		for (int i = 0; i < N; i++)
		{
			p[i] *= beta;
			p[i] += s[i];
		}

		num_iteration++;

		if (num_iteration > max_iteration) break;
		if (res_new < DBL_EPSILON)
		{
			cout << "Converge!!" << endl;
			break;
		}
	}


	cout << "Iteration Number : " << num_iteration << endl;
	cout << "Time             : " << (double)(clock() - before) / CLOCKS_PER_SEC << "\n";
	cout << "Residual         : " << sqrt(res_new) << endl;
	cout << "--------  End : CG  -------- " << endl;
	cout << endl;
}

inline void PCGSolver::IncompleteCholeskyDecomposition(const int & i_res_input, const int & j_res_input, const CSR<double>& A, CSR<double>& L, const Field2D<int>& bc_input)
{
	const int N = A.rowNum;
	const int nz = A.valueNum;

	//Field2D<int>& index_field = bc_input;

	//const int i_start(index_field.i_start), i_end(index_field.i_end), j_start(index_field.j_start), j_end(index_field.j_end);

	//int start_ix(0);

	//int i, j;
	//LOOPS_2D(i, j, i_start, j_start, i_end, j_end)
	//{
	//	index_field(i, j) = -1;
	//}

	//LOOPS_2D(i, j, index_field.i_start_g, index_field.j_start_g, index_field.i_end_g, index_field.j_end_g)
	//{
	//	if (i < i_start || i > i_end)
	//	{
	//		index_field(i, j) = bc_input(i, j);
	//	}
	//	else
	//	{
	//		index_field(i, j) = bc_input(i, j);
	//	}
	//}

	//ofstream fout;
	//fout.open("index_field");
	//for (int j = 0; j < index_field.grid.j_res; j++)
	//{
	//for (int i = 0; i < index_field.grid.i_res; i++)
	//{
	//fout << index_field(i, j) << " ";
	//}
	//fout << endl;
	//}
	//fout.close();

	int number = 0;
	double sum, coef;
	for (int i = bc_input.iStart; i <= bc_input.iEnd; i++)
	{
		for (int j = bc_input.jStart; j <= bc_input.jEnd; j++)
		{
			sum = 0, coef = 0;

			if (bc_input(i, j - 1) > -1)
			{
				if (bc_input(i, j) == bc_input.iRes)
				{
					coef = 1 / L.values[L.indPrt[bc_input(i, j) - bc_input.iRes]] * (A(bc_input(i, j - 1), bc_input(i, j)));
				}
				else if (((bc_input(i, j) > bc_input.iRes) && (bc_input(i, j) < 2 * bc_input.iRes)) || (bc_input(i, j) % bc_input.iRes == 0))
				{
					coef = 1 / L.values[L.indPrt[bc_input(i, j) - bc_input.iRes] + 1] * (A(bc_input(i, j - 1), bc_input(i, j)));
				}
				else
				{
					if (bc_input(i + 1, j) == BC_PER)
					{
						coef = 1 / L.values[L.indPrt[bc_input(i, j) - bc_input.iRes] + 3] * (A(bc_input(i, j - 1), bc_input(i, j)));
					}
					else
					{
						coef = 1 / L.values[L.indPrt[bc_input(i, j) - bc_input.iRes] + 2] * (A(bc_input(i, j - 1), bc_input(i, j)));
					}
				}

				L.AssignValue(bc_input(i, j), bc_input(i, j - 1), coef);
				number += 1;




				if (bc_input(i + 1, j) == BC_PER)
				{
					coef = 1 / L.values[L.indPrt[bc_input(i, j) - 2] - 1] * (A(bc_input(i, j) - (bc_input.iRes - 1), bc_input(i, j)));

					L.AssignValue(bc_input(i, j), bc_input(i, j) - (bc_input.iRes - 1), coef);
					number += 1;
				}

				if (bc_input(i - 1, j) > -1)
				{
					if (bc_input(i + 1, j) == BC_PER)
					{
						coef = 1 / L.values[L.indPrt[bc_input(i, j) - 1]] * (A(bc_input(i - 1, j), bc_input(i, j)));
					}
					else
					{
						if (bc_input(i, j) == 1)
						{
							coef = 1 / L.values[L.indPrt[bc_input(i, j) - 1]] * (A(bc_input(i - 1, j), bc_input(i, j)));
						}
						else if (bc_input(i, j) < bc_input.iRes)
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

				if (bc_input(i, j) == 0)
				{
					coef = sqrt(A(bc_input(i, j), bc_input(i, j)));
				}
				else
				{
					sum = 0;

					for (int k = L.indPrt[bc_input(i, j)]; k < number; k++)
					{
						sum += L.values[k]* L.values[k];
					}
					coef = sqrt(A(bc_input(i, j), bc_input(i, j)) - sum);
				}
				L.AssignValue(bc_input(i, j), bc_input(i, j), coef);

				number += 1;
			}
		}

	}
}

inline void PCGSolver::MultiplicationByMinverse(const int & i_res_input, const int & j_res_input, const CSR<double>& M, VectorND<double>& x, const VectorND<double>& b)
{
	VTN y(x.iLength);

	double one_over_M_start = 1 / M(0, 0);

	y[0] = b[0] * one_over_M_start;

	double sum;
	double one_over_Mii;
	int number(0), num_2(0);
#pragma omp parallel for private(sum, one_over_Mii)
	for (int i = 1; i < y.iLength; i++)
	{
		sum = 0;
		for (int k = M.indPrt[i]; k < (M.indPrt[i + 1] - 1); k++)
		{
			sum += M.values[k] * y[M.columns[k]];
		}

		one_over_Mii = 1 / M(i, i);
		y[i] = (b[i] - sum)*one_over_Mii;
	}

	VTN summation(x.iLength);


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
}
