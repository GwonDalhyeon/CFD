#pragma once
#include "CG.h"

class CGSolver
{
public:
	CGSolver();
	~CGSolver();

	static void SolverFull(const Array2D<double> & A, const VectorND<double> & b,
		VectorND<double> & x);

	static void SolverSparse(const Array2D<double> & A, const VectorND<double> & b,
		VectorND<double> & x);

	static void SolverSparse(const int & rowNum, const VectorND<double> & A,
		const VectorND<int> & row, const VectorND<int> & col,
		const VectorND<double> & b,	VectorND<double> & x);

	static VectorND<double> SolverCSR(const CSR<double>& A, const VectorND<double>& b);
	static VectorND<double> SolverCSR(const CSR<double>& A, const VectorND<double>& b, const double & tol);

	static void SparseA(const Array2D<double> & A, VectorND<double> & a,
		VectorND<int> & row, VectorND<int> & col, int & nonzeroNum);
private:

};

CGSolver::CGSolver()
{
}

CGSolver::~CGSolver()
{
}

inline void CGSolver::SolverFull(const Array2D<double>& A, const VectorND<double>& b,
	VectorND<double>& x)
{
	assert(A.iRes == A.jRes);
	assert(A.iRes == b.iLength);
	assert(A.iRes == x.iLength);
	clock_t before;
	double  result;
	before = clock();
	cout << "Start : CG" << endl;

	CGLib::r8ge_cg(A.iRes, A.values, b.values, x.values);

	result = (double)(clock() - before) / CLOCKS_PER_SEC;
	cout << "time : " << result << "\n";
	//printf("�ɸ��ð��� %5.2f �Դϴ�.\n", result);

	double *r;
	r = CGLib::r8ge_res(A.iRes, A.iRes, A.values, x.values, b.values);
	double r_norm = CGLib::r8vec_norm(A.iRes, r);
	cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";

	delete[] r;
}

inline void CGSolver::SolverSparse(const Array2D<double>& A, const VectorND<double>& b,
	VectorND<double>& x)
{
	assert(A.iRes == A.jRes);
	assert(A.iRes == b.iLength);
	assert(A.iRes == x.iLength);

	VectorND<double> a;
	VectorND<int> row;
	VectorND<int> col;
	int nonzeroNum;


	SparseA(A, a, row, col, nonzeroNum);

	clock_t before;
	double  result;
	before = clock();
	cout << "Start : CG" << endl;

	CGLib::r8sp_cg(A.iRes, nonzeroNum, row.values, col.values, a.values, b.values, x.values);

	result = (double)(clock() - before) / CLOCKS_PER_SEC;
	cout << "time : " << result << "\n";
	//printf("�ɸ��ð��� %5.2f �Դϴ�.\n", result);

	double *r;
	r = CGLib::r8ge_res(A.iRes, A.iRes, A.values, x.values, b.values);
	double r_norm = CGLib::r8vec_norm(A.iRes, r);
	cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";


	delete[] r;
}

inline void CGSolver::SolverSparse(const int & rowNum, const VectorND<double>& A,
	const VectorND<int> & row, const VectorND<int> & col,
	const VectorND<double>& b, VectorND<double>& x)
{
	clock_t before;
	double  result;
	before = clock();
	cout << "Start : CG" << endl;

	CGLib::r8sp_cg(rowNum, A.iLength, row.values, col.values, A.values, b.values, x.values);

	result = (double)(clock() - before) / CLOCKS_PER_SEC;
	cout << "time : " << result << "\n";

	double * resB = CGLib::r8sp_mv(b.iLength, b.iLength, A.iLength, row.values, col.values, A.values, x.values);

	double *r = new double[b.iLength];
#pragma omp parallel for
	for (int i = 0; i < b.iLength; i++)
	{
		r[i] = 0;
		r[i] = b.values[i] - resB[i];
	}
	double r_norm = CGLib::r8vec_norm(b.iLength, r);
	cout << "  Norm of residual ||Ax-b|| = " << r_norm << "\n";


	delete[] r;
}

inline VectorND<double> CGSolver::SolverCSR(const CSR<double>& A, const VectorND<double>& b)
{
	clock_t before;
	double  result;
	before = clock();
	cout << "Start : CG " << endl;

	int num = A.rowNum;
	double tolerance = 10e-5; 1000 * DBL_EPSILON;

	VectorND<double> rOld(num);
	VectorND<double> p(num);
	VectorND<double> rNew(num);
	VectorND<double> x(num);

	int j;

	rOld = b;
	p = b;
	x = 0;

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double temp = 0;
	double residual;
	double residualOld = rOld.magnitude2();

	for (int k = 0; k < 2 * num; k++)
	{
		temp1 = 0;
		temp2 = 0;
		//#pragma omp parallel for private (j) reduction(+:temp2)
		for (int i = 0; i < num; i++)
		{
			//A.indPrt
			for (int n = A.indPrt[i]; n < A.indPrt[i + 1]; n++)
			{
				j = int(A.columns[n]);
				temp2 = temp2 + p[i] * A.values[n] * p[j];
			}
		}
		alpha = residualOld / temp2;

		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] + alpha*p[i];
			temp = 0;
#pragma omp parallel for private (j) reduction(+:temp) 
			for (int n = A.indPrt[i]; n < A.indPrt[i + 1]; n++)
			{
				j = int(A.columns[n]);
				temp = temp + A.values[n] * p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = rNew.magnitude2();

		temp = sqrt(abs(residual));

		//cout << k << " " << temp << endl;
		if (temp < tolerance)
		{
			cout << "CG iterataion : " << k << " " << temp << endl;
			result = (double)(clock() - before) / CLOCKS_PER_SEC;
			cout << "time : " << result << "\n";
			cout << "End : CG " << endl;
			cout << endl;
			return x;
		}

		beta = residual / residualOld;
#pragma omp parallel for
		for (int i = 0; i < p.iLength; i++)
		{
			p[i] = rNew[i] + beta*p[i];
		}
		rOld = rNew;

		residualOld = residual;

	}

	result = (double)(clock() - before) / CLOCKS_PER_SEC;
	cout << "time : " << result << "\n";
	cout << "End : CG " << endl;
	cout << endl;

	return x;
}

inline VectorND<double> CGSolver::SolverCSR(const CSR<double>& A, const VectorND<double>& b, const double & tol)
{
	clock_t before;
	double  result;
	before = clock();
	cout << "Start : CG " << endl;

	int num = A.rowNum;
	double tolerance = tol;

	VectorND<double> rOld(num);
	VectorND<double> p(num);
	VectorND<double> rNew(num);
	VectorND<double> x(num);

	int j;

	rOld = b;
	p = b;
	x = 0;

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double temp = 0;
	double residual;
	double residualOld = rOld.magnitude2();

	for (int k = 0; k < 2 * num; k++)
	{
		temp1 = 0;
		temp2 = 0;
#pragma omp parallel for private (j) reduction(+:temp2)
		for (int i = 0; i < num; i++)
		{
			//A.indPrt
			for (int n = A.indPrt[i]; n < A.indPrt[i + 1]; n++)
			{
				j = int(A.columns[n]);
				temp2 = temp2 + p[i] * A.values[n] * p[j];
			}
		}
		alpha = residualOld / temp2;

		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] + alpha*p[i];
			temp = 0;
#pragma omp parallel for private (j) reduction(+:temp) 
			for (int n = A.indPrt[i]; n < A.indPrt[i + 1]; n++)
			{
				j = int(A.columns[n]);
				temp = temp + A.values[n] * p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = rNew.magnitude2();

		temp = sqrt(abs(residual));

		//cout << k << " " << temp << endl;
		if (temp < tolerance)
		{
			cout << "CG iterataion : " << k << " " << temp << endl;
			result = (double)(clock() - before) / CLOCKS_PER_SEC;
			cout << "time : " << result << "\n";
			cout << "End : CG " << endl;
			cout << endl;
			return x;
		}

		beta = residual / residualOld;
#pragma omp parallel for
		for (int i = 0; i < p.iLength; i++)
		{
			p[i] = rNew[i] + beta*p[i];
		}
		rOld = rNew;

		residualOld = residual;

	}



	return x;
}

inline void CGSolver::SparseA(const Array2D<double> & A, VectorND<double> & a,
	VectorND<int> & row, VectorND<int> & col, int & nonzeroNum)
{
	clock_t before;
	double  result;
	before = clock();
	cout << "Start : Make sparse matrix." << endl;

	VectorND<double> tempA(0, A.iRes*A.jRes);
	VectorND<int> tempRow(0, A.iRes*A.jRes);
	VectorND<int> tempCol(0, A.iRes*A.jRes);

	nonzeroNum = 0;
	for (int i = A.iStart; i <= A.iEnd; i++)
	{
		for (int j = A.jStart; j <= A.jEnd; j++)
		{
			if (abs(A(i, j)) > 0)
			{
				tempA(nonzeroNum) = A(i, j);
				tempRow(nonzeroNum) = i - A.iStart;
				tempCol(nonzeroNum) = j - A.jStart;
				nonzeroNum++;
			}
		}
	}

	a = VectorND<double>(0, nonzeroNum);
	row = VectorND<int>(0, nonzeroNum);
	col = VectorND<int>(0, nonzeroNum);
#pragma omp parallel for
	for (int i = 0; i < nonzeroNum; i++)
	{
		a(i) = tempA(i);
		row(i) = tempRow(i);
		col(i) = tempCol(i);
	}

	result = (double)(clock() - before) / CLOCKS_PER_SEC;
	cout << "time : " << result << "\n";
	cout<<"End : Make sparse matrix." << endl;
	cout << endl;
}