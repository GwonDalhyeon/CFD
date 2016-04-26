//#pragma once
//
//#include "BiCGLib.h"
//
//
//class BiCGSolver
//{
//public:
//	BiCGSolver();
//	~BiCGSolver();
//
//	static void SolverFull(const Array2D<double> & A, const VectorND<double> & b,
//		VectorND<double> & x);
//
//	static void SolverSparse(const Array2D<double> & A, const VectorND<double> & b,
//		VectorND<double> & x);
//
//	static void SolverSparse(const int & rowNum, const VectorND<double> & A,
//		const VectorND<int> & row, const VectorND<int> & col,
//		const VectorND<double> & b, VectorND<double> & x);
//
//	static void SparseA(const Array2D<double> & A, VectorND<double> & a,
//		VectorND<int> & row, VectorND<int> & col, int & nonzeroNum);
//
//private:
//
//};
//
//BiCGSolver::BiCGSolver()
//{
//}
//
//BiCGSolver::~BiCGSolver()
//{
//}
//
//inline void BiCGSolver::SparseA(const Array2D<double>& A, VectorND<double>& a, VectorND<int>& row, VectorND<int>& col, int & nonzeroNum)
//{
//	clock_t before;
//	double  result;
//	before = clock();
//	cout << "Start : Make sparse matrix." << endl;
//
//	VectorND<double> tempA(0, A.iRes*A.jRes);
//	VectorND<int> tempRow(0, A.iRes*A.jRes);
//	VectorND<int> tempCol(0, A.iRes*A.jRes);
//
//	nonzeroNum = 0;
//	for (int i = A.iStart; i <= A.iEnd; i++)
//	{
//		for (int j = A.jStart; j <= A.jEnd; j++)
//		{
//			if (abs(A(i, j)) > 0)
//			{
//				tempA(nonzeroNum) = A(i, j);
//				tempRow(nonzeroNum) = i - A.iStart;
//				tempCol(nonzeroNum) = j - A.jStart;
//				nonzeroNum++;
//			}
//		}
//	}
//
//	a = VectorND<double>(0, nonzeroNum);
//	row = VectorND<int>(0, nonzeroNum);
//	col = VectorND<int>(0, nonzeroNum);
//#pragma omp parallel for
//	for (int i = 0; i < nonzeroNum; i++)
//	{
//		a(i) = tempA(i);
//		row(i) = tempRow(i);
//		col(i) = tempCol(i);
//	}
//
//	result = (double)(clock() - before) / CLOCKS_PER_SEC;
//	cout << "time : " << result << "\n";
//	cout << "End : Make sparse matrix." << endl;
//	cout << endl;
//}
