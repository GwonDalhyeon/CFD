#pragma once


#pragma once
#include "CommonDef.h"
#include "VectorND.h"
#include "CSROld.h"
#include "CSR.h"
#include "Field2D.h"
#include "Array2D.h"

template <class TT>
VectorND<TT> CG(const CSR<TT>& A, VectorND<TT> b);

template <class TT>
double* CG(const CSR<TT>& A, double* b);

double* CG(int num, double* A, double* b);

//void incompleteCholesky(int num, double* A);

double* PCG(int num, double* A, double* b);

template <class TT>
static  Array2D<TT> GaussSeidel(const Array2D<TT>& A, const Array2D<TT>& x, Array2D<TT>& b);


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template <class TT>
VectorND<TT> CG(const CSR<TT>& A, VectorND<TT> b)
{
	int num = A.rowNum;
	double tolerance = 1000 * DBL_EPSILON;

	VectorND<TT> rOld(num);
	VectorND<TT> p(num);
	VectorND<TT> rNew(num);
	VectorND<TT> x(num);

	int j;

	rOld = b;
	p = b;
	x = 0;

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double temp = 0;
	double residual;
	TT residualOld = rOld.magnitude2();

	for (int k = 0; k < 2 * num; k++)
	{
		temp1 = 0;
		temp2 = 0;
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
			for (int n = A.indPrt[i]; n < A.indPrt[i + 1]; n++)
			{
				j = int(A.columns[n]);
				temp = temp + A.values[n] * p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = rNew.magnitude2();

		temp = sqrt(abs(residual));

		cout << k << " " << temp << endl;
		if (temp < tolerance)
		{
			cout << "CG iterataion : " << k << endl;
			cout << endl;
			return x;
		}

		beta = residual / residualOld;

		for (int i = 0; i < p.iLength; i++)
		{
			p[i] = rNew[i] + beta*p[i];
		}
		rOld = rNew;

		residualOld = residual;

	}

	return x;

}

template <class TT>
double* CG(const CSR<TT>& A, double* b)
{
	int num = A.colNum;
	double tolerance = 1000 * DBL_EPSILON;
	double* rOld = new double[num];
	double* p = new double[num];
	double* rNew = new double[num];
	double* x = new double[num];
	int j;

	for (int i = 0; i < num; i++)
	{
		rOld[i] = b[i];
		p[i] = b[i];
		x[i] = 0;
	}

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double temp = 0;
	double residual;
	double residualOld = 0;
	for (int i = 0; i < num; i++)
	{
		residualOld = residualOld + rOld[i] * rOld[i];
	}



	for (int k = 0; k < 2 * num; k++)
	{
		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			for (int n = A.indptr[i]; n < A.indptr[i + 1]; n++)
			{
				j = int(A.col[n]);
				temp2 = temp2 + p[i] * A.val[n] * p[j];
			}
		}
		alpha = residualOld / temp2;

		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] + alpha*p[i];
			temp = 0;
			for (int n = A.indptr[i]; n < A.indptr[i + 1]; n++)
			{
				j = int(A.col[n]);
				temp = temp + A.val[n] * p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = 0;
		for (int i = 0; i < num; i++)
		{
			residual = residual + rNew[i] * rNew[i];
			//cout<<rNew[i]<<endl;
		}

		temp = sqrt(abs(residual));
		if (k % 10 == 0)
		{
			cout << ".";
		}
		cout << k << " " << temp << endl;
		if (temp < tolerance)
		{
			delete rNew, rOld, p;
			cout << "CG iterataion : " << k << endl;
			cout << endl;
			return x;
		}

		beta = residual / residualOld;

		for (int i = 0; i < num; i++)
		{
			p[i] = rNew[i] + beta*p[i];
			rOld[i] = rNew[i];
		}
		residualOld = residual;

	}
	delete[] rNew, rOld, p;
	return x;

}

template<class TT>
inline Array2D<TT> GaussSeidel(const Array2D<TT>& A, const Array2D<TT>& x, Array2D<TT>& b)
{
	VectorND<TT> reshapeX = VectorND<TT>(x.ijRes);
	VectorND<TT> reshapeB = VectorND<TT>(b.ijRes);
	CSR<double> Acsr = CSR<double>(A);

	return Array2D<TT>();
}

double* CG(int num, double* A, double* b)
{
	double tolerance = 1000 * DBL_EPSILON;
	double* rOld = new double[num];
	double* p = new double[num];
	double* rNew = new double[num];
	double* x = new double[num];

	for (int i = 0; i < num; i++)
	{
		rOld[i] = b[i];
		p[i] = b[i];
		x[i] = 0;
	}

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double temp = 0;
	double residual;
	double residualOld = 0;
	for (int i = 0; i < num; i++)
	{
		residualOld = residualOld + rOld[i] * rOld[i];
	}

	for (int k = 0; k < 2 * num; k++)
	{
		temp1 = 0;
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			for (int j = 0; j < num; j++)
			{
				temp2 = temp2 + p[i] * A[i*num + j] * p[j];
			}
		}
		alpha = residualOld / temp2;

		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] + alpha*p[i];
			temp = 0;
			for (int j = 0; j < num; j++)
			{
				temp = temp + A[i*num + j] * p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = 0;
		for (int i = 0; i < num; i++)
		{
			residual = residual + rNew[i] * rNew[i];
			//cout<<rNew[i]<<endl;
		}

		temp = sqrt(abs(residual));

		//cout << k << " " << temp << endl;
		if (temp < tolerance)
		{
			delete rNew, rOld, p;
			cout << "CG iterataion : " << k << endl;
			cout << endl;
			return x;
		}

		beta = residual / residualOld;

		for (int i = 0; i < num; i++)
		{
			p[i] = rNew[i] + beta*p[i];
			rOld[i] = rNew[i];
		}
		residualOld = residual;

	}
	delete[] rNew, rOld, p;
	return x;

}

double* incompleteCholesky(int num, double* A)
{
	double* L = new double[num*num];
	for (int i = 0; i < num; i++)
	{
		for (int j = 0; j < num; j++)
		{
			L[i*num + j] = 0;
		}
	}

	double temp;
	for (int i = 0; i < num; i++)
	{
		temp = 0;
		if (A[i*num + i] != 0)
		{
			for (int j = 0; j <= i - 1; j++)
			{
				temp = temp + L[i*num + j] * L[i*num + j];
			}
			L[i*num + i] = sqrt(A[i*num + i] - temp);
		}

		for (int j = i + 1; j < num; j++)
		{
			if (A[i*num + j] != 0 && L[i*num + i] != 0)
			{
				temp = 0;
				for (int k = 0; k <= i - 1; k++)
				{
					temp = temp + L[i*num + k] * L[j*num + k];
				}
				L[j*num + i] = (A[i*num + j] - temp) / L[i*num + i];
			}
		}
	}

	return L;
}

//void incompleteCholesky(int num, double* A)
//{
//	//double* L= new double[num*num];
//	for (int i = 0; i < num; i++)
//	{
//		for (int j = 0; j < num; j++)
//		{
//			//A[i*num + j] = 0;
//		}
//	}
//
//	double temp;
//	for (int i = 0; i < num; i++)
//	{
//		A[i*num + i] = sqrt(A[i*num + i]);
//
//		for (int j = i+1; j <= num; j++)
//		{
//			if (A[j*num + i]!=0)
//			{
//				A[j*num + i] = A[j*num + i]/A[i*num + i];
//			}
//		}
//
//		for (int j = i+1; j < num; j++)
//		{
//			for (int k = j; k <= num; k++)
//			{
//				if (A[k*num + j]!=0)
//				{
//					A[k*num + j] = A[k*num + j] - A[k*num+i]*A[j*num + i];
//				}
//			}
//		}
//		for (int j = 0; j < i; j++)
//		{
//			A[j*num + i] = 0;
//		}
//	}
//
//	//return L;
//}


double* PCG(int num, double* A, double* b)
{
	double tolerance = 1000 * DBL_EPSILON;
	double* rOld = new double[num];
	double* p = new double[num];
	double* rNew = new double[num];
	double* zOld = new double[num];
	double* zNew = new double[num];
	double* x = new double[num];

	double temp;
	double* M = new double[num*num];
	for (int j = 0; j < num; j++)
	{
		for (int i = 0; i < num; i++)
		{
			if (i == j)
			{
				M[j*num + j] = 1 / A[j*num + j];
			}
			else
			{
				M[i*num + j] = 0;
			}
		}
	}
	//double* L = incompleteCholesky(num,A);

	//for (int i = 0; i < num; i++)
	//{
	//	for (int j = 0; j < num; j++)
	//	{
	//		M[i*num + j]  = 0;
	//		for (int k = 0; k < num; k++)
	//		{
	//			M[i*num + j] = M[i*num + j] + L[i*num + k]*L[j*num + k];
	//		}
	//	}
	//}
	//ofstream asdf, zxcv;
	//asdf.open("D:\Data/L.txt");
	//zxcv.open("D:\Data/M.txt");
	//for (int i = 0; i < num; i++)
	//{
	//	for (int j = 0; j < num; j++)
	//	{
	//		asdf<<L[i*num +j]<<" ";
	//		zxcv<<M[i*num +j]<<" ";
	//	}
	//	asdf<<endl;
	//	zxcv<<endl;

	//}
	//asdf.close();
	//zxcv.close();


	for (int i = 0; i < num; i++)
	{
		rOld[i] = b[i];
		temp = 0;
		for (int j = 0; j < num; j++)
		{
			temp = temp + M[i*num + j] * rOld[j];
		}
		zOld[i] = temp;
		p[i] = zOld[i];
		x[i] = 0;

	}

	double alpha = 0;
	double beta = 0;
	double temp1 = 0, temp2 = 0;
	double residual;
	double residualOld = 0;
	double zrOld = 0;
	double zrNew = 0;

	for (int i = 0; i < num; i++)
	{
		residualOld = residualOld + rOld[i] * rOld[i];
		zrOld = zrOld + rOld[i] * zOld[i];
	}

	for (int k = 0; k < num * 2; k++)
	{
		temp2 = 0;
		for (int i = 0; i < num; i++)
		{
			for (int j = 0; j < num; j++)
			{
				temp2 = temp2 + p[i] * A[i*num + j] * p[j];
			}
		}
		alpha = zrOld / temp2;

		for (int i = 0; i < num; i++)
		{
			x[i] = x[i] + alpha*p[i];
			temp = 0;
			for (int j = 0; j < num; j++)
			{
				temp = temp + A[i*num + j] * p[j];
			}
			rNew[i] = rOld[i] - alpha*temp;
		}

		residual = 0;
		for (int i = 0; i < num; i++)
		{
			residual = residual + rNew[i] * rNew[i];
		}

		temp = sqrt(residual);
		cout << k << " " << temp << endl;
		if (temp < tolerance)
		{
			delete zNew, zOld, rNew, rOld, p;
			delete M;

			return x;
		}

		for (int i = 0; i < num; i++)
		{
			temp = 0;
			for (int j = 0; j < num; j++)
			{
				temp = temp + M[i*num + j] * rNew[j];
			}
			zNew[i] = temp;
		}

		zrNew = 0;
		for (int i = 0; i < num; i++)
		{
			zrNew = zrNew + rNew[i] * zNew[i];
		}
		beta = zrNew / zrOld;

		for (int i = 0; i < num; i++)
		{
			p[i] = zNew[i] + beta*p[i];
			rOld[i] = rNew[i];
			zOld[i] = zNew[i];
		}
		zrOld = zrNew;
		residualOld = residual;

	}

	delete[] zNew, zOld, rNew, rOld, p;
	delete[] M;

	return x;

}

