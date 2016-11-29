#pragma once

#include "CommonDef.h"
#include "VectorND.h"
#include "Array2D.h"

template <class TT>
class CSR
{
public:
	VectorND<TT> values;
	VectorND<int> columns;
	VectorND<int> indPrt;

	int colNum;
	int rowNum;
	int valueNum;

	int value_ix = 0;
	int prev_row = -1;

	CSR();
	~CSR();

	CSR(const int& rNum, const int& cNum);
	CSR(const int& rNum, const int& cNum, const int& vNum);
	CSR(const Array2D<TT>& ipArray);
	//CSR(const Field2D<TT>& ipField);

	inline void operator = (const CSR<TT>& ipCSR);
	inline TT& operator()(const int& row_input, const int& column_input) const;

	inline void AssignValue(const int& row_input, const int& column_input, const TT& values_input);

	inline void Multiply(const VectorND<TT>& x, VectorND<TT>& b) const;
	inline void ComputeResidual(const VectorND<TT>& x, const VectorND<TT>& b, VectorND<TT>& residual) const;

private:

};

template<class TT>
inline std::ostream& operator<<(std::ostream& output, const CSR<TT>& ipCSR);



template <class TT>
CSR<TT>::CSR()
{
}

template <class TT>
CSR<TT>::~CSR()
{
}

template<class TT>
inline CSR<TT>::CSR(const int & rNum, const int & cNum)
{
	rowNum = rNum;
	colNum = cNum;

	indPrt = VectorND<int>(rowNum + 1);

	value_ix = 0;
	prev_row = -1;
}

template<class TT>
inline CSR<TT>::CSR(const int & rNum, const int & cNum, const int & vNum)
{
	rowNum = rNum;
	colNum = cNum;
	valueNum = vNum;

	indPrt = VectorND<int>(rowNum + 1);

	value_ix = 0;
	prev_row = -1;
}

template<class TT>
inline CSR<TT>::CSR(const Array2D<TT>& ipArray)
{
	clock_t before;
	double  result;
	before = clock();
	cout << "Start : CSR" << endl;

	rowNum = ipArray.iRes;
	colNum = ipArray.jRes;
	indPrt = VectorND<int>(rowNum + 1);

	value_ix = 0;
	prev_row = -1;
	int qwetr = int(floor(sqrt(double(rowNum*colNum)))) * 10;
	TT* tempVal = new TT[int(floor(sqrt(double(rowNum*colNum)))) * 10];
	int* tempCol = new int[int(floor(sqrt(double(rowNum*colNum)))) * 10];

#pragma omp parallel for
	for (int i = 0; i < rowNum + 1; i++)
	{
		indPrt[i] = -1;

	}
	int tempIndex = 0;


	for (int i = 0; i < rowNum; i++)
	{
		for (int j = 0; j < colNum; j++)
		{
			if (ipArray(i + ipArray.iStart, j + ipArray.jStart) != 0)
			{
				tempVal[tempIndex] = ipArray(i + ipArray.iStart, j + ipArray.jStart);
				tempCol[tempIndex] = j;
				if (indPrt[i]<0)
				{
					indPrt[i] = tempIndex;
				}
				tempIndex += 1;
			}
		}
	}
	valueNum = tempIndex;
	indPrt[rowNum] = tempIndex;

	values = VectorND<TT>(valueNum);
	columns = VectorND<int>(valueNum);

#pragma omp parallel for
	for (int i = 0; i < valueNum; i++)
	{
		values[i] = tempVal[i];
		columns[i] = tempCol[i];
	}
	delete[] tempVal, tempCol;

	result = (double)(clock() - before) / CLOCKS_PER_SEC;
	cout << "time : " << result << "\n";
	cout << "End : Make sparse matrix." << endl;
	cout << endl;
}


template<class TT>
inline void CSR<TT>::operator=(const CSR<TT>& ipCSR)
{
	values = ipCSR.values;
	columns = ipCSR.columns;
	indPrt = ipCSR.indPrt;

	colNum = ipCSR.colNum;
	rowNum = ipCSR.rowNum;
	valueNum = ipCSR.valueNum;

	value_ix = ipCSR.value_ix;
	prev_row = ipCSR.prev_row;
}

template<class TT>
inline TT & CSR<TT>::operator()(const int & row_input, const int & column_input) const
{
	static TT values_for_test(0);

	for (int vix = indPrt[row_input]; vix < indPrt[row_input + 1]; vix++)
	{
		if (columns[vix] == column_input)
		{
			return values[vix];
		}
	}

	return values_for_test;
}

template<class TT>
inline void CSR<TT>::AssignValue(const int & row_input, const int & column_input, const TT & values_input)
{
	values[value_ix] = values_input;

	if (row_input != prev_row)
	{
		indPrt[row_input] = value_ix;
		prev_row = row_input;
	}

	columns[value_ix] = column_input;

	value_ix++;
}

template<class TT>
inline void CSR<TT>::Multiply(const VectorND<TT>& x, VectorND<TT>& b) const
{
	assert(rowNum == x.iLength);
	assert(x.iLength == b.iLength);

	//TT *bval(b.values), *xval(x.values);
	TT v = 0;
//#pragma omp parallel for private (v)
//	for (int row = 0; row < rowNum; row++)
//	{
//		v = 0;
//		for (int vix = indPrt[row]; vix < indPrt[row + 1]; vix++)
//		{
//			v += values[vix] * x.values[columns.values[vix]];
//		}
//
//		b.values[row] = v;
//	}

	int num = rowNum;
	int j;
#pragma omp parallel for private (v, j)
	for (int i = 0; i < num; i++)
	{
		v = 0;
//#pragma omp parallel for private (j) reduction(+:v) 
		for (int n = indPrt[i]; n < indPrt[i + 1]; n++)
		{
			j = int(columns[n]);
			v += values[n] * x[j];
		}
		b[i] = v;
	}
}

template<class TT>
inline void CSR<TT>::ComputeResidual(const VectorND<TT>& x, const VectorND<TT>& b, VectorND<TT>& residual) const
{
	assert(rowNum == x.iLength);
	assert(x.iLength == b.iLength);
	assert(residual.iLength == rowNum);

	//TT *bval(b.values), *xval(x.values), *rval(residual.values);
	TT v = 0;
#pragma omp parallel for private (v)
	for (int row = 0; row < rowNum; row++)
	{
		// Compute A*x
		v = 0;
		for (int vix = indPrt[row]; vix < indPrt[row + 1]; vix++)
		{
			v += values[vix] * x[columns[vix]];
		}

		// residual = b - A*x
		residual[row] = b[row] - v;
	}
}



template<class TT>
inline std::ostream& operator<<(std::ostream& output, const CSR<TT>& ipCSR)
{
	output << "CSR" << endl;
	output << "- colNum = " << ipCSR.colNum << endl;
	output << "- rowNum = " << ipCSR.rowNum << endl;
	output << "- valueNum = " << ipCSR.valueNum << endl;

	return output;
}
