#pragma once

#include "CommonDef.h"
#include "VectorND.h"
#include "Field2D.h"

template <class TT>
class CSR
{
public:
	VectorND<TT> values;   // the values of the nonzero elements
	VectorND<int> columns; // the column indices of the elements in the val vector
	VectorND<int> indPrt;  // the locations in the val vector that start a row

	int rowNum;
	int valueNum;

	int value_ix;
	int prev_row;

	CSR();
	~CSR();

	CSR(const int& _rNum, const int& _valueNum);
	CSR(const int& rNum, const int& cNum, const int& vNum);
	CSR(const CSR<TT>& ipCSR);
	CSR(const Array2D<TT>& ipArray);
	CSR(const Field2D<TT>& ipField);
	//CSR(const Field2D<TT>& ipField);

	inline void Initialize(const int& _rNum, const int& _valueNum);
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
inline CSR<TT>::CSR(const int & _rNum, const int & _valueNum)
{
	Initialize(_rNum, _valueNum);
}

template<class TT>
inline CSR<TT>::CSR(const int & rNum, const int & cNum, const int & vNum)
{
	rowNum = rNum;
	valueNum = vNum;

	indPrt = VectorND<int>(rowNum + 1);

	value_ix = 0;
	prev_row = -1;
}

template<class TT>
inline CSR<TT>::CSR(const CSR<TT>& ipCSR)
{
	rowNum = ipCSR.rowNum;
	valueNum = ipCSR.valueNum;

	values = ipCSR.values;
	columns = ipCSR.columns;
	indPrt = ipCSR.indPrt;

	value_ix = ipCSR.value_ix;
	prev_row = ipCSR.prev_row;
}

template<class TT>
inline CSR<TT>::CSR(const Array2D<TT>& ipArray)
{
	clock_t before;
	double  result;
	before = clock();
	cout << "Start : CSR" << endl;

	rowNum = ipArray.iRes;
	indPrt = VectorND<int>(rowNum + 1);

	value_ix = 0;
	prev_row = -1;

	TT* tempVal = new TT[int(floor(sqrt(double(rowNum*rowNum)))) * 10];
	int* tempCol = new int[int(floor(sqrt(double(rowNum*rowNum)))) * 10];

#pragma omp parallel for
	for (int i = 0; i < rowNum + 1; i++)
	{
		indPrt[i] = -1;

	}
	int tempIndex = 0;


	for (int i = 0; i < rowNum; i++)
	{
		for (int j = 0; j < rowNum; j++)
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
	tempVal = 0;
	tempCol = 0;

	cout << "time : " << (double)(clock() - before) / CLOCKS_PER_SEC << endl;
	cout << "End : Make sparse matrix." << endl;
	cout << endl;
}

template<class TT>
inline void CSR<TT>::Initialize(const int & _rNum, const int & _valueNum)
{
	rowNum = _rNum;
	valueNum = _valueNum;

	values = VectorND<TT>(valueNum);
	columns = VectorND<int>(valueNum);
	indPrt = VectorND<int>(rowNum + 1);

	value_ix = 0;
	prev_row = -1;

	indPrt[rowNum] = valueNum;
}

template<class TT>
inline void CSR<TT>::operator=(const CSR<TT>& ipCSR)
{
	values = ipCSR.values;
	columns = ipCSR.columns;
	indPrt = ipCSR.indPrt;

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


template<class TT>
static void IncompleteCholeskyDecomposition(const int& i_res_input, const int& j_res_input, const CSR<TT>& A, CSR<TT>& L, const Field2D<int>& bc_input)
{
	onst int N = A.N;
	const int nz = A.nz;

	Field2D<int>& index_field = bc_input;

	//const int i_start(index_field.iStartI), i_end(index_field.iEndI), j_start(index_field.jStartI), j_end(index_field.jEndI);

	//int start_ix(0);


	/*ofstream fout;
	fout.open("index_field");
	for (int j = 0; j < index_field.grid.j_res; j++)
	{
	for (int i = 0; i < index_field.iResI; i++)
	{
	fout << index_field(i, j) << " ";
	}
	fout << endl;
	}
	fout.close();*/

	TT sum, coef;
	int number(0);
//#pragma omp parallel for private (sum, coef)
	for (int i = index_field.iStart; i <= index_field.iEnd; i++)
	{
		for (int j = index_field.jStart; j <= index_field.jEnd; j++)
		{

			sum = 0, coef = 0;
			
			if (index_field(i, j - 1) > -1)
			{
				if (index_field(i, j) == index_field.iResI)
				{
					coef = 1 / L.values[L.indPrt[index_field(i, j) - index_field.iResI]] * (A(index_field(i, j - 1), index_field(i, j)));
				}
				else if (((index_field(i, j) > index_field.iResI) && (index_field(i, j) < 2 * index_field.iResI)) || (index_field(i, j) % index_field.iResI == 0))
				{
					coef = 1 / L.values[L.indPrt[index_field(i, j) - index_field.iResI] + 1] * (A(index_field(i, j - 1), index_field(i, j)));
				}
				else
				{
					if (index_field(i + 1, j) == BC_PER)
					{
						coef = 1 / L.values[L.indPrt[index_field(i, j) - index_field.iResI] + 3] * (A(index_field(i, j - 1), index_field(i, j)));
					}
					else
					{
						coef = 1 / L.values[L.indPrt[index_field(i, j) - index_field.iResI] + 2] * (A(index_field(i, j - 1), index_field(i, j)));
					}
				}

				L.AssignValue(index_field(i, j), index_field(i, j - 1), coef);
				number += 1;
			}

			if (index_field(i + 1, j) == BC_PER)
			{
				coef = 1 / L.values[L.indPrt[index_field(i, j) - 2] - 1] * (A(index_field(i, j) - (index_field.iResI - 1), index_field(i, j)));

				L.AssignValue(index_field(i, j), index_field(i, j) - (index_field.iResI - 1), coef);
				number += 1;
			}

			if (index_field(i - 1, j) > -1)
			{
				if (index_field(i + 1, j) == BC_PER)
				{
					coef = 1 / L.values[L.indPrt[index_field(i, j) - 1]] * (A(index_field(i - 1, j), index_field(i, j)));
				}
				else
				{
					if (index_field(i, j) == 1)
					{
						coef = 1 / L.values[L.indPrt[index_field(i, j) - 1]] * (A(index_field(i - 1, j), index_field(i, j)));
					}
					else if (index_field(i, j) < index_field.iResI)
					{
						coef = 1 / L.values[L.indPrt[index_field(i, j) - 1] + 1] * (A(index_field(i - 1, j), index_field(i, j)));
					}
					else
					{
						coef = 1 / L.values[L.indPrt[index_field(i, j)] - 1] * (A(index_field(i - 1, j), index_field(i, j)));
					}
				}

				L.AssignValue(index_field(i, j), index_field(i - 1, j), coef);
				number += 1;
			}

			if (index_field(i, j) == 0)
			{
				coef = sqrt(A(index_field(i, j), index_field(i, j)));
			}
			else
			{
				sum = (TT)0;

				for (int k = L.indPrt[index_field(i, j)]; k < number; k++)
				{
					sum += L.values[k]* L.values[k];
				}
				coef = sqrt(A(index_field(i, j), index_field(i, j)) - sum);
			}
			L.AssignValue(index_field(i, j), index_field(i, j), coef);

			number += 1;
		}
	}
}