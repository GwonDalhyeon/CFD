#pragma once



//#ifndef CSR_H
//#define CSR_H
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

	CSR();
	~CSR();

	CSR(const Array2D<TT>& ipArray);

	inline void operator = (const CSR<TT>& ipCSR);

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
inline CSR<TT>::CSR(const Array2D<TT>& ipArray)
{
	rowNum = ipArray.iRes;
	colNum = ipArray.jRes;
	indPrt = VectorND<int>(rowNum + 1);

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
				tempIndex = tempIndex + 1;
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

//#endif // !CSR_H

