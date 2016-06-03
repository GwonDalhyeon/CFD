#pragma once


//#ifndef Array1D_H
//#define Array1D_H
#include "CommonDef.h"
#include "Grid1D.h"
#include "VectorND.h"


template <class TT>
class Array1D
{
public:
	int iStart, iEnd;

	int iRes;
	TT* values;

	Array1D();
	~Array1D();

	Array1D(const int& ipiRes);
	//Array1D(const int& ipiStart, const int& ipiEnd, const int& ipiRes);
	Array1D(const int& ipiStart, const int& ipiRes);
	Array1D(const Array1D<TT>& ipArray);

	Array1D(const Grid1D& ipGrid);

	void initialize(const int& iS, const int& iE, const int& iLL);
	void initialValues();

	const int index(const int& i) const;

	inline TT& operator [](const int& i) const;

	inline TT& operator ()(const int& i)const;

	inline void operator =(const double& constant);

	inline void operator *=(const double& constant);

	inline void operator +=(const double& constant);

	inline void operator -=(const double& constant);

	inline void operator /=(const double& constant);

	inline void operator =(const Array1D<TT>& ipArray);

	Array1D<TT> operator + (const Array1D<TT>& ipArray) const;

	Array1D<TT> operator - (const Array1D<TT>& ipArray) const;

	Array1D<TT> operator * (const Array1D<TT>& ipArray) const;

	Array1D<TT> operator / (const Array1D<TT>& ipArray) const;

	Array1D<TT> operator + (const TT& constant) const;

	Array1D<TT> operator - (const TT& constant) const;

	Array1D<TT> operator * (const TT& constant) const;

	Array1D<TT> operator / (const TT& constant) const;

	static TT L1Norm(const Array1D<TT>& ipArray1);
	static TT L2Norm(const Array1D<TT>& ipArray1);
	static TT InnerProduct(const Array1D<TT>& ipArray1, const Array1D<TT>& ipArray2);

	// Write MATLAB Variable.
	inline void Variable(const char * varName);

private:

};

//#endif // !Array1D





template <class TT>
Array1D<TT>::Array1D()
{
	values = nullptr;
}

template <class TT>
Array1D<TT>::~Array1D()
{
	if (values != nullptr)
	{
		delete[] values;
	}
}

template<class TT>
inline Array1D<TT>::Array1D(const int & ipiRes)
{
	if (values != nullptr)
	{
		values = nullptr;
	}
	initialize(0, ipiRes - 1, ipiRes);

	assert(iRes > 0);

	values = new TT[iRes];

	initialValues();
}

template<class TT>
inline Array1D<TT>::Array1D(const int & ipiStart, const int & ipiRes)
{
	if (values != nullptr)
	{
		values = nullptr;
		//delete[] values;
	}
	initialize(ipiStart, ipiStart + ipiRes - 1, ipiRes);

	assert(iRes > 0 && iEnd == iStart + iRes - 1);
	values = new TT[iRes];

	initialValues();
}

template<class TT>
inline Array1D<TT>::Array1D(const Array1D<TT>& ipArray)
{
	if (values != nullptr)
	{
		values = nullptr;
	}

	initialize(ipArray.iStart, ipArray.iEnd, ipArray.iRes);

	values = new TT[iRes];

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		values[i] = ipArray.values[i];
	}
}

template<class TT>
inline Array1D<TT>::Array1D(const Grid1D & ipGrid)
{
	//if (values != nullptr)
	//{
	//	delete[] values;
	//}

	initialize(ipGrid.iStart, ipGrid.iEnd, ipGrid.iRes);

	values = new TT[iRes];

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		values[i] = 0;
	}
}

template<class TT>
inline void Array1D<TT>::initialize(const int & iS, const int & iE, const int & iL)
{
	iStart = iS;
	iEnd = iE;
	iRes = iL;
}

template<class TT>
inline void Array1D<TT>::initialValues()
{
#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		values[i] = 0;
	}
}

template<class TT>
const int Array1D<TT>::index(const int & i) const
{
	assert(i >= iStart && i <= iEnd);
	return (i - iStart);
}

template<class TT>
inline TT & Array1D<TT>::operator[](const int & i) const
{
	//assert(i >= iStart && i <=iEnd);
	return values[i];
}

template<class TT>
inline TT & Array1D<TT>::operator()(const int & i) const
{
	//assert(i >= iStart && i <=iEnd);
	return values[i];
}


template<class TT>
inline void Array1D<TT>::operator=(const double & constant)
{
#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		values[i] = constant;
	}
}

template<class TT>
inline void Array1D<TT>::operator*=(const double & constant)
{
#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		values[i] *= constant;
	}
}

template<class TT>
inline void Array1D<TT>::operator+=(const double & constant)
{
#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		values[i] += constant;
	}
}

template<class TT>
inline void Array1D<TT>::operator-=(const double & constant)
{
#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		values[i] -= constant;
	}
}

template<class TT>
inline void Array1D<TT>::operator/=(const double & constant)
{
	assert(constant != 0);

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		values[i] *= 1 / constant;
	}
}

template<class TT>
inline void Array1D<TT>::operator=(const Array1D<TT>& ipArray)
{
	if (values != nullptr)
	{
		delete[] values;
	}
	initialize(ipArray.iStart, ipArray.iEnd, ipArray.iRes);

	assert(iRes > 0);

	values = new TT[iRes];
#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		values[i] = ipArray.values[i];
	}
}

template<class TT>
Array1D<TT> Array1D<TT>::operator+(const Array1D<TT>& ipArray) const
{
	Array1D<TT> tempArray(ipArray.iStart, ipArray.iRes);

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		tempArray.values[i] = values[i] + ipArray.values[i];
	}
	return tempArray;
}

template<class TT>
Array1D<TT> Array1D<TT>::operator-(const Array1D<TT>& ipArray) const
{
	Array1D<TT> tempArray(ipArray.iStart, ipArray.iRes);

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		tempArray.values[i] = values[i] - ipArray.values[i];
	}
	return tempArray;
}

template<class TT>
Array1D<TT> Array1D<TT>::operator*(const Array1D<TT>& ipArray) const
{
	Array1D<TT> tempArray(ipArray.iStart, ipArray.iRes);

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		tempArray.values[i] = values[i] * ipArray.values[i];
	}
	return tempArray;
}

template<class TT>
Array1D<TT> Array1D<TT>::operator/(const Array1D<TT>& ipArray) const
{
	Array1D<TT> tempArray(ipArray.iStart, ipArray.iRes);

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		if (ipArray(i) != 0)
		{
			tempArray.values[i] = values[i] / ipArray.values[i];
		}
		else
		{
			tempArray.values[i] = 0;
		}
	}
	return tempArray;
}

template<class TT>
Array1D<TT> Array1D<TT>::operator+(const TT & constant) const
{
	Array1D<TT> tempArray = this;

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		tempArray.values[i] = values[i] + constant;
	}
	return tempArray;
}

template<class TT>
Array1D<TT> Array1D<TT>::operator-(const TT & constant) const
{
	Array1D<TT> tempArray = this;

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		tempArray.values[i] = values[i] - constant;
	}
	return tempArray;
}

template<class TT>
Array1D<TT> Array1D<TT>::operator*(const TT & constant) const
{
	Array1D<TT> tempArray = this;

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		tempArray.values[i] = values[i] - constant;
	}
	return tempArray;
}

template<class TT>
Array1D<TT> Array1D<TT>::operator/(const TT & constant) const
{
	assert(constant != 0);
	Array1D<TT> tempArray = this;

#pragma omp parallel for
	for (int i = 0; i < iRes; i++)
	{
		tempArray.values[i] = values[i] / constant;
	}
	return tempArray;
}


template<class TT>
inline TT Array1D<TT>::L1Norm(const Array1D<TT>& ipArray1)
{
	TT l1 = 0;
#pragma omp parallel for reduction(+ : l1)
	for (int i = 0; i < iRes; i++)
	{
		l1 += abs(ipArray1.values[i]);
	}
	return l1;
}

template<class TT>
inline TT Array1D<TT>::L2Norm(const Array1D<TT>& ipArray1)
{

	TT l2 = 0;
#pragma omp parallel for reduction(+ : l2)
	for (int i = 0; i < iRes; i++)
	{
		l2 += ipArray1.values[i] * ipArray1.values[i];
	}
	return sqrt(l2);
}

template<class TT>
inline TT Array1D<TT>::InnerProduct(const Array1D<TT>& ipArray1, const Array1D<TT>& ipArray2)
{
	TT inner = 0;

#pragma omp parallel for reduction(+ : inner)
	for (int i = 0; i < iRes; i++)
	{
		inner += ipArray1.values[i] * ipArray2.values[i];
	}
	return TT();
}

template<class TT>
inline void Array1D<TT>::Variable(const char * varName)
{
	MATLAB.Variable(varName, iRes, 1, values);
}



