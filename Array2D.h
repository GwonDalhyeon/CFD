#pragma once


//#ifndef Array2D_H
//#define Array2D_H
#include "CommonDef.h"
#include "Grid2D.h"
#include "Vector2D.h"


template <class TT>
class Array2D
{
public:
	union
	{
		struct { int iStart, jStart, iEnd, jEnd; };
		struct { int ijStart[2], ijEnd[2]; };
	};

	int iRes, jRes;
	int ijRes;
	TT* values;

	Array2D();
	~Array2D();

	Array2D(const int& ipiRes);
	Array2D(const int& ipiRes, const int& ipjRes);
	Array2D(const int& ipiStart, const int& ipiEnd, const int& ipiRes);
	Array2D(const int& ipiStart, const int& ipiRes, const int& ipjStart, const int& ipjRes);
	Array2D(const Array2D<TT>& ipArray);

	Array2D(const Grid2D& ipGrid);

	void initialize(const int& iS, const int& iL, const int& jS, const int& jL);
	void initialize(const int& iS, const int& iE, const int& iL, const int& jS, const int& jE, const int& jL, const __int64& ijL);
	void initialValues();

	const int index(const Vector2D<int>& ipVector) const;

	const int index(const int& i, const int& j) const;

	inline TT& operator [](const int& i) const;

	inline TT& operator ()(const int& i)const;

	inline TT& operator ()(const int& i, const int& j) const;

	inline TT& operator ()(const Vector2D<int>& ipVector) const;

	inline void operator =(const TT& constant);

	inline void operator *=(const TT& constant);

	inline void operator +=(const TT& constant);

	inline void operator -=(const TT& constant);

	inline void operator /=(const TT& constant);

	inline void operator =(const Array2D<TT>& ipArray);

	Array2D<TT> operator + (const Array2D<TT>& ipArray) const;

	Array2D<TT> operator - (const Array2D<TT>& ipArray) const;

	Array2D<TT> operator * (const Array2D<TT>& ipArray) const;

	Array2D<TT> operator / (const Array2D<TT>& ipArray) const;

	Array2D<TT> operator + (const TT& constant) const;

	Array2D<TT> operator - (const TT& constant) const;

	Array2D<TT> operator * (const TT& constant) const;

	Array2D<TT> operator / (const TT& constant) const;

	static TT L1Norm(const Array2D<TT>& ipArray1);
	static TT L2Norm(const Array2D<TT>& ipArray1);
	static TT InnerProduct(const Array2D<TT>& ipArray1, const Array2D<TT>& ipArray2);

	// Write MATLAB Variable.
	inline void Variable(const char * varName);
	inline void Delete();
private:

};

//#endif // !Array2D





template <class TT>
Array2D<TT>::Array2D()
{
	values = nullptr;
}

template <class TT>
Array2D<TT>::~Array2D()
{
	if (values != nullptr)
	{
		delete[] values;
	}
}

template<class TT>
inline Array2D<TT>::Array2D(const int & ipiRes)
{
	if (values != nullptr)
	{
		values = nullptr;
	}
	initialize(0, ipiRes - 1, ipiRes, 0, 0, 0, ipiRes);

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const int & ipiRes, const int& ipjRes)
{
	if (values != nullptr)
	{
		values = nullptr;
	}

	initialize(0, ipiRes - 1, ipiRes, 0, ipjRes - 1, ipjRes, ipiRes*ipjRes);

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const int & ipiStart, const int & ipiEnd, const int & ipiRes)
{
	if (values != nullptr)
	{
		delete[] values;
	}

	initialize(ipiStart, ipiEnd, ipiRes, 0, 0, 0, ipiRes);

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const int & ipiStart, const int & ipiRes, const int & ipjStart, const int & ipjRes)
{
	if (values != nullptr)
	{
		values = nullptr;
		//delete[] values;
	}
	initialize(ipiStart, ipiStart + ipiRes - 1, ipiRes, ipjStart, ipjStart + ipjRes - 1, ipjRes, ipiRes*ipjRes);

	initialValues();
}

template<class TT>
inline Array2D<TT>::Array2D(const Array2D<TT>& ipArray)
{
	if (values != nullptr)
	{
		values = nullptr;
	}

	initialize(ipArray.iStart, ipArray.iEnd, ipArray.iRes, ipArray.jStart, ipArray.jEnd, ipArray.jRes, ipArray.ijRes);

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		values[i] = ipArray.values[i];
	}
}

template<class TT>
inline Array2D<TT>::Array2D(const Grid2D & ipGrid)
{
	//if (values != nullptr)
	//{
	//	delete[] values;
	//}

	initialize(ipGrid.iStart, ipGrid.iEnd, ipGrid.iRes, ipGrid.jStart, ipGrid.jEnd, ipGrid.jRes, ipGrid.iRes*ipGrid.jRes);

	memset(values, 0, ijRes * sizeof(TT));
//#pragma omp parallel for
//	for (int i = 0; i < ijRes; i++)
//	{
//		values[i] = 0;
//	}
}

template<class TT>
inline void Array2D<TT>::initialize(const int & iS, const int & iL, const int & jS, const int & jL)
{
	initialize(iS, iL, iS + iL - 1, jS, jL, jS + jL - 1, (__int64) iL*jL);
}

template<class TT>
inline void Array2D<TT>::initialize(const int & iS, const int & iE, const int & iL, const int & jS, const int & jE, const int & jL, const __int64& ijL)
{
	iStart = iS;
	iEnd = iE;
	iRes = iL;
	jStart = jS;
	jEnd = jE;
	jRes = jL;
	ijRes = ijL;

	assert(iRes > 0);
	assert(iEnd == iStart + iRes - 1);
	assert(jRes > 0);
	assert(jEnd == jStart + jRes - 1);

	values = new TT[ijRes];
	initialValues();
}

template<class TT>
inline void Array2D<TT>::initialValues()
{
	memset(values, 0, ijRes * sizeof(TT));
//#pragma omp parallel for
//	for (int i = 0; i < ijRes; i++)
//	{
//		values[i] = 0;
//	}
}

template<class TT>
const int Array2D<TT>::index(const Vector2D<int>& ipVector) const
{
	assert(ipVector[0] >= iStart && ipVector[0] <= iEnd);
	assert(ipVector[1] >= jStart && ipVector[1] <= jEnd);
	return (ipVector[0] - iStart) + (ipVector[1] - jStart)*jRes;
}

template<class TT>
const int Array2D<TT>::index(const int & i, const int & j) const
{
	assert(i >= iStart && i <= iEnd);
	assert(j >= jStart && j <= jEnd);
	return (i - iStart) + (j - jStart)*iRes;
}

template<class TT>
inline TT & Array2D<TT>::operator[](const int & i) const
{
	//assert(i >= iStart && i <=iEnd);
	return values[i];
}

template<class TT>
inline TT & Array2D<TT>::operator()(const int & i) const
{
	//assert(i >= iStart && i <=iEnd);
	return values[i];
}

template<class TT>
inline TT & Array2D<TT>::operator()(const int & i, const int & j) const
{
	assert(i >= iStart && i <= iEnd);
	assert(j >= jStart && j <= jEnd);

	return values[index(i, j)];
}

template<class TT>
inline TT & Array2D<TT>::operator()(const Vector2D<int>& ipVector) const
{
	assert(ipVector[0] >= iStart && ipVector[0] <= iEnd);
	assert(ipVector[1] >= jStart && ipVector[1] <= jEnd);

	return values[index(ipVector[0], ipVector[1])];
}

template<class TT>
inline void Array2D<TT>::operator=(const TT & constant)
{
#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		values[i] = constant;
	}
}

template<class TT>
inline void Array2D<TT>::operator*=(const TT & constant)
{
#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		values[i] *= constant;
	}
}

template<class TT>
inline void Array2D<TT>::operator+=(const TT & constant)
{
#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		values[i] += constant;
	}
}

template<class TT>
inline void Array2D<TT>::operator-=(const TT & constant)
{
#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		values[i] -= constant;
	}
}

template<class TT>
inline void Array2D<TT>::operator/=(const TT & constant)
{
	assert(constant != 0);

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		values[i] *= 1 / constant;
	}
}

template<class TT>
inline void Array2D<TT>::operator=(const Array2D<TT>& ipArray)
{
	if (values != nullptr)
	{
		delete[] values;
	}
	initialize(ipArray.iStart, ipArray.iEnd, ipArray.iRes, ipArray.jStart, ipArray.jEnd, ipArray.jRes, ipArray.ijRes);

	memcpy(values, ipArray.values, ijRes * sizeof(TT));
}

template<class TT>
Array2D<TT> Array2D<TT>::operator+(const Array2D<TT>& ipArray) const
{
	Array2D<TT> tempArray(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
			tempArray(i) = values[i] + ipArray(i);
	}
	return tempArray;
}

template<class TT>
Array2D<TT> Array2D<TT>::operator-(const Array2D<TT>& ipArray) const
{
	Array2D<TT> tempArray(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		tempArray(i) = values[i] - ipArray(i);
	}
	return tempArray;
}

template<class TT>
Array2D<TT> Array2D<TT>::operator*(const Array2D<TT>& ipArray) const
{
	Array2D<TT> tempArray(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		tempArray(i) = values[i] * ipArray(i);
	}
	return tempArray;
}

template<class TT>
Array2D<TT> Array2D<TT>::operator/(const Array2D<TT>& ipArray) const
{
	Array2D<TT> tempArray(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
			if (ipArray(i) != 0)
			{
				tempArray(i) = values[i] / ipArray(i);
			}
			else
			{
				tempArray(i) = 0;
			}
	}
	return tempArray;
}

template<class TT>
Array2D<TT> Array2D<TT>::operator+(const TT & constant) const
{
	Array2D<TT> tempArray = this;

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		tempArray(i) = values[i] + constant;
	}
	return tempArray;
}

template<class TT>
Array2D<TT> Array2D<TT>::operator-(const TT & constant) const
{
	Array2D<TT> tempArray = this;

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		tempArray(i) = values[i] - constant;
	}
	return tempArray;
}

template<class TT>
Array2D<TT> Array2D<TT>::operator*(const TT & constant) const
{
	Array2D<TT> tempArray = this;

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		tempArray(i) = values[i] - constant;
	}
	return tempArray;
}

template<class TT>
Array2D<TT> Array2D<TT>::operator/(const TT & constant) const
{
	assert(constant != 0);
	Array2D<TT> tempArray = this;

#pragma omp parallel for
	for (int i = 0; i < ijRes; i++)
	{
		tempArray(i) = values[i] / constant;
	}
	return tempArray;
}


template<class TT>
inline TT Array2D<TT>::L1Norm(const Array2D<TT>& ipArray1)
{
	TT l1 = 0;
#pragma omp parallel for reduction(+ : l1)
	for (int i = 0; i < ijRes; i++)
	{
		l1 += abs(ipArray1[i]);
	}
	return l1;
}

template<class TT>
inline TT Array2D<TT>::L2Norm(const Array2D<TT>& ipArray1)
{

	TT l2 = 0;
#pragma omp parallel for reduction(+ : l2)
	for (int i = 0; i < ijRes; i++)
	{
		l2 += ipArray1[i]* ipArray1[i];
	}
	return sqrt(l2);
}

template<class TT>
inline TT Array2D<TT>::InnerProduct(const Array2D<TT>& ipArray1, const Array2D<TT>& ipArray2)
{
	TT inner = 0;

#pragma omp parallel for reduction(+ : inner)
	for (int i = 0; i < ijRes; i++)
	{
		inner += ipArray1[i] * ipArray2[i];
	}
	return TT();
}

template<class TT>
inline void Array2D<TT>::Variable(const char * varName)
{
	MATLAB.Variable(varName, iRes, jRes, values);
}

template<class TT>
inline void Array2D<TT>::Delete()
{
	if (values != nullptr)
	{
		delete[] values;
	}
}



