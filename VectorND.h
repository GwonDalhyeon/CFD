#pragma once

//#ifndef VectorND_H
//#define VectorND_H
#include"CommonDef.h"
#include "ToMATLAB.h";

template <class TT>
class VectorND
{
public:

	TT* values;
	int iStart;
	int iEnd;
	int iLength;

	VectorND();
	~VectorND();
	VectorND(const int & ipLength);
	VectorND(const int & ipStart, const int & ipLength);
	VectorND(const VectorND<TT>& ipVector);
	VectorND(const int & ipLength, const TT* ipValues);
	VectorND(const int & ipStart, const int & ipLength, const TT* ipValues);

	const int index(const int& i) const;

	inline TT& operator [](const int& i) const;

	inline TT& operator ()(const int& i) const;

	inline void operator = (const TT& constant);

	inline void operator = (const VectorND<TT>& ipVector);

	inline void operator += (const VectorND<TT>& ipVector);

	inline void operator -= (const VectorND<TT>& ipVector);

	inline void operator /= (const VectorND<TT>& ipVector);

	inline void operator *= (const VectorND<TT>& ipVector);

	inline void operator += (const TT& constant);

	inline void operator -= (const TT& constant);

	inline void operator *= (const TT& constant);

	inline void operator /= (const TT& constant);

	VectorND<TT> operator + (const VectorND<TT>& ipVector);

	VectorND<TT> operator - (const VectorND<TT>& ipVector);

	VectorND<TT> operator * (const VectorND<TT>& ipVector);

	VectorND<TT> operator / (const VectorND<TT>& ipVector);

	VectorND<TT> operator + (const TT& constant);

	VectorND<TT> operator - (const TT& constant);

	VectorND<TT> operator * (const TT& constant);

	VectorND<TT> operator / (const TT& constant);

	inline void WriteFile(const string& fileName);

	TT magnitude();
	TT magnitude2();

	void normalize();

	inline void Variable(const char * varName);
	inline void Variable(const char * varName, const int & dim);

private:

};


template<class TT>
inline VectorND<TT> operator + (const TT& constant, const VectorND<TT>& ipVector);

template<class TT>
inline VectorND<TT> operator - (const TT& constant, const VectorND<TT>& ipVector);

template<class TT>
inline VectorND<TT> operator *(const TT& constant, const VectorND<TT>& ipVector);

template<class TT>
inline VectorND<TT> operator / (const TT& constant, const VectorND<TT>& ipVector);

template<class TT>
inline TT dotProduct(const VectorND<TT>& ipVector1, const VectorND<TT>& ipVector2);

template<class TT>
inline std::ostream& operator << (std::ostream& output, const VectorND<TT>& ipVector);

template<class TT>
VectorND<TT> normalize(const VectorND<TT>& ipVector);


//#endif // !VectorND


template <class TT>
VectorND<TT>::VectorND()
{
	values = nullptr;
}

template <class TT>
VectorND<TT>::~VectorND()
{
	if (iLength > 0 && values != nullptr)
	{
		delete[] values;
	}
}

template<class TT>
inline VectorND<TT>::VectorND(const int & ipLength)
{
	values = new TT[ipLength];
	iStart = 0;
	iLength = ipLength;
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] = 0;
	}
	iEnd = iStart + ipLength - 1;
}

template<class TT>
inline VectorND<TT>::VectorND(const int & ipStart, const int & ipLength)
{
	values = new TT[ipLength];
	iStart = ipStart;
	iLength = ipLength;
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] = 0;
	}
	iEnd = iStart + ipLength - 1;
}

template<class TT>
inline VectorND<TT>::VectorND(const VectorND<TT>& ipVector)
{
	if (values != nullptr)
	{
		values = nullptr;
	}
	assert(ipVector.iLength > 0);
	iStart = ipVector.iStart;
	iLength = ipVector.iLength;
	iEnd = ipVector.iEnd;
	values = new TT[iLength];

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] = ipVector.values[i];
	}
}

template<class TT>
inline VectorND<TT>::VectorND(const int & ipLength, const TT * ipValues)
{
	if (values != nullptr)
	{
		values = nullptr;
	}
	assert(ipVector.iLength > 0);
	values = new TT[ipLength];
	iStart = 0;
	iLength = ipLength;
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] = ipValues[i];
	}
	iEnd = iStart + ipLength - 1;
}

template<class TT>
inline VectorND<TT>::VectorND(const int & ipStart, const int & ipLength, const TT * ipValues)
{
	if (values != nullptr)
	{
		values = nullptr;
	}
	assert(ipVector.iLength > 0);
	values = new TT[ipLength];
	iStart = ipStart;
	iLength = ipLength;
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] = ipValues[i];
	}
	iEnd = iStart + ipLength - 1;
}

template<class TT>
inline const int VectorND<TT>::index(const int & i) const
{
	assert(i >= iStart || i <= iEnd);
	return (i - iStart);
}

template<class TT>
inline TT & VectorND<TT>::operator[](const int & i) const
{
	assert(i >= iStart || i <= iEnd);
	return values[index(i)];
}

template<class TT>
inline TT & VectorND<TT>::operator()(const int & i) const
{
	assert(i >= iStart || i <= iEnd);
	return values[index(i)];
}

template<class TT>
inline void VectorND<TT>::operator=(const TT & constant)
{
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] = constant;
		//values[i] = constant;
	}
}

template<class TT>
inline void VectorND<TT>::operator=(const VectorND<TT>& ipVector)
{
	iStart = ipVector.iStart;
	iLength = ipVector.iLength;
	iEnd = ipVector.iEnd;

	if (values != nullptr)
	{
		delete[] values;
	}

	values = new TT[iLength];

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] = ipVector.values[i];
		//values[i] = ipVector(i);
	}
}

template<class TT>
inline void VectorND<TT>::operator+=(const VectorND<TT>& ipVector)
{
	assert(iLength == ipVector.iLength);
	assert(iStart == ipVector.iStart);
	assert(iEnd == ipVector.iEnd);
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] += ipVector.values[i];
	}
}

template<class TT>
inline void VectorND<TT>::operator-=(const VectorND<TT>& ipVector)
{
	assert(iLength == ipVector.iLength);
	assert(iStart == ipVector.iStart);
	assert(iEnd == ipVector.iEnd);
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] -= ipVector.values[i];
	}
}

template<class TT>
inline void VectorND<TT>::operator/=(const VectorND<TT>& ipVector)
{
	assert(iLength == ipVector.iLength);
	assert(iStart == ipVector.iStart);
	assert(iEnd == ipVector.iEnd);
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] /= ipVector.values[i];
	}
}

template<class TT>
inline void VectorND<TT>::operator*=(const VectorND<TT>& ipVector)
{
	assert(iLength == ipVector.iLength);
	assert(iStart == ipVector.iStart);
	assert(iEnd == ipVector.iEnd);
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] *= ipVector.values[i];
	}
}

template<class TT>
inline void VectorND<TT>::operator+=(const TT & constant)
{
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] += constant;
	}
}

template<class TT>
inline void VectorND<TT>::operator-=(const TT & constant)
{
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] -= constant;
	}
}

template<class TT>
inline void VectorND<TT>::operator*=(const TT & constant)
{
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] *= constant;
	}
}

template<class TT>
inline void VectorND<TT>::operator/=(const TT & constant)
{
	assert(constant > 0 || constant < 0);
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		values[i] /= constant;
	}
}

template<class TT>
VectorND<TT> VectorND<TT>::operator+(const VectorND<TT>& ipVector)
{
	assert(iLength == ipVector.iLength);
	assert(iStart == ipVector.iStart);
	assert(iEnd == ipVector.iEnd);

	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = values[i] + ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
VectorND<TT> VectorND<TT>::operator-(const VectorND<TT>& ipVector)
{
	assert(iLength == ipVector.iLength);
	assert(iStart == ipVector.iStart);
	assert(iEnd == ipVector.iEnd);

	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = values[i] - ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
VectorND<TT> VectorND<TT>::operator*(const VectorND<TT>& ipVector)
{
	assert(iLength == ipVector.iLength);
	assert(iStart == ipVector.iStart);
	assert(iEnd == ipVector.iEnd);

	VectorND<TT> returnVT(ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = values[i] * ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
VectorND<TT> VectorND<TT>::operator/(const VectorND<TT>& ipVector)
{
	assert(iLength == ipVector.iLength);
	assert(iStart == ipVector.iStart);
	assert(iEnd == ipVector.iEnd);

	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		assert(ipVector.values[i]>0 && ipVector.values[i]<0);
		returnVT.values[i] = values[i] / ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
VectorND<TT> VectorND<TT>::operator+(const TT & constant)
{
	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = values[i] + constant;
	}
	return returnVT;
}

template<class TT>
VectorND<TT> VectorND<TT>::operator-(const TT & constant)

{
	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = values[i] - constant;
	}
	return returnVT;
}

template<class TT>
VectorND<TT> VectorND<TT>::operator*(const TT & constant)
{
	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = values[i] * constant;
	}
	return returnVT;
}

template<class TT>
VectorND<TT> VectorND<TT>::operator/(const TT & constant)
{
	assert(constant != 0);	
	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = values[i] / constant;
	}
	return returnVT;
}

template<class TT>
inline void VectorND<TT>::WriteFile(const string & fileName)
{
	ofstream solutionFile;
	solutionFile.open("D:\\Data/" + fileName + ".txt", ios::binary);
	for (int i = 0; i < iLength; i++)
	{
		solutionFile << values[i] << endl;
	}
	solutionFile.close();
}



template<class TT>
inline TT VectorND<TT>::magnitude2()
{
	TT mag2 = 0;
#pragma omp parallel for reduction(+:mag2)
	for (int i = 0; i < iLength; i++)
	{
		mag2 = mag2 + values[i] * values[i];
	}
	return TT(mag2);
}

template<class TT>
inline TT VectorND<TT>::magnitude()
{
	return TT(sqrt(magnitude2()));
}

template<class TT>
inline void VectorND<TT>::normalize()
{
	*this /= magnitude();
}

template<class TT>
inline void VectorND<TT>::Variable(const char * varName)
{
	MATLAB.Variable(varName, iLength, 1, values);
}



template<class TT>
inline VectorND<TT> operator+(const TT & constant, const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = constant + ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
inline VectorND<TT> operator-(const TT & constant, const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = constant - ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
inline VectorND<TT> operator*(const TT & constant, const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);
#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		returnVT.values[i] = constant * ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
inline VectorND<TT> operator/(const TT & constant, const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVT(ipVector.iStart, ipVector.iLength);

#pragma omp parallel for
	for (int i = 0; i < iLength; i++)
	{
		assert(ipVector.values[i] != 0);
		returnVT.values[i] = constant / ipVector.values[i];
	}
	return returnVT;
}

template<class TT>
inline TT dotProduct(const VectorND<TT>& ipVector1, const VectorND<TT>& ipVector2)
{
	assert(ipVector1.iLength == ipVector2.iLength);

	double dotPro = 0;
#pragma omp parallel for reduction(+:dotPro)
	for (int i = 0; i < iLength; i++)
	{
		dotPro = dotPro + ipVector1.values[i] * ipVector2.values[i];
	}

	return dotPro;
}

template<class TT>
inline std::ostream & operator<<(std::ostream & output, const VectorND<TT>& ipVector)
{
	for (int i = ipVector.iStart; i <= ipVector.iEnd; i++)
	{
		output << ipVector.index(i) << " " << ipVector(i) << endl;;
	}
	return output;
}

template<class TT>
VectorND<TT> normalize(const VectorND<TT>& ipVector)
{
	VectorND<TT> returnVector = ipVector;
	returnVector.normalize();
	return returnVector;
}
