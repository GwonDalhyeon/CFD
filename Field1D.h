#pragma once


#include "CommonDef.h"
#include "Vector2D.h"
#include "VectorND.h"
#include "Grid1D.h"
#include "Array1D.h"
#include "CombineStructure.h"

template<class TT>
class Field1D
{
public:
	Grid1D grid;
	Grid1D ghostGrid;
	Array1D<TT> dataArray;
	Array1D<TT> ghostDataArray;

	int ghostWidth; // Default value is 1.
	int iRes;
	int iStart, iEnd;
	double xMin, xMax;
	double xLength;
	double dx;
	double twodx;
	double dx2;
	double oneOverdx;
	double oneOver2dx;
	double oneOverdx2;

	Field1D();
	~Field1D();


	// Non-Ghost Grid
	Field1D(const Grid1D& ipGrid);
	Field1D(const double& ipXMin, const double& ipXmax, const int& ipiRes);
	Field1D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes);

	void initialize(const Grid1D& ipGrid);
	void initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes);


	// Ghost Grid
	Field1D(const int& ipGhostWidth, const Grid1D& ipGrid);
	Field1D(const int& ipGhostWidth, const double& ipXMin, const double& ipXmax, const int& ipiRes);
	Field1D(const int& ipGhostWidth, const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes);

	void initialize(const int& ipGhostWidth, const Grid1D& ipGrid);




	inline TT& operator [](const int& i) const;

	const int index(const int& i) const;

	inline TT& operator ()(const int& i) const;

	inline TT operator ()(const double& x) const;

	inline void operator = (const Field1D<TT>& ipField);

	inline void operator *=(const double& constant);

	inline void operator +=(const double& constant);

	inline void operator -=(const double& constant);

	inline void operator /=(const double& constant);

	Field1D<TT> operator + (const Field1D<TT>& ipField) const;

	Field1D<TT> operator - (const Field1D<TT>& ipField) const;

	Field1D<TT> operator * (const Field1D<TT>& ipField) const;

	Field1D<TT> operator / (const Field1D<TT>& ipField) const;

	Field1D<TT> operator + (const TT& constant) const;

	Field1D<TT> operator - (const TT& constant) const;

	Field1D<TT> operator * (const TT& constant) const;

	Field1D<TT> operator / (const TT& constant) const;

	//inline TT& operator ()(const Vector2D <double>& ipVector) const
	//{
	//	assert(ipVector.x >= xMin && ipVector.x <= xMax);
	//	assert(ipVector.y >= yMin && ipVector.y <= yMax);

	//	return interpolation(ipVector.x, ipVector.y);
	//}


	// Write MATLAB Variable
	inline void Variable(const char * varName);

	inline void WriteFile(const string& fileName);

	inline TT interpolation(const double& x) const;

	inline void FillGhostCell();

	inline VI containedCell(const double& x) const;

	inline VT gradient(const int& i);

	inline TT minmod(const TT& constant1, const TT constant2) const;

	// Derivative
	inline TT dxxPhi(const int& i) const;
	inline TT dxPhi(const int& i) const;
	inline TT dxPlusPhi(const int& i) const;
	inline TT dxMinusPhi(const int& i) const;

	inline TT dxPlusPhiSubcell(const int& i) const;
	inline TT dxMinusPhiSubcell(const int& i) const;

	static TT L1Norm(const Field1D<TT>& ipField1);
	static TT L2Norm(const Field1D<TT>& ipField1);
	static TT InnerProduct(const Field1D<TT>& ipField1, const Field1D<TT>& ipField2);

	static Field1D<double> Gradient(const Field1D<double>& ipField);
private:

};



//#endif // !Field1D_H



template<class TT>
Field1D<TT>::Field1D()
{
}

template<class TT>
Field1D<TT>::~Field1D()
{
}

template<class TT>
inline Field1D<TT>::Field1D(const Grid1D & ipGrid)
{
	grid = ipGrid;
	initialize(1, grid);
}

template<class TT>
inline Field1D<TT>::Field1D(const double & ipXMin, const double & ipXmax, const int & ipiRes)
{
	grid.Initialize(ipXMin, ipXmax, 0, ipiRes);
	initialize(1, grid);
}

template<class TT>
inline Field1D<TT>::Field1D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes)
{
	grid.Initialize(ipXMin, ipXmax, ipiStart, ipiRes);
	initialize(1, grid);
}

template<class TT>
inline void Field1D<TT>::initialize(const Grid1D & ipGrid)
{
	initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes);
	dataArray = Array1D<TT>(grid);
}

template<class TT>
inline void Field1D<TT>::initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes)
{
	iRes = ipiRes;
	iStart = ipiStart;
	iEnd = iStart + iRes - 1;
	xMin = ipXMin;
	xMax = ipXmax;
	xLength = xMax - xMin;
	dx = xLength / (double)(iRes - 1);
	twodx = 2.0 * dx;
	dx2 = dx*dx;
	oneOverdx = 1.0 / dx;
	oneOver2dx = 1.0 / twodx;
	oneOverdx2 = oneOverdx*oneOverdx;
}

template<class TT>
inline Field1D<TT>::Field1D(const int & ipGhostWidth, const Grid1D & ipGrid)
{
	grid = ipGrid;

	initialize(ipGhostWidth, grid);
}

template<class TT>
inline Field1D<TT>::Field1D(const int & ipGhostWidth, const double & ipXMin, const double & ipXmax, const int & ipiRes)
{
	grid.Initialize(ipXMin, ipXmax, 0, ipiRes);

	initialize(ipGhostWidth, grid);
}

template<class TT>
inline Field1D<TT>::Field1D(const int & ipGhostWidth, const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes)
{
	grid.Initialize(ipXMin, ipXmax, ipiStart, ipiRes);

	initialize(ipGhostWidth, grid);
}

template<class TT>
inline void Field1D<TT>::initialize(const int & ipGhostWidth, const Grid1D & ipGrid)
{
	ghostWidth = ipGhostWidth;

	initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes);
	dataArray = Array1D<TT>(ipGrid);

	double widthX = double(ghostWidth)*(ipGrid.xMax - ipGrid.xMin) / double(ipGrid.iRes - 1);

	ghostGrid.Initialize(ipGrid.xMin - widthX, ipGrid.xMax + widthX, ipGrid.iStart - ghostWidth, ipGrid.iRes + 2 * ghostWidth);
	ghostDataArray = Array1D<TT>(ghostGrid);
}



template<class TT>
inline TT & Field1D<TT>::operator[](const int & i) const
{
	assert(i >= 0 && i < iRes);
	return dataArray(i);
}

template<class TT>
const int Field1D<TT>::index(const int & i) const
{
	if (i >= iStart && i <= iEnd)
	{
		return dataArray.index(i);

	}
	else if (i >= iStart - ghostWidth && i <= iEnd + ghostWidth)
	{
		return ghostDataArray.index(i);
	}
	else
	{
		assert(i >= iStart);
		assert(i <= iEnd);
		assert(i >= iStart - ghostWidth);
		assert(i <= iEnd + ghostWidth);
	}
}

template<class TT>
inline TT & Field1D<TT>::operator()(const int & i) const
{
	if (i >= iStart && i <= iEnd)
	{
		return dataArray(i);
	}
	else if (i >= iStart - ghostWidth && i <= iEnd + ghostWidth)
	{
		return ghostDataArray(i);
	}
	else
	{
		cout << "i=" << i << endl;
		cout << "iStart = " << iStart << endl;
		cout << "iEnd   = " << iEnd << endl;
		assert(i >= iStart);
		assert(i <= iEnd);
		assert(i >= iStart - ghostWidth);
		assert(i <= iEnd + ghostWidth);
	}
}

template<class TT>
inline TT Field1D<TT>::operator()(const double & x) const
{
	assert(x >= xMin && x <= xMax);
	//TT& a=x;// = interpolation(x, y);
	int cell = containedCell(x);
	TT value0, value1;
	value0 = dataArray(cell);
	value1 = dataArray(cell.i + 1);

	double distance0, distance1;
	distance00 = abs(grid(cell) - x);
	distance10 = abs(grid(cell.i + 1) - x);

	return ((dataArray(cell)*distance00 + dataArray(cell.i + 1)*distance1) / (distance0 + distance1));
	//return dataArray(1,1);
}

template<class TT>
inline void Field1D<TT>::operator=(const Field1D<TT>& ipField)
{
	grid = ipField.grid;
	ghostGrid = ipField.ghostGrid;
	dataArray = ipField.dataArray;
	ghostDataArray = ipField.ghostDataArray;

	ghostWidth = ipField.ghostWidth;
	iRes = ipField.iRes;
	iStart = ipField.iStart;
	iEnd = iStart + iRes - 1;
	xMin = ipField.xMin;
	xMax = ipField.xMax;
	xLength = xMax - xMin;
	dx = xLength / (double)(iRes - 1);
	twodx = 2.0 * dx;
	dx2 = dx*dx;
	oneOverdx = 1.0 / dx;
	oneOver2dx = 1.0 / twodx;
	oneOverdx2 = oneOverdx*oneOverdx;
}



template<class TT>
inline void Field1D<TT>::operator*=(const double & constant)
{
#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		dataArray(i) *= constant;
	}
}

template<class TT>
inline void Field1D<TT>::operator+=(const double & constant)
{
#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		dataArray(i) += constant;
	}
}

template<class TT>
inline void Field1D<TT>::operator-=(const double & constant)
{
#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		dataArray(i) -= constant;
	}
}

template<class TT>
inline void Field1D<TT>::operator/=(const double & constant)
{
	assert(constant != 0);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
			dataArray(i) *= 1 / constant;
	}
}

template<class TT>
Field1D<TT> Field1D<TT>::operator+(const Field1D<TT>& ipField) const
{
	Field1D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
			tempField(i) = dataArray(i) + ipField(i);
	}
	return tempField;
}

template<class TT>
Field1D<TT> Field1D<TT>::operator-(const Field1D<TT>& ipField) const
{
	FField1D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
			tempField(i) = dataArray(i) - ipField(i);
	}
	return tempField;
}

template<class TT>
Field1D<TT> Field1D<TT>::operator*(const Field1D<TT>& ipField) const
{
	Field1D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
			tempField(i) = dataArray(i) * ipField(i);
	}
	return tempField;
}

template<class TT>
Field1D<TT> Field1D<TT>::operator/(const Field1D<TT>& ipField) const
{
	Field1D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
			if (ipField(i) != 0)
			{
				tempField(i) = dataArray(i) / ipField(i);
			}
	}
	return tempField;
}

template<class TT>
Field1D<TT> Field1D<TT>::operator+(const TT & constant) const
{
	Field1D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
			tempField(i) = constant + dataArray(i);
	}
	return tempField;
}

template<class TT>
Field1D<TT> Field1D<TT>::operator-(const TT & constant) const
{
	Field1D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
			tempField(i) = dataArray(i) - constant;
	}
	return tempField;
}

template<class TT>
Field1D<TT> Field1D<TT>::operator*(const TT & constant) const
{
	Field1D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
			tempField(i) = dataArray(i) - constant;
	}
	return tempField;
}

template<class TT>
Field1D<TT> Field1D<TT>::operator/(const TT & constant) const
{
	assert(constant != 0);
	Field1D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
			tempField(i) = dataArray(i) / constant;
	}
	return tempField;
}

template<class TT>
inline void Field1D<TT>::Variable(const char * varName)
{
	MATLAB.Variable(varName, iRes, 1, dataArray.values);
}


template<class TT>
inline void Field1D<TT>::WriteFile(const string & fileName)
{
	ofstream solutionFile;
	solutionFile.open("D:\\Data/" + fileName + ".txt", ios::binary);
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
			solutionFile << setprecision(16) << i << " "  << grid(i) << " " << dataArray(i) << endl;
	}
	solutionFile.close();
}

template<class TT>
inline TT Field1D<TT>::interpolation(const double & x) const
{
	assert(x >= xMin && x <= xMax);

	int cell = containedCell(x);

	double distance0, distance1;
	distance0 = double(grid(cell) - x);
	distance1 = double(grid(cell.i + 1) - x);

	return ((dataArray(cell)*distance0 + dataArray(cell.i + 1)*distance1) / (distance0 + distance1));
}

//template<class TT>
//inline TT Field1D<TT>::interpolation(const double x) const
//{
//	return interpolation(x);
//}

template<class TT>
inline void Field1D<TT>::FillGhostCell()
{
	if (ghostWidth <= 0)
	{
		return;
	}

	for (int i = iStart; i <= iEnd; i++)
	{
		ghostDataArray(i) = dataArray(i);
	}

	// Left
	for (int i = iStart - 1; i <= iStart - ghostWidth; i--)
	{
		ghostDataArray(i) = ghostDataArray(i + 1) + oneOverdx*(ghostDataArray(i + 1) - ghostDataArray(i + 2));
	}

	// Right
	for (int i = iEnd + 1; i <= iEnd + ghostWidth; i++)
	{
		ghostDataArray(i) = ghostDataArray(i - 1) + oneOverdx*(ghostDataArray(i - 1) - ghostDataArray(i - 2));
	}
}

template<class TT>
inline VI Field1D<TT>::containedCell(const double & x) const
{
	return int (floor((x - xMin)*oneOverdx));
}

template<class TT>
inline TT Field1D<TT>::minmod(const TT & constant1, const TT constant2) const
{
	if (constant1*constant2<0)
	{
		return 0;
	}
	else if (abs(constant1) >= abs(constant2))
	{
		return constant2;
	}
	else
	{
		return constant1;
	}
	return TT();
}

template<class TT>
inline TT Field1D<TT>::dxxPhi(const int & i) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (dataArray(i + 1) - 2 * dataArray(i) + dataArray(i - 1))*grid.oneOverdx2;
	}
	else if (i == grid.iStart)
	{
		return (dataArray(i) - 2 * dataArray(i + 1) + dataArray(i + 2))*grid.oneOverdx2;
	}
	else
	{
		return (dataArray(i - 2) - 2 * dataArray(i - 1) + dataArray(i))*grid.oneOverdx2;
	}
}

template<class TT>
inline TT Field1D<TT>::dxPhi(const int & i) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (dataArray(i + 1) - dataArray(i - 1))*grid.oneOver2dx;// -grid.dx / 2.0*minmod(dxxPhi(i, j), dxxPhi(i + 1, j));
	}
	else if (i == grid.iStart)
	{
		return dxPlusPhi(i);
	}
	else
	{
		return dxMinusPhi(i);
	}
}

template<class TT>
inline TT Field1D<TT>::dxPlusPhi(const int & i) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);

	if (i < grid.iEnd)
	{
		return (dataArray(i + 1) - dataArray(i))*grid.oneOverdx;
	}
	else
	{
		return (dataArray(i) - dataArray(i - 1))*grid.oneOverdx;
	}
}

template<class TT>
inline TT Field1D<TT>::dxMinusPhi(const int & i) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);

	if (i > grid.iStart)
	{
		return (dataArray(i) - dataArray(i - 1))*grid.oneOverdx;
	}
	else
	{
		return (dataArray(i + 1) - dataArray(i))*grid.oneOverdx;
	}
}


template<class TT>
inline TT Field1D<TT>::dxPlusPhiSubcell(const int & i) const
{
	return TT();
}

template<class TT>
inline TT Field1D<TT>::dxMinusPhiSubcell(const int & i) const
{
	return TT();
}

template<class TT>
inline TT Field1D<TT>::L1Norm(const Field1D<TT>& ipField1)
{

	return Array1D<TT>::L1Norm(ipField1.dataArray);
}

template<class TT>
inline TT Field1D<TT>::L2Norm(const Field1D<TT>& ipField1)
{
	return Array1D<TT>::L2Norm(ipField1.dataArray);
}

template<class TT>
inline TT Field1D<TT>::InnerProduct(const Field1D<TT>& ipField1, const Field1D<TT>& ipField2)
{
	return Array1D<TT>::InnerProduct(ipField1.dataArray, ipField2.dataArray);
}

//template<class TT>
//inline Field1D<VT> Field1D<TT>::Gradient(const Field1D<double>& ipField)
//{
//	Field1D<VT> gradient(ipField.grid);
//#pragma omp parallel for
//	for (int i = gradient.iStart; i <= gradient.iEnd; i++)
//	{
//			gradient(i) = ipField.dxPhi(i)
//	}
//	return gradient;
//}