#pragma once

#include "CombineStructure.h"


template<class TT>
class Field2D
{
public:
	Grid2D grid;
	Grid2D ghostGrid;
	Array2D<TT> dataArray;
	Array2D<TT> ghostDataArray;

	int ghostWidth; // Default value is 1.
	int iRes, jRes;
	int iStart, jStart, iEnd, jEnd;
	double xMin, yMin, xMax, yMax;
	double xLength, yLength;
	double dx, dy;
	double twodx, twody;
	double dx2, dy2;
	double oneOverdx, oneOverdy;
	double oneOver2dx, oneOver2dy;
	double oneOverdx2, oneOverdy2;

	Field2D();
	~Field2D();


	// Non-Ghost Grid
	Field2D(const Grid2D& ipGrid);
	Field2D(const double& ipXMin, const double& ipXmax, const int& ipiRes, const double& ipYMin, const double& ipYmax, const int& ipjRes);
	Field2D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);

	void initialize(const Grid2D& ipGrid);
	void initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);


	// Ghost Grid
	Field2D(const int& ipGhostWidth, const Grid2D& ipGrid);
	Field2D(const int& ipGhostWidth, const double& ipXMin, const double& ipXmax, const int& ipiRes, const double& ipYMin, const double& ipYmax, const int& ipjRes);
	Field2D(const int& ipGhostWidth, const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);

	void initialize(const int& ipGhostWidth, const Grid2D& ipGrid);
	//void initialize(const int& ipGhostWidth, double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);




	inline TT& operator [](const int& i) const;

	const int index(const Vector2D<int>& ipVector) const;

	const int index(const int& i, const int& j) const;

	inline TT& operator ()(const int& i) const;

	inline TT& operator ()(const Vector2D<int>& ipVector) const;

	inline TT& operator ()(const int& i, const int& j) const;

	inline TT operator ()(const double& x, const double& y) const;

	inline void operator = (const Field2D<TT>& ipField);

	inline void operator *=(const double& constant);

	inline void operator +=(const double& constant);

	inline void operator -=(const double& constant);

	inline void operator /=(const double& constant);

	Field2D<TT> operator + (const Field2D<TT>& ipField) const;

	Field2D<TT> operator - (const Field2D<TT>& ipField) const;

	Field2D<TT> operator * (const Field2D<TT>& ipField) const;

	Field2D<TT> operator / (const Field2D<TT>& ipField) const;

	Field2D<TT> operator + (const TT& constant) const;

	Field2D<TT> operator - (const TT& constant) const;

	Field2D<TT> operator * (const TT& constant) const;

	Field2D<TT> operator / (const TT& constant) const;

	//inline TT& operator ()(const Vector2D <double>& ipVector) const
	//{
	//	assert(ipVector.x >= xMin && ipVector.x <= xMax);
	//	assert(ipVector.y >= yMin && ipVector.y <= yMax);

	//	return interpolation(ipVector.x, ipVector.y);
	//}

	
	// Write MATLAB Variable
	inline void Variable(const char * varName);

	inline void WriteFile(const string& fileName);

	inline TT interpolation(const double& x, const double& y) const;
	inline TT interpolation(const Vector2D<double>& ipVector) const;
	
	inline void FillGhostCell();

	inline Vector2D<int> containedCell(const double& x, const double& y) const;

	inline Vector2D<double> Gradient(const int& i, const int& j);
	static Field2D<Vector2D<double>> Gradient(const Field2D<double>& ipField);

	inline TT minmod(const TT& constant1, const TT constant2) const;

	// Derivative
	inline TT dxxPhi(const int& i, const int& j) const;
	inline TT dxPhi(const int& i, const int& j) const;
	inline TT dxPlusPhi(const int& i, const int& j) const;
	inline TT dxMinusPhi(const int& i, const int& j) const;

	inline TT dyyPhi(const int& i, const int& j) const;
	inline TT dyPhi(const int& i, const int& j) const;
	inline TT dyPlusPhi(const int& i, const int& j) const;
	inline TT dyMinusPhi(const int& i, const int& j) const;

	inline TT dxyPhi(const int& i, const int& j) const;

	inline TT dxPlusPhiSubcell(const int& i, const int& j) const;
	inline TT dxMinusPhiSubcell(const int& i, const int& j) const;
	inline TT dyPlusPhiSubcell(const int& i, const int& j) const;
	inline TT dyMinusPhiSubcell(const int& i, const int& j) const;
	
	static TT L1Norm(const Field2D<TT>& ipField1);
	static TT L2Norm(const Field2D<TT>& ipField1);
	static TT InnerProduct(const Field2D<TT>& ipField1, const Field2D<TT>& ipField2);

	inline TT Divegence(const int& i, const int& j) const;
	static Field2D<double> Divegence(const Field2D<Vector2D<double>>& ipField);
	inline Array2D<double> Hessian(const int& i, const int& j) const;
	//static Array2D<double> Hessian(const Field2D<double>& ipField);
private:

};



//#endif // !Field2D_H



template<class TT>
Field2D<TT>::Field2D()
{
}

template<class TT>
Field2D<TT>::~Field2D()
{
}

template<class TT>
inline Field2D<TT>::Field2D(const Grid2D & ipGrid)
{
	grid = ipGrid;
	initialize(1, grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const double & ipXMin, const double & ipXmax, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjRes)
{
	grid.initialize(ipXMin, ipXmax, 0, ipiRes, ipYMin, ipYmax, 0, ipjRes);
	initialize(1, grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	grid.initialize(ipXMin, ipXmax, ipiStart, ipiRes, ipYMin, ipYmax, ipjStart, ipjRes);
	initialize(1, grid);
}

template<class TT>
inline void Field2D<TT>::initialize(const Grid2D & ipGrid)
{
	initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes, ipGrid.yMin, ipGrid.yMax, ipGrid.jStart, ipGrid.jRes);
	dataArray = Array2D<TT>(grid);
}

template<class TT>
inline void Field2D<TT>::initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	iRes = ipiRes;
	jRes = ipjRes;
	iStart = ipiStart;
	jStart = ipjStart;
	iEnd = iStart + iRes - 1;
	jEnd = jStart + jRes - 1;
	xMin = ipXMin;
	yMin = ipYMin;
	xMax = ipXmax;
	yMax = ipYmax;
	xLength = xMax - xMin;
	yLength = yMax - yMin;
	dx = xLength / (double)(iRes - 1);
	dy = yLength / (double)(jRes - 1);
	twodx = 2.0 * dx;
	twody = 2.0*dy;
	dx2 = dx*dx;
	dy2 = dy*dy;
	oneOverdx = 1.0 / dx;
	oneOverdy = 1.0 / dy;
	oneOver2dx = 1.0 / twodx;
	oneOver2dy = 1.0 / twody;
	oneOverdx2 = oneOverdx*oneOverdx;
	oneOverdy2 = oneOverdy*oneOverdy;
}

template<class TT>
inline Field2D<TT>::Field2D(const int & ipGhostWidth, const Grid2D & ipGrid)
{
	grid = ipGrid;

	initialize(ipGhostWidth, grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const int & ipGhostWidth, const double & ipXMin, const double & ipXmax, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjRes)
{
	grid.initialize(ipXMin, ipXmax, 0, ipiRes, ipYMin, ipYmax, 0, ipjRes);

	initialize(ipGhostWidth, grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const int & ipGhostWidth, const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	grid.initialize(ipXMin, ipXmax, ipiStart, ipiRes, ipYMin, ipYmax, ipjStart, ipjRes);

	initialize(ipGhostWidth, grid);
}

template<class TT>
inline void Field2D<TT>::initialize(const int & ipGhostWidth, const Grid2D & ipGrid)
{
	ghostWidth = ipGhostWidth;

	initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes, ipGrid.yMin, ipGrid.yMax, ipGrid.jStart, ipGrid.jRes);
	dataArray = Array2D<TT>(ipGrid);

	double widthX = double(ghostWidth)*(ipGrid.xMax - ipGrid.xMin) / double(ipGrid.iRes - 1);
	double widthY = double(ghostWidth)*(ipGrid.yMax - ipGrid.yMin) / double(ipGrid.jRes - 1);

	ghostGrid.initialize(ipGrid.xMin - widthX, ipGrid.xMax + widthX, ipGrid.iStart - ghostWidth, ipGrid.iRes + 2 * ghostWidth, ipGrid.yMin - widthY, ipGrid.yMax + widthY, ipGrid.jStart - ghostWidth, ipGrid.jRes + 2 * ghostWidth);
	ghostDataArray = Array2D<TT>(ghostGrid);
}



template<class TT>
inline TT & Field2D<TT>::operator[](const int & i) const
{
	assert(i >= 0 && i < iRes*jRes);
	return dataArray(i);
}

template<class TT>
const int Field2D<TT>::index(const Vector2D<int>& ipVector) const
{
	assert(ipVector[0] >= iStart && ipVector[0] <= iEnd);
	assert(ipVector[1] >= jStart && ipVector[1] <= jEnd);
	return dataArray.index(ipVector);
}

template<class TT>
const int Field2D<TT>::index(const int & i, const int & j) const
{
	if ((i >= iStart && i <= iEnd)&& (j >= jStart && j <= jEnd))
	{
		return dataArray.index(i, j);

	}
	else if ((i >= iStart - ghostWidth && i <= iEnd + ghostWidth) && (j >= jStart - ghostWidth && j <= jEnd + ghostWidth))
	{
		return ghostDataArray.index(i, j);
	}
	else
	{
		assert(i >= iStart);
		assert(i <= iEnd);
		assert(j >= jStart);
		assert(j <= jEnd);
		assert(i >= iStart - ghostWidth);
		assert(i <= iEnd + ghostWidth);
		assert(j >= jStart - ghostWidth);
		assert(j <= jEnd + ghostWidth);
	}
}

template<class TT>
inline TT & Field2D<TT>::operator()(const int & i) const
{
	assert(i >= iStart && i <= iEnd);
	return dataArray(i);
}

template<class TT>
inline TT & Field2D<TT>::operator()(const Vector2D<int>& ipVector) const
{
	assert(ipVector[0] >= iStart && ipVector[0] <= iEnd);
	assert(ipVector[1] >= jStart && ipVector[1] <= jEnd);
	return dataArray(ipVector[0], ipVector[1]);
}

template<class TT>
inline TT & Field2D<TT>::operator()(const int & i, const int & j) const
{
	if ((i >= iStart && i <= iEnd) && (j >= jStart && j <= jEnd))
	{
		return dataArray(i, j);
	}
	else if ((i >= iStart - ghostWidth && i <= iEnd + ghostWidth) && (j >= jStart - ghostWidth && j <= jEnd + ghostWidth))
	{
		return ghostDataArray(i, j);
	}
	else
	{
		cout << "(i,j)=(" << i << "," << j << ")" << endl;
		cout << "iStart = " << iStart<< endl;
		cout << "iEnd   = " << iEnd << endl;
		cout << "jStart = " << jStart << endl;
		cout << "jEnd   = " << jEnd << endl;
		cout << "jEnd   = " << jEnd << endl;
		assert(i >= iStart);
		assert(i <= iEnd);
		assert(j >= jStart);
		assert(j <= jEnd);
		assert(i >= iStart - ghostWidth);
		assert(i <= iEnd + ghostWidth);
		assert(j >= jStart - ghostWidth);
		assert(j <= jEnd + ghostWidth);
	}
}

template<class TT>
inline TT Field2D<TT>::operator()(const double & x, const double & y) const
{
	assert(x >= xMin && x <= xMax);
	assert(y >= yMin && y <= yMax);
	//TT& a=x;// = interpolation(x, y);
	Vector2D<double> xy(x, y);
	Vector2D<int> cell = containedCell(x, y);
	TT value00, value10, value01, value11;
	value00 = dataArray(cell);
	value10 = dataArray(cell.i + 1, cell.j);
	value01 = dataArray(cell.i, cell.j + 1);
	value11 = dataArray(cell.i + 1, cell.j + 1);

	double distance00, distance10, distance01, distance11;
	distance00 = (grid(cell) - xy).magnitude();
	distance10 = (grid(cell.i + 1, cell.j) - xy).magnitude();
	distance01 = (grid(cell.i, cell.j + 1) - xy).magnitude();
	distance11 = (grid(cell.i + 1, cell.j + 1) - xy).magnitude();

	return ((dataArray(cell)*distance00 + dataArray(cell.i + 1, cell.j)*distance10 + dataArray(cell.i, cell.j + 1)*distance01 + dataArray(cell.i + 1, cell.j + 1)*distance11) / (distance00 + distance01 + distance10 + distance11));
	//return dataArray(1,1);
}

template<class TT>
inline void Field2D<TT>::operator=(const Field2D<TT>& ipField)
{
	grid = ipField.grid;
	ghostGrid = ipField.ghostGrid;
	dataArray = ipField.dataArray;
	ghostDataArray = ipField.ghostDataArray;

	ghostWidth = ipField.ghostWidth;
	iRes = ipField.iRes;
	jRes = ipField.jRes;
	iStart = ipField.iStart;
	jStart = ipField.jStart;
	iEnd = iStart + iRes - 1;
	jEnd = jStart + jRes - 1;
	xMin = ipField.xMin;
	yMin = ipField.yMin;
	xMax = ipField.xMax;
	yMax = ipField.yMax;
	xLength = xMax - xMin;
	yLength = yMax - yMin;
	dx = xLength / (double)(iRes - 1);
	dy = yLength / (double)(jRes - 1);
	twodx = 2.0 * dx;
	twody = 2.0*dy;
	dx2 = dx*dx;
	dy2 = dy*dy;
	oneOverdx = 1.0 / dx;
	oneOverdy = 1.0 / dy;
	oneOver2dx = 1.0 / twodx;
	oneOver2dy = 1.0 / twody;
	oneOverdx2 = oneOverdx*oneOverdx;
	oneOverdy2 = oneOverdy*oneOverdy;
}



template<class TT>
inline void Field2D<TT>::operator*=(const double & constant)
{
#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			dataArray(i,j) *= constant;
		}
	}
}

template<class TT>
inline void Field2D<TT>::operator+=(const double & constant)
{
#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			dataArray(i, j) += constant;
		}
	}
}

template<class TT>
inline void Field2D<TT>::operator-=(const double & constant)
{
#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			dataArray(i, j) -= constant;
		}
	}
}

template<class TT>
inline void Field2D<TT>::operator/=(const double & constant)
{
	assert(constant != 0);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			dataArray(i, j) *= 1 / constant;
		}
	}
}

template<class TT>
Field2D<TT> Field2D<TT>::operator+(const Field2D<TT>& ipField) const
{
	Field2D<TT> tempField(ipField.ghostWidth,ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			tempField(i, j) = dataArray(i, j) + ipField(i, j);

		}
	}
	return tempField;
}

template<class TT>
Field2D<TT> Field2D<TT>::operator-(const Field2D<TT>& ipField) const
{
	FField2D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			tempField(i, j) = dataArray(i, j) - ipField(i, j);
		}
	}
	return tempField;
}

template<class TT>
Field2D<TT> Field2D<TT>::operator*(const Field2D<TT>& ipField) const
{
	Field2D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			tempField(i, j) = dataArray(i, j) * ipField(i, j);
		}
	}
	return tempField;
}

template<class TT>
Field2D<TT> Field2D<TT>::operator/(const Field2D<TT>& ipField) const
{
	Field2D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			if (ipField(i, j) != 0)
			{
				tempField(i, j) = dataArray(i, j) / ipField(i, j);
			}
		}
	}
	return tempField;
}

template<class TT>
Field2D<TT> Field2D<TT>::operator+(const TT & constant) const
{
	Field2D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			tempField(i, j) = constant + dataArray(i, j);
		}
	}
	return tempField;
}

template<class TT>
Field2D<TT> Field2D<TT>::operator-(const TT & constant) const
{
	Field2D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			tempField(i, j) = dataArray(i, j) - constant;
		}
	}
	return tempField;
}

template<class TT>
Field2D<TT> Field2D<TT>::operator*(const TT & constant) const
{
	Field2D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			tempField(i, j) = dataArray(i, j) - constant;
		}
	}
	return tempField;
}

template<class TT>
Field2D<TT> Field2D<TT>::operator/(const TT & constant) const
{
	assert(constant != 0); 
	Field2D<TT> tempField(ipField.ghostWidth, ipField.grid);

#pragma omp parallel for
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			tempField(i, j) = dataArray(i, j) / constant;
		}
	}
	return tempField;
}

template<class TT>
inline void Field2D<TT>::Variable(const char * varName)
{
	MATLAB.Variable(varName, iRes, jRes, dataArray.values);
}


template<class TT>
inline void Field2D<TT>::WriteFile(const string & fileName)
{
	ofstream solutionFile;
	solutionFile.open("D:\\Data/" + fileName + ".txt", ios::binary);
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			solutionFile << setprecision(16) << i << " " << j << " " << grid(i, j) << " " << dataArray(i, j) << endl;
		}
	}
	solutionFile.close();
}

template<class TT>
inline TT Field2D<TT>::interpolation(const double & x, const double & y) const
{
	assert(x >= xMin && x <= xMax);
	assert(y >= yMin && y <= yMax);

	Vector2D<double> xy(x, y);
	Vector2D<int> cell = containedCell(x, y);

	double distance00, distance10, distance01, distance11;
	distance00 = (grid(cell) - xy).magnitude();
	distance10 = (grid(cell.i + 1, cell.j) - xy).magnitude();
	distance01 = (grid(cell.i, cell.j + 1) - xy).magnitude();
	distance11 = (grid(cell.i + 1, cell.j + 1) - xy).magnitude();

	return ((dataArray(cell)*distance00 + dataArray(cell.i + 1, cell.j)*distance10 + dataArray(cell.i, cell.j + 1)*distance01 + dataArray(cell.i + 1, cell.j + 1)*distance11) / (distance00 + distance01 + distance10 + distance11));
}

template<class TT>
inline TT Field2D<TT>::interpolation(const Vector2D<double>& ipVector) const
{
	return interpolation(ipVector.x, ipVector.y);
}

template<class TT>
inline void Field2D<TT>::FillGhostCell()
{
	if (ghostWidth<=0)
	{
		return;
	}


	for (int j = jStart; j <=jEnd; j++)
	{
		for (int i = iStart; i <= iStart + 1; i++)
		{
			ghostDataArray(i, j) = dataArray(i, j);
		}
		for (int i = iEnd - 1; i <= iEnd; i++)
		{
			ghostDataArray(i, j) = dataArray(i, j);
		}

		// Left
		for (int i = iStart-1; i <= iStart-ghostWidth; i--)
		{
			ghostDataArray(i, j) = ghostDataArray(i + 1, j) + oneOverdx*(ghostDataArray(i + 1, j) - ghostDataArray(i + 2, j));
		}

		// Right
		for (int i = iEnd + 1; i <= iEnd + ghostWidth; i++)
		{
			ghostDataArray(i, j) = ghostDataArray(i - 1, j) + oneOverdx*(ghostDataArray(i - 1, j) - ghostDataArray(i - 2, j));
		}
	}

	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jStart + 1; j++)
		{
			ghostDataArray(i, j) = dataArray(i, j);
		}
		for (int j = jEnd - 1; j <= jEnd; j++)
		{
			ghostDataArray(i, j) = dataArray(i, j);
		}

		// Bottom
		for (int j = jStart - 1; j <= jStart - ghostWidth; j--)
		{
			ghostDataArray(i, j) = ghostDataArray(i, j + 1) + oneOverdy*(ghostDataArray(i, j + 1) - ghostDataArray(i, j + 2));
		}

		// Top
		for (int j = jEnd + 1; j <= jEnd + ghostWidth; j++)
		{
			ghostDataArray(i, j) = ghostDataArray(i, j - 1) + oneOverdy*(ghostDataArray(i, j - 1) - ghostDataArray(i, j - 2));
		}
	}

	// Corner
	ghostDataArray(iStart - 1, jStart - 1) = ghostDataArray(iStart, jStart) + oneOverdx*(ghostDataArray(iStart, jStart) - ghostDataArray(iStart + 1, jStart)) + oneOverdy*(ghostDataArray(iStart, jStart) - ghostDataArray(iStart, jStart + 1));
	ghostDataArray(iStart - 1, jEnd + 1)   = ghostDataArray(iStart, jEnd)   + oneOverdx*(ghostDataArray(iStart, jEnd) - ghostDataArray(iStart + 1, jEnd))     + oneOverdy*(ghostDataArray(iStart, jEnd) - ghostDataArray(iStart, jEnd - 1));
	ghostDataArray(iEnd + 1, jStart - 1)   = ghostDataArray(iEnd, jStart)   + oneOverdx*(ghostDataArray(iEnd, jStart) - ghostDataArray(iEnd - 1, jStart))     + oneOverdy*(ghostDataArray(iEnd, jStart) - ghostDataArray(iEnd, jStart + 1));
	ghostDataArray(iEnd + 1, jEnd + 1)     = ghostDataArray(iEnd, jEnd)     + oneOverdx*(ghostDataArray(iEnd, jEnd) - ghostDataArray(iEnd - 1, jEnd))         + oneOverdy*(ghostDataArray(iEnd, jEnd) - ghostDataArray(iEnd, jEnd - 1));
}

template<class TT>
inline Vector2D<int> Field2D<TT>::containedCell(const double & x, const double & y) const
{
	return Vector2D<int>(floor((x - xMin)*oneOverdx), floor((y - yMin)*oneOverdy));
}

template<class TT>
inline Vector2D<double> Field2D<TT>::Gradient(const int & i, const int & j)
{
	return Vector2D<double>(dxPhi(i, j), dyPhi(i, j));
}

template<class TT>
inline Field2D<Vector2D<double>> Field2D<TT>::Gradient(const Field2D<double>& ipField)
{
	Field2D<Vector2D<double>> gradient(ipField.grid);
#pragma omp parallel for
	for (int i = gradient.iStart; i <= gradient.iEnd; i++)
	{
		for (int j = gradient.jStart; j <= gradient.jEnd; j++)
		{
			gradient(i, j) = Vector2D<double>(ipField.dxPhi(i, j), ipField.dyPhi(i, j));
		}
	}
	return gradient;
}

template<class TT>
inline TT Field2D<TT>::minmod(const TT & constant1, const TT constant2) const
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
inline TT Field2D<TT>::dxxPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (dataArray(i + 1, j) - 2 * dataArray(i, j) + dataArray(i - 1, j))*grid.oneOverdx2;
	}
	else if (i == grid.iStart)
	{
		return (dataArray(i, j) - 2 * dataArray(i + 1, j) + dataArray(i + 2, j))*grid.oneOverdx2;
	}
	else
	{
		return (dataArray(i - 2, j) - 2 * dataArray(i - 1, j) + dataArray(i, j))*grid.oneOverdx2;
	}
}

template<class TT>
inline TT Field2D<TT>::dxPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (dataArray(i + 1, j) - dataArray(i - 1, j))*grid.oneOver2dx;// -grid.dx / 2.0*minmod(dxxPhi(i, j), dxxPhi(i + 1, j));
	}
	else if (i == grid.iStart)
	{
		return dxPlusPhi(i, j);
	}
	else
	{
		return dxMinusPhi(i, j);
	}
}

template<class TT>
inline TT Field2D<TT>::dxPlusPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i < grid.iEnd)
	{
		return (dataArray(i + 1, j) - dataArray(i, j))*grid.oneOverdx;
	}
	else
	{
		return (dataArray(i, j) - dataArray(i - 1, j))*grid.oneOverdx;
	}
}

template<class TT>
inline TT Field2D<TT>::dxMinusPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart)
	{
		return (dataArray(i, j) - dataArray(i - 1, j))*grid.oneOverdx;
	}
	else
	{
		return (dataArray(i + 1, j) - dataArray(i, j))*grid.oneOverdx;
	}
}

template<class TT>
inline TT Field2D<TT>::dyyPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (dataArray(i, j + 1) - 2 * dataArray(i, j) + dataArray(i, j - 1))*grid.oneOverdy2;
	}
	else if (j == grid.jStart)
	{
		return (dataArray(i, j) - 2 * dataArray(i, j + 1) + dataArray(i, j + 2))*grid.oneOverdy2;
	}
	else
	{
		return (dataArray(i, j - 2) - 2 * dataArray(i, j - 1) + dataArray(i, j))*grid.oneOverdy2;
	}
}

template<class TT>
inline TT Field2D<TT>::dyPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (dataArray(i, j + 1) - dataArray(i, j - 1))*grid.oneOver2dy;// -grid.dy / 2.0*minmod(dyyPhi(i, j), dyyPhi(i, j + 1));
	}
	else if (j == grid.jStart)
	{
		return dyPlusPhi(i, j);
	}
	else
	{
		return dyMinusPhi(i, j);
	}
}

template<class TT>
inline TT Field2D<TT>::dyPlusPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j < grid.jEnd)
	{
		return (dataArray(i, j + 1) - dataArray(i, j))*grid.oneOverdy;
	}
	else
	{
		return (dataArray(i, j) - dataArray(i, j - 1))*grid.oneOverdy;
	}
}

template<class TT>
inline TT Field2D<TT>::dyMinusPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart)
	{
		return (dataArray(i, j) - dataArray(i, j - 1))*grid.oneOverdy;
	}
	else
	{
		return (dataArray(i, j + 1) - dataArray(i, j))*grid.oneOverdy;
	}
}

template<class TT>
inline TT Field2D<TT>::dxyPhi(const int & i, const int & j) const
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (dxPhi(i, j + 1) - dxPhi(i, j - 1))*grid.oneOver2dy;
	}
	else if (j == grid.iStart)
	{
		return (dxPhi(i, j + 1) - dxPhi(i, j))*grid.oneOverdy;
	}
	else
	{
		return (dxPhi(i, j) - dxPhi(i, j - 1))*grid.oneOverdy;
	}
}

template<class TT>
inline TT Field2D<TT>::dxPlusPhiSubcell(const int & i, const int & j) const
{
	return TT();
}

template<class TT>
inline TT Field2D<TT>::dxMinusPhiSubcell(const int & i, const int & j) const
{
	return TT();
}

template<class TT>
inline TT Field2D<TT>::dyPlusPhiSubcell(const int & i, const int & j) const
{
	return TT();
}

template<class TT>
inline TT Field2D<TT>::dyMinusPhiSubcell(const int & i, const int & j) const
{
	return TT();
}

template<class TT>
inline TT Field2D<TT>::L1Norm(const Field2D<TT>& ipField1)
{

	return Array2D<TT>::L1Norm(ipField1.dataArray);
}

template<class TT>
inline TT Field2D<TT>::L2Norm(const Field2D<TT>& ipField1)
{
	return Array2D<TT>::L2Norm(ipField1.dataArray);
}

template<class TT>
inline TT Field2D<TT>::InnerProduct(const Field2D<TT>& ipField1, const Field2D<TT>& ipField2)
{
	return Array2D<TT>::InnerProduct(ipField1.dataArray, ipField2.dataArray);
}



template<class TT>
inline TT Field2D<TT>::Divegence(const int & i, const int & j) const
{
	return dxPhi(i, j) + dyPhi(i, j);
}

template<class TT>
inline Field2D<double> Field2D<TT>::Divegence(const Field2D<Vector2D<double>>& ipField)
{
	Field2D<double> divergence(ipField.grid);
	//FV dxField(ipField.grid);
	//FV dyField(ipField.grid);

#pragma omp parallel for
	for (int i = divergence.iStart; i <= divergence.iEnd; i++)
	{
		for (int j = divergence.jStart; j <= divergence.jEnd; j++)
		{
			divergence(i, j) = 0;// ipField.dxPhi(i, j);// +ipField.dyPhi(i, j);
		}
	}
	return divergence;
}

template<class TT>
inline Array2D<double> Field2D<TT>::Hessian(const int & i, const int & j) const
{
	Array2D<double> hessian(2, 2);
	hessian(0, 0) = dxxPhi(i, j);
	hessian(0, 1) = dxyPhi(i, j);
	hessian(1, 0) = dxyPhi(i, j);
	hessian(1, 1) = dyyPhi(i, j);
	return hessian;
}




// Typedef
//typedef double			T;
typedef Vector2D<double> VT;
typedef Vector2D<int> VI;
typedef VectorND<double> VTN;
typedef VectorND<int> VIN;
typedef Array2D<double> AD;
typedef Array2D<int> AI;
typedef Field2D<double> FD;
typedef Field2D<Vector2D<double>> FV;

//typedef complex<T>		Tcomp;
//typedef complex<int>    Icomp;
