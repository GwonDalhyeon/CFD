#pragma once

#include "CombineStructure.h"


template<class TT>
class Field2D
{
public:
	Grid2D grid;
	Grid2D ghostGrid;
	Grid2D innerGrid;
	Array2D<TT> dataArray;
	Array2D<TT> ghostDataArray;
	Array2D<TT> dataArrayOld;
	Array2D<TT> ghostDataArrayOld;
	Array2D<int> BC;
	// TVDRK3 Variable
	Array2D<TT> K1;
	Array2D<TT> K2;
	Array2D<TT> K3;
	// WENO Derivation Variable
	Array2D<TT> dfdxM;
	Array2D<TT> dfdxP;
	Array2D<TT> dfdyM;
	Array2D<TT> dfdyP;
	Array2D<Vector2D<double>> gradient;

	int ghostWidth; // Default value is 1.
	
	int& iRes = grid.iRes;				int& jRes = grid.jRes;
	int& iStart = grid.iStart;			int& jStart = grid.jStart;
	int& iEnd = grid.iEnd;				int& jEnd = grid.jEnd;
	double& xMin = grid.xMin;			double& yMin = grid.yMin;
	double& xMax = grid.xMax;			double& yMax = grid.yMax;
	
	int& iResG = ghostGrid.iRes;		int& jResG = ghostGrid.jRes;
	int& iStartG = ghostGrid.iStart;	int& jStartG = ghostGrid.jStart; 
	int& iEndG = ghostGrid.iEnd;		int& jEndG = ghostGrid.jEnd;
	double& xMinG = ghostGrid.xMin;		double& yMinG = ghostGrid.yMin;
	double& xMaxG = ghostGrid.xMax;		double& yMaxG = ghostGrid.yMax;
	
	int& iResI = innerGrid.iRes;		int& jResI = innerGrid.jRes;
	int& iStartI = innerGrid.iStart;	int& jStartI = innerGrid.jStart; 
	int& iEndI = innerGrid.iEnd;		int& jEndI = innerGrid.jEnd;
	double& xMinI = innerGrid.xMin;		double& yMinI = innerGrid.yMin;
	double& xMaxI = innerGrid.xMax;		double& yMaxI = innerGrid.yMax;
	
	double& xLength = grid.xLength;			double& yLength = grid.yLength;
	double& dx = grid.dx;					double& dy = grid.dy;
	double& twodx = grid.twodx;				double& twody = grid.twody;
	double& dx2 = grid.dx2;					double& dy2 = grid.dy2;
	double& oneOverdx = grid.oneOverdx;		double& oneOverdy = grid.oneOverdy;
	double& oneOver2dx = grid.oneOver2dx;	double& oneOver2dy = grid.oneOver2dy;
	double& oneOverdx2 = grid.oneOverdx2;	double& oneOverdy2 = grid.oneOverdy2;

	int num_all_full_cells, nnz;
	Field2D();
	~Field2D();


	// Non-Ghost Grid
	Field2D(const Grid2D& ipGrid);
	Field2D(const double& ipXMin, const double& ipXmax, const int& ipiRes, const double& ipYMin, const double& ipYmax, const int& ipjRes);
	Field2D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);

	void Initialize(const Grid2D& ipGrid);
	void Initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);


	// Ghost Grid
	Field2D(const int& ipGhostWidth, const Grid2D& ipGrid);
	Field2D(const int& ipGhostWidth, const double& ipXMin, const double& ipXmax, const int& ipiRes, const double& ipYMin, const double& ipYmax, const int& ipjRes);
	Field2D(const int& ipGhostWidth, const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);

	void Initialize(const int& ipGhostWidth, const Grid2D& ipGrid);
	//void Initialize(const int& ipGhostWidth, double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes);




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

	//// Write MATLAB Variable
	inline void Variable(const char * varName);

	inline void WriteFile(const string& fileName);

	inline TT interpolation(const double& x, const double& y) const;
	inline TT interpolation(const Vector2D<double>& ipVector) const;
	
	inline void FillGhostCell();

	inline Vector2D<int> containedCell(const double& x, const double& y) const;

	inline Vector2D<double> Gradient(const int& i, const int& j);
	inline void Gradient();

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

	inline void CountNonZero();
	inline void CountNonZero(int& computedPtNum, int& computedEltNum);

	inline void SaveOld();
	inline void SaveOld(Array2D<TT>& copyArray);

	inline void Delete();
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
	dataArray.Delete();
	ghostDataArray.Delete();
	dataArrayOld.Delete();
	ghostDataArrayOld.Delete();
	BC.Delete();
	K1.Delete();
	K2.Delete();
	K3.Delete();
	dfdxM.Delete();
	dfdxP.Delete();
	dfdyM.Delete();
	dfdyP.Delete();
	gradient.Delete();
}

template<class TT>
inline Field2D<TT>::Field2D(const Grid2D & ipGrid)
	:grid(ipGrid)
{
	grid = ipGrid;
	Initialize(1, grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const double & ipXMin, const double & ipXmax, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjRes)
{
	grid.Initialize(ipXMin, ipXmax, 0, ipiRes, ipYMin, ipYmax, 0, ipjRes);
	Initialize(1, grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	grid.Initialize(ipXMin, ipXmax, ipiStart, ipiRes, ipYMin, ipYmax, ipjStart, ipjRes);
	Initialize(1, grid);
}

template<class TT>
inline void Field2D<TT>::Initialize(const Grid2D & ipGrid)
{
	Initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes, ipGrid.yMin, ipGrid.yMax, ipGrid.jStart, ipGrid.jRes);
	dataArray = Array2D<TT>(grid);
	dataArrayOld = Array2D<TT>(grid);
	ghostDataArray = Array2D<TT>(grid);
	ghostDataArrayOld = Array2D<TT>(grid);
	BC = Array2D<int>(grid);
	if (sizeof(TT) == 8)
	{
		K1 = Array2D<TT>(grid);
		K2 = Array2D<TT>(grid);
		K3 = Array2D<TT>(grid);
		dfdxM = Array2D<TT>(grid);
		dfdxP = Array2D<TT>(grid);
		dfdyM = Array2D<TT>(grid);
		dfdyP = Array2D<TT>(grid);
	}
}

template<class TT>
inline void Field2D<TT>::Initialize(const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
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

	Initialize(ipGhostWidth, grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const int & ipGhostWidth, const double & ipXMin, const double & ipXmax, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjRes)
{
	grid.Initialize(ipXMin, ipXmax, 0, ipiRes, ipYMin, ipYmax, 0, ipjRes);

	Initialize(ipGhostWidth, grid);
}

template<class TT>
inline Field2D<TT>::Field2D(const int & ipGhostWidth, const double & ipXMin, const double & ipXmax, const int & ipiStart, const int & ipiRes, const double & ipYMin, const double & ipYmax, const int & ipjStart, const int & ipjRes)
{
	grid.Initialize(ipXMin, ipXmax, ipiStart, ipiRes, ipYMin, ipYmax, ipjStart, ipjRes);

	Initialize(ipGhostWidth, grid);
}

template<class TT>
inline void Field2D<TT>::Initialize(const int & ipGhostWidth, const Grid2D & ipGrid)
{
	ghostWidth = ipGhostWidth;

	Initialize(ipGrid.xMin, ipGrid.xMax, ipGrid.iStart, ipGrid.iRes, ipGrid.yMin, ipGrid.yMax, ipGrid.jStart, ipGrid.jRes);

	double widthX = double(ghostWidth)*(ipGrid.xMax - ipGrid.xMin) / double(ipGrid.iRes - 1);
	double widthY = double(ghostWidth)*(ipGrid.yMax - ipGrid.yMin) / double(ipGrid.jRes - 1);

	ghostGrid.Initialize(ipGrid.xMin - widthX, ipGrid.xMax + widthX, ipGrid.iStart - ghostWidth, ipGrid.iRes + 2 * ghostWidth, ipGrid.yMin - widthY, ipGrid.yMax + widthY, ipGrid.jStart - ghostWidth, ipGrid.jRes + 2 * ghostWidth);
	ghostDataArray = Array2D<TT>(ghostGrid);

	innerGrid.Initialize(ipGrid.xMin + dx, ipGrid.xMax - dx, ipGrid.iStart + 1, ipGrid.iRes - 2, ipGrid.yMin + dy, ipGrid.yMax - dy, ipGrid.jStart + 1, ipGrid.jRes - 2);

	dataArray = Array2D<TT>(grid);
	dataArrayOld = Array2D<TT>(grid);
	ghostDataArray = Array2D<TT>(ghostGrid);
	ghostDataArrayOld = Array2D<TT>(ghostGrid);
	BC = Array2D<int>(grid);
	if (sizeof(TT) == 8)
	{
		K1 = Array2D<TT>(grid);
		K2 = Array2D<TT>(grid);
		K3 = Array2D<TT>(grid);
		dfdxM = Array2D<TT>(grid);
		dfdxP = Array2D<TT>(grid);
		dfdyM = Array2D<TT>(grid);
		dfdyP = Array2D<TT>(grid);
	}
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
	else if ((i >= iStartG && i <= iEndG) && (j >= jStartG && j <= jEndG))
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
	else if ((i >= iStartG && i <= iEndG) && (j >= jStartG && j <= jEndG))
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
	innerGrid = ipField.innerGrid;
	dataArray = ipField.dataArray;
	ghostDataArray = ipField.ghostDataArray;
	dataArrayOld = ipField.dataArrayOld;
	ghostDataArrayOld = ipField.ghostDataArrayOld;
	BC = ipField.BC;
	if (sizeof(TT) == 8)
	{
		K1 = ipField.K1;
		K2 = ipField.K2;
		K3 = ipField.K3;
		dfdxM = ipField.dfdxM;
		dfdxP = ipField.dfdxP;
		dfdyM = ipField.dfdyM;
		dfdyP = ipField.dfdyP;
	}

	ghostWidth = ipField.ghostWidth;
	//iRes = ipField.iRes;
	//jRes = ipField.jRes;
	//iStart = ipField.iStart;
	//jStart = ipField.jStart;
	//iEnd = iStart + iRes - 1;
	//jEnd = jStart + jRes - 1;
	//xMin = ipField.xMin;
	//yMin = ipField.yMin;
	//xMax = ipField.xMax;
	//yMax = ipField.yMax;
	//xLength = xMax - xMin;
	//yLength = yMax - yMin;
	//dx = xLength / (double)(iRes - 1);
	//dy = yLength / (double)(jRes - 1);
	//twodx = 2.0 * dx;
	//twody = 2.0*dy;
	//dx2 = dx*dx;
	//dy2 = dy*dy;
	//oneOverdx = 1.0 / dx;
	//oneOverdy = 1.0 / dy;
	//oneOver2dx = 1.0 / twodx;
	//oneOver2dy = 1.0 / twody;
	//oneOverdx2 = oneOverdx*oneOverdx;
	//oneOverdy2 = oneOverdy*oneOverdy;
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
inline void Field2D<TT>::Gradient()
{
	gradient = Array2D<Vector2D<double>>(grid);
#pragma omp parallel for
	for (int i = gradient.iStart; i <= gradient.iEnd; i++)
	{
		for (int j = gradient.jStart; j <= gradient.jEnd; j++)
		{
			gradient(i, j).i = dxPhi(i, j);
			gradient(i, j).j = dyPhi(i, j);
		}
	}
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

template<class TT>
inline void Field2D<TT>::CountNonZero()
{
	num_all_full_cells = 0;
	nnz = 0;
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			if (BC(i, j) < 0)
			{
				continue;
			}
			num_all_full_cells++;
			nnz++;
			if (i > iStart)
			{
				if (BC(i - 1, j) > -1)
				{
					nnz++;
				}
			}
			if (i < iEnd)
			{
				if (BC(i + 1, j) > -1)
				{
					nnz++;
				}
			}
			if (j > jStart)
			{
				if (BC(i, j - 1) > -1)
				{
					nnz++;
				}
			}
			if (j < jEnd)
			{
				if (BC(i, j + 1) > -1)
				{
					nnz++;
				}
			}

		}
	}
}

template<class TT>
inline void Field2D<TT>::CountNonZero(int & computedPtNum, int & computedEltNum)
{
	int computedPtNum = 0;
	int computedEltNum = 0;
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			computedPtNum++;
			computedEltNum++;
			if (BC(i, j) < 0)
			{
				continue;
			}

			if (i > iStart)
			{
				if (BC(i - 1, j) > -1)
				{
					computedEltNum++;
				}
			}
			if (i < iEnd)
			{
				if (BC(i + 1, j) > -1)
				{
					computedEltNum++;
				}
			}
			if (j > jStart)
			{
				if (BC(i, j - 1) > -1)
				{
					computedEltNum++;
				}
			}
			if (j < jEnd)
			{
				if (BC(i, j + 1) > -1)
				{
					computedEltNum++;
				}
			}
		}

	}
}

template<class TT>
inline void Field2D<TT>::SaveOld()
{
	dataArrayOld = dataArray;
}

template<class TT>
inline void Field2D<TT>::Delete()
{
	dataArray.Delete();
	ghostDataArray.Delete();
	dataArrayOld.Delete();
	ghostDataArrayOld.Delete();
	K1.Delete();
	K2.Delete();
	K3.Delete();
	dfdxM.Delete();
	dfdxP.Delete();
	dfdyM.Delete();
	dfdyP.Delete();
	gradient.Delete();
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
