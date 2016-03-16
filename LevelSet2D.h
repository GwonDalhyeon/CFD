#pragma once


//#ifndef LevelSet2D_H
//#define LevelSet2D_H
#include "CommonDef.h"
#include "Vector2D.h"
#include "VectorND.h"
#include "Grid2D.h"
#include "Field2D.h"



class LevelSet2D
{
public:
	Grid2D grid;

	Field2D<double> phi;
	Field2D<Vector2D<double>> normal;
	Field2D<Vector2D<double>> unitNormal;
	Field2D<Vector2D<double>> tangential;
	Field2D<double> meanCurvature;


	LevelSet2D();
	~LevelSet2D();

	LevelSet2D(const Grid2D& ipGrid);

	const int index(const int& i, const int& j) const;

	const int index(const Vector2D<int>& ipVector) const;

	inline double& operator ()(const int& i, const int& j) const;

	inline double& operator ()(const Vector2D<int>& ipVector) const;

	inline double operator ()(const double& x, const double& y) const;

	inline double operator ()(const Vector2D<double>& ipVector) const;

	inline void computeNormal();
	inline Vector2D<double> computeNormal(const int& i, const int& j);
	inline Vector2D<double> computeNormal(const Vector2D<int> ipVector);

	inline void computeUnitNormal();
	inline Vector2D<double> computeUnitNormal(const int& i, const int& j);
	inline Vector2D<double> computeUnitNormal(const Vector2D<int> ipVector);

	inline void computeMeanCurvature();
	inline double computeMeanCurvature(const int& i, const int& j);
	inline double computeMeanCurvature(const Vector2D<int> & ipVector);

	inline Vector2D<double> gradient(const int& i, const int& j);

	inline double interpolation(const double& x, const double& y);
	inline double interpolation(const Vector2D<double>& ipVector);
	inline double interpolation(const int& i, const int& j);
	//inline double interpolation(const Vector2D<int>& ipVector);

	// Derivative
	inline double dxxPhi(const int& i, const int& j);
	inline double dxPhi(const int& i, const int& j);
	inline double dxPlusPhi(const int& i, const int& j);
	inline double dxMinusPhi(const int& i, const int& j);

	inline double dyyPhi(const int& i, const int& j);
	inline double dyPhi(const int& i, const int& j);
	inline double dyPlusPhi(const int& i, const int& j);
	inline double dyMinusPhi(const int& i, const int& j);

	inline double dxyPhi(const int& i, const int& j);

private:

};


//#endif // !LevelSet2D



LevelSet2D::LevelSet2D()
{
}

LevelSet2D::~LevelSet2D()
{
}


inline LevelSet2D::LevelSet2D(const Grid2D & ipGrid)
{
	grid = ipGrid;
	phi = Field2D<double>(ipGrid);
	normal = Field2D<Vector2D<double>>(ipGrid);
	unitNormal = Field2D<Vector2D<double>>(ipGrid);
	tangential = Field2D<Vector2D<double>>(ipGrid);
	meanCurvature = Field2D<double>(ipGrid);
}

const int LevelSet2D::index(const int & i, const int & j) const
{
	assert(i >= phi.iStart && i <= phi.iEnd);
	assert(j >= phi.jStart && j <= phi.jEnd);
	return phi.index(i, j);
}

const int LevelSet2D::index(const Vector2D<int>& ipVector) const
{
	assert(ipVector[0] >= phi.iStart && ipVector[0] <= phi.iEnd);
	assert(ipVector[1] >= phi.jStart && ipVector[1] <= phi.jEnd);
	return phi.index(ipVector);
}

inline double & LevelSet2D::operator()(const int & i, const int & j) const
{
	assert(i >= phi.iStart && i <= phi.iEnd);
	assert(j >= phi.jStart && j <= phi.jEnd);
	return phi(i, j);
}


inline double & LevelSet2D::operator()(const Vector2D<int>& ipVector) const
{
	assert(ipVector[0] >= phi.iStart && ipVector[0] <= phi.iEnd);
	assert(ipVector[1] >= phi.jStart && ipVector[1] <= phi.jEnd);
	return phi(ipVector);
}



inline double LevelSet2D::operator()(const double & x, const double & y) const
{
	assert(x >= grid.xMin && x <= grid.xMax);
	assert(y >= grid.yMin && y <= grid.yMax);

	return phi(x, y);
}

inline double LevelSet2D::operator()(const Vector2D<double>& ipVector) const
{
	assert(ipVector[0] >= grid.xMin && ipVector[0] <= grid.xMax);
	assert(ipVector[1] >= grid.yMin && ipVector[1] <= grid.yMax);
	return phi(ipVector.x, ipVector.y);
}

inline void LevelSet2D::computeNormal()
{
#pragma omp parallel for
	for (int j = phi.jStart; j <= phi.jEnd; j++)
	{
		for (int i = phi.iStart; i <= phi.iEnd; i++)
		{
			normal.dataArray(i, j) = computeNormal(i, j);
		}
	}
}

inline Vector2D<double> LevelSet2D::computeNormal(const int & i, const int & j)
{
	Vector2D<double> normal;

	if (j > phi.jStart && j < phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j))*phi.oneOver2dx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1))*phi.oneOver2dy;
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1))*phi.oneOver2dy;
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1))*phi.oneOver2dy;
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}
	}
	else if (j == phi.jStart)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j))*phi.oneOver2dx;;
			normal.values[1] = (phi(i, j + 1) - phi(i, j))*phi.oneOverdy;
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j))*phi.oneOverdy;
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j + 1) - phi(i, j))*phi.oneOverdy;
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else if (j == phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j))*phi.oneOver2dx;;
			normal.values[1] = (phi(i, j) - phi(i, j - 1))*phi.oneOverdy;
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j) - phi(i, j - 1))*phi.oneOverdy;
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j))*phi.oneOverdx;
			normal.values[1] = (phi(i, j) - phi(i, j - 1))*phi.oneOverdy;
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else
	{
		assert(j >= phi.jStart && j <= phi.jEnd);
	}


	return normal;
}

inline Vector2D<double> LevelSet2D::computeNormal(const Vector2D<int> ipVector)
{
	return computeNormal(ipVector[0], ipVector[1]);
}



inline void LevelSet2D::computeUnitNormal()
{
#pragma omp parallel for
	for (int j = phi.jStart; j <= phi.jEnd; j++)
	{
		for (int i = phi.iStart; i <= phi.iEnd; i++)
		{
			unitNormal.dataArray(i, j) = computeUnitNormal(i, j);
		}
	}
}

inline Vector2D<double> LevelSet2D::computeUnitNormal(const int & i, const int & j)
{
	Vector2D<double> normal;

	if (j > phi.jStart && j < phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j)) / (abs(phi(i + 1, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1)) / (abs(phi(i, j + 1) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j)) / (abs(phi(i + 1, j) - phi(i, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1)) / (abs(phi(i, j + 1) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j)) / (abs(phi(i, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j - 1)) / (abs(phi(i, j + 1) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}
	}
	else if (j == phi.jStart)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j)) / (abs(phi(i + 1, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j)) / (abs(phi(i, j + 1) - phi(i, j)) + DBL_EPSILON);
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j)) / (abs(phi(i + 1, j) - phi(i, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j)) / (abs(phi(i, j + 1) - phi(i, j)) + DBL_EPSILON);
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j)) / (abs(phi(i, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j + 1) - phi(i, j)) / (abs(phi(i, j + 1) - phi(i, j)) + DBL_EPSILON);
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else if (j == phi.jEnd)
	{

		if (i > phi.iStart && i < phi.iEnd)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i - 1, j)) / (abs(phi(i + 1, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j) - phi(i, j - 1)) / (abs(phi(i, j) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iStart)
		{
			normal.values[0] = (phi(i + 1, j) - phi(i, j)) / (abs(phi(i + 1, j) - phi(i, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j) - phi(i, j - 1)) / (abs(phi(i, j) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else if (i == phi.iEnd)
		{
			normal.values[0] = (phi(i, j) - phi(i - 1, j)) / (abs(phi(i, j) - phi(i - 1, j)) + DBL_EPSILON);
			normal.values[1] = (phi(i, j) - phi(i, j - 1)) / (abs(phi(i, j) - phi(i, j - 1)) + DBL_EPSILON);
		}
		else
		{
			cout << "Level set unitNormal error." << endl;
		}

	}
	else
	{
		assert(j >= phi.jStart && j <= phi.jEnd);
	}


	return normal;
}

inline Vector2D<double> LevelSet2D::computeUnitNormal(const Vector2D<int> ipVector)
{
	return computeUnitNormal(ipVector[0], ipVector[1]);
}

inline void LevelSet2D::computeMeanCurvature()
{
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.iStart; j <= grid.jEnd; j++)
		{
			meanCurvature(i, j) = computeMeanCurvature(i, j);
		}
	}
}

inline double LevelSet2D::computeMeanCurvature(const int & i, const int & j)
{
	return -(dxxPhi(i, j)*dyPhi(i, j)*dyPhi(i, j) - 2.0*dxyPhi(i, j)*dxPhi(i, j)*dyPhi(i, j) + dyyPhi(i, j)*dxPhi(i, j)*dxPhi(i, j)) / pow(dxPhi(i, j)*dxPhi(i, j) + dyPhi(i, j)*dyPhi(i, j) + DBL_EPSILON, 3.0 / 2.0);
}

inline double LevelSet2D::computeMeanCurvature(const Vector2D<int> & ipVector)
{
	return computeMeanCurvature(ipVector.i, ipVector.j);
}

inline Vector2D<double> LevelSet2D::gradient(const int & i, const int & j)
{
	return Vector2D<double>(dxPhi(i, j), dyPhi(i, j));
}

inline double LevelSet2D::interpolation(const double & x, const double & y)
{
	if (grid.xMin <= x && grid.xMax >= x &&grid.yMin <= y && grid.yMax >= y)
	{
		Vector2D<double> xy(x, y);
		Vector2D<int> cell = phi.containedCell(x, y);

		double distance00, distance10, distance01, distance11;
		distance00 = (grid.point(cell.i, cell.j) - xy).magnitude();
		distance10 = (grid.point(cell.i + 1, cell.j) - xy).magnitude();
		distance01 = (grid.point(cell.i, cell.j + 1) - xy).magnitude();
		distance11 = (grid.point(cell.i + 1, cell.j + 1) - xy).magnitude();
		if (distance00 < grid.dx / 2)
		{
			return phi(cell);
		}
		if (distance10 < grid.dx / 2)
		{
			return phi(cell.i + 1, cell.j);
		}
		if (distance01 < grid.dx / 2)
		{
			return phi(cell.i, cell.j + 1);
		}
		if (distance11<grid.dx / 2)
		{
			return phi(cell.i + 1, cell.j + 1);
		}

		return ((phi(cell)*distance00 + phi(cell.i + 1, cell.j)*distance10 + phi(cell.i, cell.j + 1)*distance01 + phi(cell.i + 1, cell.j + 1)*distance11) / (distance00 + distance01 + distance10 + distance11));
	}
	else if (grid.xMin <= x && grid.xMax >= x &&grid.yMin > y)
	{
		double value1 = interpolation(x, grid.yMin);
		double value2 = interpolation(x, grid.yMin + grid.dy);

		return (value1 - value2)*grid.oneOverdy*(grid.yMin - y) + value1;
	}
	else if (grid.xMin <= x && grid.xMax >= x &&grid.yMax < y)
	{
		double value1 = interpolation(x, grid.yMax);
		double value2 = interpolation(x, grid.yMax - grid.dy);

		return (value1 - value2)*grid.oneOverdy*(y - grid.yMax) + value1;
	}
	else if (grid.xMin > x  &&grid.yMin <= y && grid.yMax >= y)
	{
		double value1 = interpolation(grid.xMin, y);
		double value2 = interpolation(grid.xMin + grid.dx, y);

		return (value1 - value2)*grid.oneOverdx*(grid.xMin - x) + value1;
	}
	else if (grid.xMax < x &&grid.yMin <= y && grid.yMax >= y)
	{
		double value1 = interpolation(grid.xMax, y);
		double value2 = interpolation(grid.xMax - grid.dx, y);

		return (value1 - value2)*grid.oneOverdx*(x - grid.xMax) + value1;
	}
	else
	{
		double tempX;
		double tempY;
		if (x < grid.xMin)
		{
			tempX = grid.xMin;
		}
		else
		{
			tempX = grid.xMax;
		}

		if (y < grid.yMin)
		{
			tempY = grid.yMin;
		}
		else
		{
			tempY = grid.yMax;
		}

		double value1 = interpolation(tempX, y);
		double value2 = interpolation(x, tempY);

		return (value1 + value2) / 2.0;
	}
}

inline double LevelSet2D::interpolation(const Vector2D<double>& ipVector)
{
	return interpolation(ipVector.x, ipVector.y);
}

inline double LevelSet2D::interpolation(const int & i, const int & j)
{
	return interpolation(grid.point(i, j));
}

inline double LevelSet2D::dxxPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (phi(i + 1, j) - 2.0 * phi(i, j) + phi(i - 1, j))*grid.oneOverdx2;
	}
	else if (i == grid.iStart)
	{
		double tempPhi = interpolation(i - 1, j);
		return (phi(i + 1, j) - 2.0 * phi(i, j) + tempPhi)*grid.oneOverdx2;
		//return (phi(i, j) - 2 * phi(i + 1, j) + phi(i + 2, j))*grid.oneOver2dx;
	}
	else
	{
		double tempPhi = interpolation(i + 1, j);
		return (tempPhi - 2.0 * phi(i, j) + phi(i - 1, j))*grid.oneOverdx2;
		//return (phi(i - 2, j) - 2 * phi(i - 1, j) + phi(i, j))*grid.oneOver2dx;
	}
}

inline double LevelSet2D::dxPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart && i < grid.iEnd)
	{
		return (phi(i + 1, j) - phi(i - 1, j))*grid.oneOver2dx;
	}
	else if (i == grid.iStart)
	{
		double tempPhi = interpolation(i - 1, j);
		return (phi(i + 1, j) - tempPhi)*grid.oneOver2dx;
		//return dxPlusPhi(i, j);
	}
	else
	{
		double tempPhi = interpolation(i + 1, j);
		return (tempPhi - phi(i - 1, j))*grid.oneOver2dx;
		//return dxMinusPhi(i, j);
	}
}

inline double LevelSet2D::dxPlusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i < grid.iEnd)
	{
		return (phi(i + 1, j) - phi(i, j))*grid.oneOverdx;
	}
	else
	{
		double tempPhi = interpolation(i + 1, j);
		return (tempPhi - phi(i, j))*grid.oneOverdx;
		//return (phi(i, j) - phi(i - 1, j))*grid.oneOverdx;
	}
}

inline double LevelSet2D::dxMinusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (i > grid.iStart)
	{
		return (phi(i, j) - phi(i - 1, j))*grid.oneOverdx;
	}
	else
	{
		double tempPhi = interpolation(i - 1, j);
		return (phi(i, j) - tempPhi)*grid.oneOverdx;
		//return (phi(i + 1, j) - phi(i, j))*grid.oneOverdx;
	}
}

inline double LevelSet2D::dyyPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (phi(i, j + 1) - 2 * phi(i, j) + phi(i, j - 1))*grid.oneOverdy2;
	}
	else if (j == grid.jStart)
	{
		return (phi(i, j) - 2 * phi(i, j + 1) + phi(i, j + 2))*grid.oneOverdy2;
	}
	else
	{
		return (phi(i, j - 2) - 2 * phi(i, j - 1) + phi(i, j))*grid.oneOverdy2;
	}
}

inline double LevelSet2D::dyPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart && j < grid.jEnd)
	{
		return (phi(i, j + 1) - phi(i, j - 1))*grid.oneOver2dy;
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

inline double LevelSet2D::dyPlusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j < grid.jEnd)
	{
		return (phi(i, j + 1) - phi(i, j))*grid.oneOverdy;
	}
	else
	{
		return (phi(i, j) - phi(i, j - 1))*grid.oneOverdy;
	}
}

inline double LevelSet2D::dyMinusPhi(const int & i, const int & j)
{
	assert(i >= grid.iStart && i <= grid.iEnd);
	assert(j >= grid.jStart && j <= grid.jEnd);

	if (j > grid.jStart)
	{
		return (phi(i, j) - phi(i, j - 1))*grid.oneOverdy;
	}
	else
	{
		return (phi(i, j + 1) - phi(i, j))*grid.oneOverdy;
	}
}

inline double LevelSet2D::dxyPhi(const int & i, const int & j)
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
