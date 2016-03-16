#pragma once
#include "GridInfomation.h"

//#ifndef LevelSet_H
//#define LevelSet_H


class LevelSet
{
public:
	GridInfo grid;
	double* phi;

	LevelSet();
	LevelSet(const GridInfo& inputGrid);
	~LevelSet();

	int index(int i, int j);
	int index(int i, int j, int k);

	double unitNormal(int node);
	void unitNormal(int node1, int node2, double* normal);

	LevelSet& operator=(const LevelSet& originLevelSet);
private:

};

inline LevelSet::LevelSet()
{
}

LevelSet::LevelSet(const GridInfo& inputGrid)
{
	grid = inputGrid;
	//grid = GridInfo(inputGrid.X0, inputGrid.X1, inputGrid.numX);

	if (grid.dimension == 1)
	{
		phi = new double[grid.numX];
	}
	if (grid.dimension == 2)
	{
		phi = new double[grid.numX*grid.numY];
	}
	if (grid.dimension == 3)
	{
		phi = new double[grid.numX*grid.numY*grid.numZ];
	}
}

LevelSet::~LevelSet()
{
	//delete[]grid;
	delete[] phi;

}

inline int LevelSet::index(int i, int j)
{
	return i + j*grid.numX;
}

inline int LevelSet::index(int i, int j, int k)
{
	return i + j*grid.numX + k*grid.numX*grid.numY;
}

inline double LevelSet::unitNormal(int node)
{
	int i = node;
	if (i<grid.numX - 1 && i>0)
	{
		return (phi[i + 1] - phi[i - 1]) / abs(phi[i + 1] - phi[i - 1]);
	}
	else if (i == 0)
	{
		return (phi[i + 1] - phi[i]) / abs(phi[i + 1] - phi[i]);
	}
	else if (i == grid.numX - 1)
	{
		return (phi[i] - phi[i - 1]) / abs(phi[i] - phi[i - 1]);
	}
	return 0;
}

inline void LevelSet::unitNormal(int node1, int node2, double* normal)
{
	int i = node1;
	int j = node2;
	normal[0] = 0;
	normal[1] = 0;

	if (i<grid.numX - 1 && i>0)
	{
		if (j >= 0 && j<grid.numY)
		{
			normal[0] = (phi[index(i + 1, j)] - phi[index(i - 1, j)]) / (abs(phi[index(i + 1, j)] - phi[index(i - 1, j)]) + DBL_EPSILON);
		}
		else if (j < 0)
		{
			normal[0] = (phi[index(i + 1, 0)] - phi[index(i - 1, 0)]) / (abs(phi[index(i + 1, 0)] - phi[index(i - 1, 0)]) + DBL_EPSILON);
		}
		else if (j >= grid.numY)
		{
			normal[0] = (phi[index(i + 1, grid.numY - 1)] - phi[index(i - 1, grid.numY - 1)]) / (abs(phi[index(i + 1, grid.numY - 1)] - phi[index(i - 1, grid.numY - 1)]) + DBL_EPSILON);
		}
	}
	else if (i == 0)
	{
		if (j<grid.numY && j >= 0)
		{
			normal[0] = (phi[index(i + 1, j)] - phi[index(i, j)]) / (abs(phi[index(i + 1, j)] - phi[index(i, j)]) + DBL_EPSILON);
		}
		else if (j < 0)
		{
			normal[0] = (phi[index(i + 1, 0)] - phi[index(i, 0)]) / (abs(phi[index(i + 1, 0)] - phi[index(i, 0)]) + DBL_EPSILON);
		}
		else if (j >= grid.numY)
		{
			normal[0] = (phi[index(i + 1, grid.numY - 1)] - phi[index(i, grid.numY - 1)]) / (abs(phi[index(i + 1, grid.numY - 1)] - phi[index(i, grid.numY - 1)]) + DBL_EPSILON);
		}
	}
	else if (i == grid.numX - 1)
	{
		if (j >= 0 && j<grid.numY)
		{
			normal[0] = (phi[index(i, j)] - phi[index(i - 1, j)]) / (abs(phi[index(i, j)] - phi[index(i - 1, j)]) + DBL_EPSILON);
		}
		else if (j < 0)
		{
			normal[0] = (phi[index(i, 0)] - phi[index(i - 1, 0)]) / (abs(phi[index(i, 0)] - phi[index(i - 1, 0)]) + DBL_EPSILON);
		}
		else if (j >= grid.numY)
		{
			normal[0] = (phi[index(i, grid.numY - 1)] - phi[index(i - 1, grid.numY - 1)]) / (abs(phi[index(i, grid.numY - 1)] - phi[index(i - 1, grid.numY - 1)]) + DBL_EPSILON);
		}
	}

	if (j<grid.numY - 1 && j>0)
	{
		if (i<grid.numX  && i >= 0)
		{
			normal[1] = (phi[index(i, j + 1)] - phi[index(i, j - 1)]) / (abs(phi[index(i, j + 1)] - phi[index(i, j - 1)]) + DBL_EPSILON);
		}
		else if (i<0)
		{
			normal[1] = (phi[index(0, j + 1)] - phi[index(0, j - 1)]) / (abs(phi[index(0, j + 1)] - phi[index(0, j - 1)]) + DBL_EPSILON);
		}
		else if (i >= grid.numX)
		{
			normal[1] = (phi[index(grid.numX - 1, j + 1)] - phi[index(grid.numX - 1, j - 1)]) / (abs(phi[index(grid.numX - 1, j + 1)] - phi[index(grid.numX - 1, j - 1)]) + DBL_EPSILON);
		}
	}
	else if (j == 0)
	{
		if (i<grid.numX  && i >= 0)
		{
			normal[1] = (phi[index(i, j + 1)] - phi[index(i, j)]) / (abs(phi[index(i, j + 1)] - phi[index(i, j)]) + DBL_EPSILON);
		}
		else if (i<0)
		{
			normal[1] = (phi[index(0, j + 1)] - phi[index(0, j)]) / (abs(phi[index(0, j + 1)] - phi[index(0, j)]) + DBL_EPSILON);
		}
		else if (i >= grid.numX)
		{
			normal[1] = (phi[index(grid.numX - 1, j + 1)] - phi[index(grid.numX - 1, j)]) / (abs(phi[index(grid.numX - 1, j + 1)] - phi[index(grid.numX - 1, j)]) + DBL_EPSILON);
		}
	}
	else if (j == grid.numY - 1)
	{
		if (i<grid.numX  && i >= 0)
		{
			normal[1] = (phi[index(i, j)] - phi[index(i, j - 1)]) / (abs(phi[index(i, j)] - phi[index(i, j - 1)]) + DBL_EPSILON);
		}
		else if (i<0)
		{
			normal[1] = (phi[index(0, j)] - phi[index(0, j - 1)]) / (abs(phi[index(0, j)] - phi[index(0, j - 1)]) + DBL_EPSILON);
		}
		else if (i >= grid.numX)
		{
			normal[1] = (phi[index(grid.numX - 1, j)] - phi[index(grid.numX - 1, j - 1)]) / (abs(phi[index(grid.numX - 1, j)] - phi[index(grid.numX - 1, j - 1)]) + DBL_EPSILON);
		}
	}
}

inline LevelSet & LevelSet::operator=(const LevelSet & originLevelSet)
{
	grid = originLevelSet.grid;

	if (phi != nullptr)
	{
		delete[] phi;
	}

	if (grid.dimension == 1)
	{
		phi = new double[grid.numX];
		for (int i = 0; i < grid.numX; i++)
		{
			phi[i] = originLevelSet.phi[i];
		}
	}
	else if (grid.dimension == 2)
	{
		phi = new double[grid.numX*grid.numY];
		for (int i = 0; i < grid.numX*grid.numY; i++)
		{
			phi[i] = originLevelSet.phi[i];
		}
	}
	else
	{
		phi = new double[grid.numX*grid.numY*grid.numZ];
		for (int i = 0; i < grid.numX*grid.numY*grid.numZ; i++)
		{
			phi[i] = originLevelSet.phi[i];
		}
	}


	return *this;
}




//#endif // !LevelSet_H
