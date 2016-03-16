#pragma once


#include "AdvectionMethod2D.h"


class EulerianFluidSolver2D
{
public:
	Grid2D grid;
	Grid2D cellGrid;

	LevelSet2D levelSet;

	Field2D<double> pressure;

	Field2D<Vector2D<double>> velocity;

	double reynoldNum;
	double dt;
	double cflCondition;

	int ghostWidth;

	EulerianFluidSolver2D();
	~EulerianFluidSolver2D();

	void InitialCondition(const int& example);
	void FluidSolver(const int& example);
private:

};




EulerianFluidSolver2D::EulerianFluidSolver2D()
{
}

EulerianFluidSolver2D::~EulerianFluidSolver2D()
{
}

inline void EulerianFluidSolver2D::InitialCondition(const int & example)
{
	if (example==1)
	{
		cout << "Cavity Flow" << endl;

		grid = Grid2D(0, 1.0, 101, 0, 1, 101);
		cellGrid = Grid2D(grid.xMin + grid.dx/2, grid.xMax - grid.dx/2, grid.iRes - 1, grid.yMin + grid.dy/2, grid.yMax - grid.dy/2, grid.jRes - 1);

		levelSet = LevelSet2D(grid);
		pressure = Field2D<double>(cellGrid);
		velocity = Field2D<Vector2D<double>>(grid);

		reynoldNum = 500;
		cflCondition = 0.1;
		

	}

	if (example==2)
	{

	}
}

inline void EulerianFluidSolver2D::FluidSolver(const int & example)
{
	InitialCondition(example);
}
