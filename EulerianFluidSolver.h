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

	int maxIteration;
	int writeOutputIteration;

	EulerianFluidSolver2D();
	~EulerianFluidSolver2D();

	void InitialCondition(const int& example);
	void FluidSolver(const int& example, const int& timeSteppingOrder);

	void TimeAdvanceForwardEuler(Field2D<Vector2D<double>>& ipVelocity);

	void OutputResult(const int& iter);
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
		
		maxIteration = 1000;
		writeOutputIteration = 1;
	}

	if (example==2)
	{

	}
}

inline void EulerianFluidSolver2D::FluidSolver(const int & example, const int& timeSteppingOrder)
{
	InitialCondition(example);

	//OutputResult(0);

	for (int i = 0; i < maxIteration; i++)
	{
		if (timeSteppingOrder==1)
		{
			TimeAdvanceForwardEuler(velocity);
		}

		if (i%writeOutputIteration == 0)
		{
			//OutputResult(i);
		}
	}
}

inline void EulerianFluidSolver2D::TimeAdvanceForwardEuler(Field2D<Vector2D<double>>& ipVelocity)
{
	ipVelocity.FillGhostCell();

	Field2D<Vector2D<double>> originVelocity = ipVelocity;

//#pragma omp parallel for
	for (int i = ipVelocity.iStart; i <= ipVelocity.iEnd; i++)
	{
		for (int j = ipVelocity.jStart; j <= ipVelocity.jEnd; j++)
		{
			double du2dx = (originVelocity(i + 1, j)(0)*originVelocity(i + 1, j)(0) - originVelocity(i - 1, j)(0)*originVelocity(i - 1, j)(0))*originVelocity.oneOver2dx;
			double duvdy = (originVelocity(i, j + 1)(0)*originVelocity(i, j + 1)(1) - originVelocity(i, j - 1)(0)*originVelocity(i, j - 1)(1))*originVelocity.oneOver2dy;
			double d2udx2 = (originVelocity(i + 1, j)(0) - 2 * originVelocity(i, j)(0) + originVelocity(i - 1, j)(0))*originVelocity.oneOverdx2;
			double d2vdy2 = (originVelocity(i, j + 1)(1) - 2 * originVelocity(i, j)(1) + originVelocity(i, j - 1)(1))*originVelocity.oneOverdy2;
			double temp = originVelocity(i, j)(0) + dt*(-du2dx - duvdy + 1 / reynoldNum*(d2udx2 + d2vdy2));

			ipVelocity(i, j)(0);
			cout << temp << endl;
			cout << ipVelocity(i, j)(0) << endl;
			cout << endl;
		}
	}

}

inline void EulerianFluidSolver2D::OutputResult(const int & iter)
{
	cout << "Write results" << endl;

	ofstream solutionFile1;
	solutionFile1.open("D:\\Data/velocity" + to_string(iter) + ".txt", ios::binary);
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			solutionFile1 << i << " " << j << " " << grid(i, j) << " " << velocity(i, j) << endl;
		}
	}
	solutionFile1.close();

	ofstream solutionFile2;
	solutionFile2.open("D:\\Data/pressure" + to_string(iter) + ".txt", ios::binary);
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			solutionFile2 << i << " " << j << " " << grid(i, j) << " " << pressure(i, j) << endl;
		}
	}
	solutionFile2.close();
}
