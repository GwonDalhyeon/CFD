#pragma once

#include "AdvectionMethod2D.h"
#include "Grid2D.h"
#include "Field2D.h"
#include "LevelSet2D.h"
#include "LinearSolver.h"


class VortexSheet
{
public:
	Grid2D grid;

	LevelSet2D levelSet;

	Field2D<double> P;
	Field2D<double> streamFunction;
	Field2D<double> velocityX;
	Field2D<double> velocityY;

	double dt;
	double cflCondition;

	int ghostWidth;

	int maxIteration;
	int writeOutputIteration;

	VortexSheet();
	~VortexSheet();

	inline void InitialCondition(const int& example);
	inline void VortexSolver(const int& example, const int& timeSteppingOrder);

	inline void GenerateLinearSystem(Array2D<double>& matrixA);
	inline void GenerateLinearSystem(const Field2D<double>& u, VectorND<double>& vectorB);
	inline void Stream2Velocity();
	inline double AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2);
	inline double DeltaFt(const double& ip);
private:

};

VortexSheet::VortexSheet()
{
}

VortexSheet::~VortexSheet()
{
}

inline void VortexSheet::InitialCondition(const int & example)
{
	if (example==1)
	{
		grid = Grid2D(-1, 1, 101, -1, 1, 101);
		levelSet = LevelSet2D(grid);
		P = Field2D<double>(grid);
		streamFunction = Field2D<double>(grid);
		velocityX = Field2D<double>(grid);
		velocityY = Field2D<double>(grid);

		cflCondition = 0.1;

		maxIteration = 1000;
		writeOutputIteration = 1;

#pragma omp parallel for 
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = grid(i, j).y + 0.05*sin(PI*grid(i, j).x);
				P(i, j) = DeltaFt(levelSet(i, j));
			}
		}

	}
	else if (example ==2)
	{
		grid = Grid2D(-1, 1, 101, -1, 1, 101);
		levelSet = LevelSet2D(grid);
		P = Field2D<double>(grid);
		streamFunction = Field2D<double>(grid);
		velocityX = Field2D<double>(grid);
		velocityY = Field2D<double>(grid);
		
		cflCondition = 0.1;
		
		maxIteration = 1000;
		writeOutputIteration = 1;

		double eps = grid.dx;

#pragma omp parallel for 
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = grid(i, j).y / (1 - 0.75*sin(PI*grid(i, j).x));
				if (abs(levelSet(i,j))<8*grid.dx)
				{
					P(i, j) = -PI / (2 * eps*eps)*sin(PI*levelSet(i, j) / eps);
				}
				else
				{
					P(i, j) = 0;
				}
			}
		}
	}
}

inline void VortexSheet::VortexSolver(const int & example, const int & timeSteppingOrder)
{
	string fileName;
	InitialCondition(example);
	fileName = "phi0";
	levelSet.phi.WriteFile(fileName);

	Array2D<double> poissonMatrix(1, (grid.iRes - 2)*(grid.jRes - 2), 1, (grid.iRes - 2)*(grid.jRes - 2));
	GenerateLinearSystem(poissonMatrix);
	cout << "Start CSR." << endl;
	CSR<double> poissonCSR(poissonMatrix);
	cout << "End CSR." << endl;

	VectorND<double> vectorB((grid.iRes - 2)*(grid.jRes - 2));


	VectorND<double> stream((grid.iRes - 2)*(grid.jRes - 2));

	int idx;
	int innerIRes = grid.iRes - 2;


	for (int i = 1; i <= maxIteration; i++)
	{

		GenerateLinearSystem(P, vectorB);
		stream = CG<double>(poissonMatrix, vectorB);


#pragma omp parallel for private(idx)
		for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
		{
			for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
			{
				idx = (i - 1) + (j - 1)*innerIRes;
				P(i, j) = stream(idx);
			}
		}

		Stream2Velocity();
		dt = AdaptiveTimeStep(velocityX, velocityY);
		AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, velocityX, velocityY, dt);
		if (i%writeOutputIteration==0)
		{
			fileName = "velocityX" + to_string(i);
			velocityX.WriteFile(fileName);

			fileName = "velocityY" + to_string(i);
			velocityY.WriteFile(fileName);

			fileName = "phi" + to_string(i);
			levelSet.phi.WriteFile(fileName);

			fileName = "stream" + to_string(i);
			P.WriteFile(fileName);
		}
	}
}



inline void VortexSheet::GenerateLinearSystem(Array2D<double>& matrixA)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int index, leftIndex, rightIndex, bottomIndex, topIndex;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

#pragma omp parallel for private(index, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			index = (i - 1)*innerIRes*innerJRes + (i - 1) + (j - 1)*innerIRes*(innerIRes*innerJRes + 1);
			leftIndex = (i - 1)*innerIRes*innerJRes + (i - 1 - 1) + (j - 1)*innerIRes*(innerIRes*innerJRes + 1);
			rightIndex = (i - 1)*innerIRes*innerJRes + i + (j - 1)*innerIRes*(innerIRes*innerJRes + 1);
			bottomIndex = (i - 1)*innerIRes*innerJRes + (i - 1) + (j - 1)*innerIRes*innerIRes*innerJRes + (j - 1 - 1)*innerIRes;
			topIndex = (i - 1)*innerIRes*innerJRes + (i - 1) + (j - 1)*innerIRes*innerIRes*innerJRes + (j)*innerIRes;

			matrixA(index) = +2 / grid.dx2 + 2 / grid.dy2;

			if (i>grid.iStart + 1)
			{
				matrixA(leftIndex) = 1 / grid.dx2;
			}
			if (i<grid.iEnd - 1)
			{
				matrixA(rightIndex) = 1 / grid.dx2;
			}
			if (j>grid.jStart + 1)
			{
				matrixA(bottomIndex) = 1 / grid.dy2;
			}
			if (j<grid.jEnd - 1)
			{
				matrixA(topIndex) = 1 / grid.dy2;
			}
		}
	}
	cout << "End Generate Linear System : matrix A" << endl;
}


inline void VortexSheet::GenerateLinearSystem(const Field2D<double>& u, VectorND<double>& vectorB)
{
	int index;
	int innerIRes = grid.iRes - 2;


#pragma omp parallel for private(index)
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			index = (i - 1) + (j - 1)*innerIRes;

			vectorB(index) = -P(i, j);

			//if (i == grid.iStart + 1)
			//{
			//	vectorB(index) += -lambda* u(i - 1, j) / grid.dx2;
			//}
			//if (i == grid.iEnd - 1)
			//{
			//	vectorB(index) += -lambda* u(i + 1, j) / grid.dx2;
			//}
			//if (j == grid.jStart + 1)
			//{
			//	vectorB(index) += -lambda* u(i, j - 1) / grid.dy2;
			//}
			//if (j == grid.jEnd - 1)
			//{
			//	vectorB(index) += -lambda* u(i, j + 1) / grid.dy2;
			//}
		}
	}
}

inline void VortexSheet::Stream2Velocity()
{
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			velocityX(i, j) = -streamFunction.dyPhi(i, j);
			velocityY(i, j) = streamFunction.dxPhi(i, j);
		}
	}
}


inline double VortexSheet::AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2)
{
	double maxVel1 = 0;
	double maxVel2 = 0;

	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			if (abs(velocity1(i, j)) > maxVel1)
			{
				maxVel1 = abs(velocity1(i, j));
			}
			if (abs(velocity2(i, j)) > maxVel2)
			{
				maxVel2 = abs(velocity2(i, j));
			}
		}
	}
	return cflCondition*(grid.dx / maxVel1 + grid.dy / maxVel2);
}



inline double VortexSheet::DeltaFt(const double & ip)
{
	double eps =min(grid.dx,grid.dy);
	return (1 + cos(PI*ip / eps)) / (2 * eps);
}
