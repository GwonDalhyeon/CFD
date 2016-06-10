#pragma once


#include "AdvectionMethod2D.h"


class EulerianFluidSolver2D
{
public:
	Grid2D gridU;
	Grid2D gridV;
	Grid2D gridP;

	Field2D<double> U; // x velocity
	Field2D<double> V; // y velocity
	Field2D<double> P; // Pressure

	LevelSet2D levelSet;


	double reynoldNum;
	double dt;
	double cflCondition;

	int ghostWidth;

	int maxIteration;
	int writeOutputIteration;

	EulerianFluidSolver2D();
	~EulerianFluidSolver2D();

	void InitialCondition(const int& example);
	void FluidSolver(const int& example);

	void TimeAdvanceForwardEuler(Field2D<Vector2D<double>>& ipVelocity);
	void AdvectionTerm(const Field2D<double>& U, const Field2D<double>& V, Field2D<double>& TermU, Field2D<double>& TermV);
	
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
		cout << "*************************" << endl;
		cout << "    Cavity Flow" << endl;
		cout << "*************************" << endl;

		gridP = Grid2D(0, 1.0, 101, 0, 1, 101);
		gridU = Grid2D(gridP.xMin - gridP.dx/2, gridP.xMax + gridP.dx/2, gridP.iRes + 1,
			gridP.yMin, gridP.yMax, gridP.jRes);
		gridV = Grid2D(gridP.xMin, gridP.xMax, gridP.iRes, 
			gridP.yMin - gridP.dy / 2, gridP.yMax + gridP.dy / 2, gridP.jRes + 1);

		P = Field2D<double>(gridP);
		U = Field2D<double>(gridU);
		V = Field2D<double>(gridV);

		// initial condition
#pragma omp parallel for
		for (int i = gridU.iStart; i <= gridU.iEnd; i++)
		{
			U(i, gridU.jEnd) = 1;
		}

		reynoldNum = 500;
		cflCondition = 0.5;
		
		maxIteration = 1000;
		writeOutputIteration = 10;


	}

	if (example==2)
	{

	}
}

inline void EulerianFluidSolver2D::FluidSolver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	const char* cmd;

	InitialCondition(example);
	gridP.Variable("Xp", "Yp");
	gridU.Variable("Xu", "Yu");
	gridV.Variable("Xv", "Yv");


	
	for (int i = 0; i < 0; i++)
	{

		if (writeFile && i%writeOutputIteration == 0)
		{
			fileName = "pressure" + to_string(i);
			P.WriteFile(fileName);
			fileName = "xVelocity" + to_string(i);
			U.WriteFile(fileName);
			fileName = "yVelocity" + to_string(i);
			V.WriteFile(fileName);
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

inline void EulerianFluidSolver2D::AdvectionTerm(const Field2D<double>& U, const Field2D<double>& V, Field2D<double>& TermU, Field2D<double>& TermV)
{
	Field2D<double> dUdxM(U.grid);
	Field2D<double> dUdxP(U.grid);
	Field2D<double> dUdyM(U.grid);
	Field2D<double> dUdyP(U.grid);
	AdvectionMethod2D<double>::ENO3rdDerivation(U, dUdxM, dUdxP, dUdyM, dUdyP);
	Field2D<double> dVdxM(V.grid);
	Field2D<double> dVdxP(V.grid);
	Field2D<double> dVdyM(V.grid);
	Field2D<double> dVdyP(V.grid);
	AdvectionMethod2D<double>::ENO3rdDerivation(V, dVdxM, dVdxP, dVdyM, dVdyP);

	int Uistart = U.iStart + 1;
	int Uiend = U.iEnd - 1;
	int Ujstart = U.jStart + 1;
	int Ujend = U.jEnd - 1;
	int Vistart = V.iStart + 1;
	int Viend = V.iEnd - 1;
	int Vjstart = V.jStart + 1;
	int Vjend = V.jEnd - 1;

	for (int i = Uistart; i <= Uiend; i++)
	{
		for (int j = Ujstart; j <= Ujend; j++)
		{

		}
	}
}

