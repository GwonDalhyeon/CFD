#pragma once


#include "AdvectionMethod2D.h"


class EulerianFluidSolver2D
{
public:
	Grid2D gridP;
	Grid2D gridPinner;
	Grid2D gridU;
	Grid2D gridUinner;
	Grid2D gridV;
	Grid2D gridVinner;
	

	Field2D<double> U; // x velocity
	Field2D<double> Ustar;
	Field2D<double> V; // y velocity
	Field2D<double> Vstar;
	Field2D<double> P; // Pressure
	Field2D<double> Rho;
	Field2D<double> Mu;

	LevelSet2D levelSet;


	double reynoldNum;
	double dt;
	
	double cflCondition;

	int ghostWidth;

	int maxIteration;
	int writeOutputIteration;

	EulerianFluidSolver2D();
	~EulerianFluidSolver2D();

	inline void InitialCondition(const int& example);
	inline void FluidSolver(const int& example);

	inline void GenerateLinearSystem(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystem(VectorND<double>& vectorB, const double & scaling);
	inline void TVDRK3TimeAdvection();
	inline void AdvectionTerm(const Field2D<double>& U, const Field2D<double>& V, Field2D<double>& TermU, Field2D<double>& TermV);
	inline void DiffusionTerm(const Field2D<double>& U, const Field2D<double>& V, Field2D<double>& TermU, Field2D<double>& TermV);

	inline double AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2);

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

		int numP = 11;
		double ddx = 0.01;
		gridP = Grid2D(0, ddx*double(numP-1), numP, 0, ddx*double(numP - 1), numP);
		gridPinner = Grid2D(gridP.xMin + gridP.dx, gridP.xMax - gridP.dx, 1, gridP.iRes - 2,
			gridP.yMin + gridP.dy, gridP.yMax - gridP.dy, 1, gridP.jRes - 2);
		gridU = Grid2D(gridP.xMin - gridP.dx/2, gridP.xMax + gridP.dx/2, gridP.iRes + 1,
			gridP.yMin, gridP.yMax, gridP.jRes);
		gridUinner = Grid2D(gridU.xMin + gridU.dx, gridU.xMax - gridU.dx, 1, gridU.iRes - 2,
			gridU.yMin + gridU.dy, gridU.yMax - gridU.dy, 1, gridU.jRes - 2);
		gridV = Grid2D(gridP.xMin, gridP.xMax, gridP.iRes, 
			gridP.yMin - gridP.dy / 2, gridP.yMax + gridP.dy / 2, gridP.jRes + 1);
		gridVinner = Grid2D(gridV.xMin + gridV.dx, gridV.xMax - gridV.dx, 1, gridV.iRes - 2,
			gridV.yMin + gridV.dy, gridV.yMax - gridV.dy, 1, gridV.jRes - 2);
		P = Field2D<double>(gridP);
		U = Field2D<double>(gridU);
		V = Field2D<double>(gridV);

		Rho = Field2D<double>(gridP);

		// initial condition
#pragma omp parallel for
		for (int i = gridU.iStart; i <= gridU.iEnd; i++)
		{
			U(i, gridU.jEnd) = 1;
		}
		Ustar = U;
		Vstar = V;

#pragma omp parallel for
		for (int i = Rho.iStart; i <= Rho.iEnd; i++)
		{
			for (int j = Rho.jStart; j <= Rho.jEnd; j++)
			{
				Rho(i, j) = 1;
			}
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

	Array2D<double> poissonMatrix(1, gridPinner.iRes*gridPinner.jRes, 1, gridPinner.iRes*gridPinner.jRes);
	VectorND<double> vectorB(gridPinner.iRes*gridPinner.jRes);
	VectorND<double> tempP(gridPinner.iRes*gridPinner.jRes);

	GenerateLinearSystem(poissonMatrix, -gridP.dx*gridP.dx);
	poissonMatrix.Variable("poisson");

	int solver = 2;
	// CG solver 1
	CSR<double> poissonCSR(poissonMatrix);

	// CG solver 2
	VectorND<double> a;
	VectorND<int> row;
	VectorND<int> col;
	int nonzeroNum;
	CGSolver::SparseA(poissonMatrix, a, row, col, nonzeroNum);

	for (int i = 0; i < 1; i++)
	{
		////////////////////////////////////////////
		//     Projection Method 1 : advection    //
		////////////////////////////////////////////
		TVDRK3TimeAdvection();
		U.Variable("U");
		V.Variable("V");
		////////////////////////////////////////////
		//     Projection Method 2 : Poisson Eq   //
		////////////////////////////////////////////
		GenerateLinearSystem(vectorB, -gridP.dx*gridP.dx);
		vectorB.Variable("vectorB");
		if (solver == 1)
		{
			tempP = CGSolver::SolverCSR(poissonCSR, vectorB, gridP.dx*gridP.dy);

		}
		else if (solver == 2)
		{
			CGSolver::SolverSparse(poissonMatrix.iRes, a, row, col, vectorB, tempP);
		}
		tempP.Variable("tempP");

		int index;
#pragma omp parallel for private(index)
		for (int i = gridPinner.iStart; i <= gridPinner.iEnd; i++)
		{
			for (int j = gridPinner.jStart; j <= gridPinner.jEnd; j++)
			{
				index = (i - gridPinner.iStart) + (j - gridPinner.jStart)*gridPinner.iRes;
				P(i, j) = tempP(index);
			}
		}
#pragma omp parallel for
		for (int i = gridP.iStart; i <= gridP.iEnd; i++)
		{
			P(i, P.jStart) = P(i, P.jStart + 1);
			P(i, P.jEnd) = P(i, P.jEnd - 1);
		}
#pragma omp parallel for
		for (int j = gridP.jStart; j <= gridP.jEnd; j++)
		{
			P(P.iStart, j) = P(P.iStart + 1, j);
			P(P.iEnd, j) = P(P.iEnd - 1, j);
		}
		P.Variable("P");
		//////////////////////////////////////////
		//     Projection Method 3 : New U,V    //
		//////////////////////////////////////////


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

inline void EulerianFluidSolver2D::GenerateLinearSystem(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int innerIStart = gridPinner.iStart;
	int innerIEnd = gridPinner.iEnd;
	int innerJStart = gridPinner.jStart;
	int innerJEnd = gridPinner.jEnd;
	int innerIRes = gridPinner.iRes;
	int innerJRes = gridPinner.jRes;

	int index, leftIndex, rightIndex, bottomIndex, topIndex;
#pragma omp parallel for private(index, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int j = innerJStart; j <= innerJEnd; j++)
	{
		for (int i = innerIStart; i <= innerIEnd; i++)
		{
			index = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			leftIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart - 1)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			rightIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart + 1)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			bottomIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart - 1)*innerIRes;
			topIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart + 1)*innerIRes;
			// Boundary condition.
			if (j == innerJStart)
			{
				if (i == innerIStart)
				{
					//leftIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIEnd - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2)+0.5;
					//matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					//rightIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIStart - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					//matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				
			}
			else if (j>innerJStart && j<innerJEnd)
			{
				if (i == innerIStart)
				{
					//leftIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIEnd - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					//matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;

				}
				else if (i == innerIEnd)
				{
					//rightIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIStart - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					//matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;

				}
				else
				{
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				
			}
			else if (j == innerJEnd)
			{
				if (i == innerIStart)
				{
					//leftIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIEnd - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					//matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					//rightIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIStart - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					//matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				
			}
		}
	}
}

inline void EulerianFluidSolver2D::GenerateLinearSystem(VectorND<double>& vectorB, const double & scaling)
{
	int innerIStart = gridPinner.iStart;
	int innerIEnd = gridPinner.iEnd;
	int innerJStart = gridPinner.jStart;
	int innerJEnd = gridPinner.jEnd;
	int innerIRes = gridPinner.iRes;
	int innerJRes = gridPinner.jRes;

	int index;

#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;

			vectorB(index) = Rho(i, j) / dt*((U(i + 1, j) - U(i, j))*gridU.oneOverdx + (V(i, j + 1) - V(i, j))*gridV.oneOverdy);
			
			if (i == innerIStart)
			{
				vectorB(index) += -P(i - 1, j)*P.oneOverdx2;
			}
			if (i == innerIEnd)
			{
				vectorB(index) += -P(i + 1, j)*P.oneOverdx2;
			}
			if (j == innerJStart)
			{
				vectorB(index) += -P(i, j - 1)*P.oneOverdy2;

			}
			if (j == innerJEnd)
			{
				vectorB(index) += -P(i, j + 1)*P.oneOverdy2;

			}
			vectorB(index) *= scaling;
		}
	}
}

inline void EulerianFluidSolver2D::TVDRK3TimeAdvection()
{
	Field2D<double> originU = U;
	Field2D<double> originV = V;

	Field2D<double> K1U(gridUinner);
	Field2D<double> K1V(gridVinner);
	Field2D<double> K2U(gridUinner);
	Field2D<double> K2V(gridVinner);
	Field2D<double> K3U(gridUinner);
	Field2D<double> K3V(gridVinner);

	Field2D<double> advectionU(gridUinner);
	Field2D<double> advectionV(gridVinner);

	Field2D<double> diffusionU(gridUinner);
	Field2D<double> diffusionV(gridVinner);

	dt = AdaptiveTimeStep(U, V);

	AdvectionTerm(U, V, advectionU, advectionV);
	DiffusionTerm(U, V, diffusionU, diffusionV);

#pragma omp parallel for
	for (int i = K1U.iStart; i <= K1U.iEnd; i++)
	{
		for (int j = K1U.jStart; j <= K1U.jEnd; j++)
		{
			K1U(i, j) = dt*(-advectionU(i, j) + 1 / reynoldNum*diffusionU(i, j));
			U(i, j) = originU(i, j) + K1U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = K1V.iStart; i <= K1V.iEnd; i++)
	{
		for (int j = K1V.jStart; j <= K1V.jEnd; j++)
		{
			K1V(i, j) = dt*(-advectionV(i, j) + 1 / reynoldNum*diffusionV(i, j));
			V(i, j) = originV(i, j) + K1V(i, j);
		}
	}

	AdvectionTerm(U, V, advectionU, advectionV);
	DiffusionTerm(U, V, diffusionU, diffusionV);

#pragma omp parallel for
	for (int i = K1U.iStart; i <= K1U.iEnd; i++)
	{
		for (int j = K1U.jStart; j <= K1U.jEnd; j++)
		{
			K2U(i, j) = dt*(-advectionU(i, j) + 1 / reynoldNum*diffusionU(i, j));
			U(i, j) = 3 / 4 * originU(i, j) + 1 / 4 * U(i, j) + 1 / 4 * K2U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = K1V.iStart; i <= K1V.iEnd; i++)
	{
		for (int j = K1V.jStart; j <= K1V.jEnd; j++)
		{
			K2V(i, j) = dt*(-advectionV(i, j) + 1 / reynoldNum*diffusionV(i, j));
			V(i, j) = 3 / 4 * originV(i, j) + 1 / 4 * V(i, j) + 1 / 4 * K2V(i, j);
		}
	}

	AdvectionTerm(U, V, advectionU, advectionV);
	DiffusionTerm(U, V, diffusionU, diffusionV);

#pragma omp parallel for
	for (int i = K1U.iStart; i <= K1U.iEnd; i++)
	{
		for (int j = K1U.jStart; j <= K1U.jEnd; j++)
		{
			K3U(i, j) = dt*(-advectionU(i, j) + 1 / reynoldNum*diffusionU(i, j));
			U(i, j) = 1 / 3 * originU(i, j) + 2 / 3 * U(i, j) + 2 / 3 * K3U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = K1V.iStart; i <= K1V.iEnd; i++)
	{
		for (int j = K1V.jStart; j <= K1V.jEnd; j++)
		{
			K2V(i, j) = dt*(-advectionV(i, j) + 1 / reynoldNum*diffusionV(i, j));
			V(i, j) = 1 / 3 * originV(i, j) + 2 / 3 * V(i, j) + 2 / 3 * K3V(i, j);
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


	double Ux, Uy;
	double aveV;
#pragma omp parallel for private(aveV, Ux, Uy)
	for (int i = TermU.iStart; i <= TermU.iEnd; i++)
	{
		for (int j = TermU.jStart; j <= TermU.jEnd; j++)
		{
			aveV = (V(i - 1, j) + V(i, j) + V(i - 1, j + 1) + V(i, j + 1)) / 4;
			if (U(i,j)>0)
			{
				Ux = dUdxM(i, j);
			}
			else
			{
				Ux = dUdxP(i, j);
			}

			if (aveV>0)
			{
				Uy = dUdyM(i, j);
			}
			else
			{
				Uy = dUdyP(i, j);
			}
			TermU(i, j) = U(i, j)*Ux + aveV*Uy;
		}
	}

	double aveU;
	double Vx, Vy;
#pragma omp parallel for private(aveU, Vx, Vy)
	for (int i = TermV.iStart; i <= TermV.iEnd; i++)
	{
		for (int j = TermV.jStart; j <= TermV.jEnd; j++)
		{
			aveU = (U(i, j - 1) + U(i, j) + U(i, j) + U(i + 1, j)) / 4;
			if (aveU>0)
			{
				Vx = dVdxM(i, j);
			}
			else
			{
				Vx = dVdxP(i, j);
			}

			if (V(i,j)>0)
			{
				Vy = dVdyM(i, j);
			}
			else
			{
				Vy = dVdyP(i, j);
			}
			TermV(i, j) = aveU*Vx + V(i, j)*Vy;
		}
	}

}

inline void EulerianFluidSolver2D::DiffusionTerm(const Field2D<double>& U, const Field2D<double>& V, Field2D<double>& TermU, Field2D<double>& TermV)
{
	Field2D<double> tempV(TermU.grid);
#pragma omp parallel for
	for (int i = TermU.iStart; i <= TermU.iEnd; i++)
	{
		for (int j = TermU.jStart; j <= TermU.jEnd; j++)
		{
			tempV(i, j) = (V(i - 1, j) + V(i, j) + V(i - 1, j + 1) + V(i, j + 1)) / 4;
		}
	}

#pragma omp parallel for
	for (int i = TermU.iStart; i <= TermU.iEnd; i++)
	{
		for (int j = TermU.jStart; j <= TermU.jEnd; j++)
		{
			TermU(i, j) = (U(i - 1, j) - 2 * U(i, j) + U(i + 1, j))*U.oneOver2dx;
			TermU(i, j) += (tempV(i, j - 1) - 2 * tempV(i, j) + tempV(i, j + 1))*U.oneOver2dy;
		}
	}

	Field2D<double> tempU(TermV.grid);
#pragma omp parallel for 
	for (int i = TermV.iStart; i <= TermV.iEnd; i++)
	{
		for (int j = TermV.jStart; j <= TermV.jEnd; j++)
		{
			tempU(i, j) = (U(i, j - 1) + U(i + 1, j - 1) + U(i, j) + U(i + 1, j)) / 4;
		}
	}

#pragma omp parallel for 
	for (int i = TermV.iStart; i <= TermV.iEnd; i++)
	{
		for (int j = TermV.jStart; j <= TermV.jEnd; j++)
		{
			TermV(i, j) = (tempU(i - 1, j) - 2 * tempU(i, j) + tempU(i + 1, j))*V.oneOver2dx;
			TermV(i, j) += (V(i, j - 1) - 2 * V(i, j) + V(i, j + 1))*V.oneOver2dy;
		}
	}
}

inline double EulerianFluidSolver2D::AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2)
{
	double maxVel = 0;
	for (int i = velocity1.iStart; i <= velocity1.iEnd; i++)
	{
		for (int j = velocity1.jStart; j <= velocity1.jEnd; j++)
		{
			if (abs(velocity1(i,j))>maxVel)
			{
				maxVel = abs(velocity1(i, j));
			}
		}
	}
	for (int i = velocity2.iStart; i <= velocity2.iEnd; i++)
	{
		for (int j = velocity2.jStart; j <= velocity2.jEnd; j++)
		{
			if (abs(velocity2(i, j))>maxVel)
			{
				maxVel = abs(velocity2(i, j));
			}
		}
	}
	return cflCondition*min(gridU.dx, gridV.dy) / maxVel;
}

