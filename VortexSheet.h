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

	int innerIStart;
	int innerJStart;
	int innerIEnd;
	int innerJEnd;
	int innerIRes;
	int innerJRes;

	int ghostWidth;

	int maxIteration;
	int writeOutputIteration;

	VortexSheet();
	~VortexSheet();

	inline void InitialCondition(const int& example);
	inline void VortexSolver(const int& example);

	inline void GenerateLinearSystem(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystem(const Field2D<double>& P, VectorND<double>& vectorB, const double & scaling);
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
		cout << "*************************" << endl;
		cout << "    Vortex Sheet in 2D" << endl;
		cout << "*************************" << endl;
		grid = Grid2D(-1, 1, 101, -1, 1, 101);
		levelSet = LevelSet2D(grid);
		P = Field2D<double>(grid);
		streamFunction = Field2D<double>(grid);
		velocityX = Field2D<double>(grid);
		velocityY = Field2D<double>(grid);

		cflCondition = 0.8;

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
		cout << "*************************" << endl;
		cout << "    Vortex Sheet Dipole" << endl;
		cout << "*************************" << endl;
		grid = Grid2D(-1, 1, 101, -1, 1, 101);
		levelSet = LevelSet2D(grid);
		P = Field2D<double>(grid);
		streamFunction = Field2D<double>(grid);
		velocityX = Field2D<double>(grid);
		velocityY = Field2D<double>(grid);
		
		cflCondition = 0.8;
		
		maxIteration = 1000;
		writeOutputIteration = 1;

		double eps = 8*grid.dx;

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

inline void VortexSheet::VortexSolver(const int & example)
{
	bool writeFile = false;
	string str;
	const char*cmd;
	double totalT = 0;
	

	InitialCondition(example);
	grid.Variable();
	//P.Variable("P");
	
	innerIStart = grid.iStart;
	innerJStart = grid.jStart + 1;
	innerIEnd = grid.iEnd - 1;
	innerJEnd = grid.jEnd - 1;
	innerIRes = innerIEnd - innerIStart + 1;
	innerJRes = innerJEnd - innerJStart + 1;
	
	Array2D<double> poissonMatrix(1, innerIRes*innerJRes, 1, innerIRes*innerJRes);

	GenerateLinearSystem(poissonMatrix, -grid.dx*grid.dy);

	CSR<double> poissonCSR(poissonMatrix);

	VectorND<double> a;
	VectorND<int> row;
	VectorND<int> col;
	int nonzeroNum;
	CGSolver::SparseA(poissonMatrix, a, row, col, nonzeroNum);

	VectorND<double> vectorB(innerIRes*innerJRes);

	VectorND<double> streamV(innerIRes*innerJRes);

	double eps = 8*grid.dx;
	int idx;	


	if (writeFile)
	{
		str = "phi0";
		levelSet.phi.WriteFile(str);

	}

	//// Write Movie 1-3
	//MATLAB.Command("writerobj = VideoWriter('WhereIsFile.avi');writerobj.FrameRate = 10;open(writerobj); ");
	//MATLAB.Command("fig = figure('units','normalized','outerposition',[0 0 1 1])");
	MATLAB.Variable("eps", eps);
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	levelSet.phi.Variable("phi0");
	MATLAB.Command("subplot(1, 3, 1)");
	MATLAB.Command("surf(X,Y,phi0)");
	MATLAB.Command("subplot(1, 3, 2)");
	if (example == 1)
	{
		MATLAB.Command("contour(X, Y, phi0, [0 0],'b');");
		MATLAB.Command("grid on");
		str = string("title(['iteration : ', num2str(") + to_string(0) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
	}
	else if (example == 2)
	{
		str = string("contour(X, Y, phi0, [") + to_string(-eps / 2) + string(",") + to_string(eps / 2) + string("],'b');");
		cmd = str.c_str();
		MATLAB.Command(cmd);
		MATLAB.Command("grid on");
		str = string("title(['iteration : ', num2str(") + to_string(0) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
	}

	MATLAB.Command("subplot(1, 3, 3)");
	velocityX.Variable("velocityX");
	velocityY.Variable("velocityY");
	MATLAB.Command("quiver(X,Y,velocityX,velocityY);");
	
	/////////////////////////
	////                /////
	////    Iteration   /////
	////                /////
	/////////////////////////
	//for (int i = 1; i <= maxIteration; i++)
	for (int i = 1; i <= 200; i++)
	{
		cout << endl;
		cout << "*************************************************" << endl;
		cout << "iteration : " << i << endl;

		GenerateLinearSystem(P, vectorB, -grid.dx*grid.dy);
		int solver = 1;
		if (solver == 1)
		{
			streamV = CGSolver::SolverCSR(poissonCSR, vectorB, grid.dx*grid.dy);
		}
		else if (solver == 2)
		{
			CGSolver::SolverSparse(poissonMatrix.iRes, a, row, col, vectorB, streamV);
		}
		streamV.Variable("streamV");

#pragma omp parallel for private(idx)
		for (int i = innerIStart; i <= innerIEnd; i++)
		{
			streamFunction(i, grid.jStart) = 0;
			streamFunction(i, grid.jEnd) = 0;
			for (int j = innerJStart; j <= innerJEnd; j++)
			{
				idx = (i - innerIStart) + (j - innerJStart)*innerIRes;
				streamFunction(i, j) = streamV(idx);

				if (i==innerIStart)
				{
					streamFunction(grid.iEnd, j) = streamFunction(grid.iStart, j);
				}
			}
		}
		
		Stream2Velocity();

		streamFunction.Variable("stream");
		velocityX.Variable("velocityX");
		velocityY.Variable("velocityY");
		
		dt = AdaptiveTimeStep(velocityX, velocityY);
		if (dt==INFINITE || dt==INFINITY ||dt==NAN)
		{
			cout << "************************************" << endl;
			cout << "             dt error!!"<<endl;
			cout << "************************************" << endl;
		}
		totalT += dt;

		// Left and Right side
		// Level set boundary condition
#pragma omp parallel for
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			if (velocityX(grid.iStart,j)>0)
			{
				levelSet(grid.iStart, j) = levelSet(grid.iEnd, j);
			}
			else
			{
				levelSet(grid.iEnd, j) = levelSet(grid.iStart, j);
			}
		}


		AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, velocityX, velocityY, dt);
		levelSet.phi.Variable("phi");

#pragma omp parallel for 
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (example == 1)
				{
					P(i, j) = DeltaFt(levelSet(i, j));
				}
				else if (example == 2)
				{
					if (abs(levelSet(i, j))<8 * grid.dx)
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
		P.Variable("P");

		MATLAB.Command("subplot(1, 3, 1)");
		MATLAB.Command("surf(X,Y,phi)");
		if (example == 1)
		{
			MATLAB.Command("subplot(1, 3, 2)");
			MATLAB.Command("contour(X, Y, phi0, [0 0],'b');");
			MATLAB.Command("hold on");
			MATLAB.Command("contour(X, Y, phi, [0 0],'r');");
			MATLAB.Command("grid on");
			MATLAB.Command("hold off");
			str = string("title(['iteration : ', num2str(") + to_string(i) + string(")]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
		}
		else if (example == 2)
		{
			MATLAB.Command("subplot(1, 3, 2)");
			str = string("contour(X, Y, phi0, [") + to_string(-eps / 2) + string(",") + to_string(eps / 2) + string("],'b');");
			cmd = str.c_str();
			MATLAB.Command(cmd);
			MATLAB.Command("hold on");
			str = string("contour(X, Y, phi, [") + to_string(-eps / 2) + string(",") + to_string(eps / 2) + string("],'r');");
			cmd = str.c_str();
			MATLAB.Command(cmd);
			MATLAB.Command("grid on");
			MATLAB.Command("hold off");
			str = string("title(['iteration : ', num2str(") + to_string(i) + string(")]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
		}
		MATLAB.Command("subplot(1, 3, 3)");
		MATLAB.Command("quiver(X,Y,velocityX,velocityY);");



		if (writeFile && i%writeOutputIteration==0)
		{
			str = "velocityX" + to_string(i);
			velocityX.WriteFile(str);

			str = "velocityY" + to_string(i);
			velocityY.WriteFile(str);

			str = "phi" + to_string(i);
			levelSet.phi.WriteFile(str);

			str = "stream" + to_string(i);
			P.WriteFile(str);
		}

		//// Write Movie 2-3
		//MATLAB.Command("f = getframe(fig);writeVideo(writerobj,f);");
	}
	//// Write Movie 3-3
	//MATLAB.Command("close(fig);close(writerobj);");
}



inline void VortexSheet::GenerateLinearSystem(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int index, leftIndex, rightIndex, bottomIndex, topIndex;

#pragma omp parallel for private(index, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int j = innerJStart; j <= innerJEnd; j++)
	{
		for (int i = innerIStart; i <= innerIEnd; i++)
		{
			index = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			leftIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart - 1) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			rightIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart + 1) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			bottomIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart - 1)*innerIRes;
			topIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart + 1)*innerIRes;


			// Boundary condition.
			if (j==innerJStart)
			{
				if (i == innerIStart)
				{
					leftIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIEnd - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
				}
				else if (i == innerIEnd)
				{
					rightIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIStart - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
				}
				matrixA(index) = scaling*(-2 * grid.oneOverdx2 - 2 * grid.oneOverdy2);
				//matrixA(index) = scaling*(-2 * grid.oneOverdx2);
				matrixA(leftIndex) = scaling * 1 * grid.oneOverdx2;
				matrixA(rightIndex) = scaling * 1 * grid.oneOverdx2;
				matrixA(topIndex) = scaling * 1 * grid.oneOverdy2;
			}
			else if (j>innerJStart && j<innerJEnd)
			{
				if (i == innerIStart)
				{
					leftIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIEnd - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
				}
				else if (i == innerIEnd)
				{
					rightIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIStart - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
				}
				matrixA(index) = scaling*(-2 * grid.oneOverdx2 - 2 * grid.oneOverdy2);
				matrixA(leftIndex) = scaling * 1 * grid.oneOverdx2;
				matrixA(rightIndex) = scaling * 1 * grid.oneOverdx2;
				matrixA(bottomIndex) = scaling * 1 * grid.oneOverdy2;
				matrixA(topIndex) = scaling * 1 * grid.oneOverdy2;
			}
			else if (j==innerJEnd)
			{
				if (i == innerIStart)
				{
					leftIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIEnd - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
				}
				else if (i == innerIEnd)
				{
					rightIndex = (i - innerIStart)*innerIRes*innerJRes + (innerIStart - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
				}
				matrixA(index) = scaling*(-2 * grid.oneOverdx2 - 2 * grid.oneOverdy2);
				matrixA(leftIndex) = scaling * 1 * grid.oneOverdx2;
				matrixA(rightIndex) = scaling * 1 * grid.oneOverdx2;
				matrixA(bottomIndex) = scaling * 1 * grid.oneOverdy2;
			}
	
		}
	}
	cout << "End Generate Linear System : matrix A" << endl;
	//matrixA.Variable("A");
}


inline void VortexSheet::GenerateLinearSystem(const Field2D<double>& P, VectorND<double>& vectorB, const double & scaling)
{
	int index;

#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;

			vectorB(index) = scaling*(-P(i, j));
		}
	}
	//vectorB.Variable("B");
}

inline void VortexSheet::Stream2Velocity()
{
	Field2D<double> wenoXMinus(grid);
	Field2D<double> wenoXPlus(grid);
	Field2D<double> wenoYMinus(grid);
	Field2D<double> wenoYPlus(grid);
	AdvectionMethod2D<double>::WENO5th(streamFunction, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);

#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			if (abs(wenoYMinus(i, j))<abs(wenoYPlus(i, j)))
			{
				velocityX(i, j) = -wenoYMinus(i, j);
			}
			else
			{
				velocityX(i, j) = -wenoYPlus(i, j);
			}

			if (i == grid.iStart)
			{
				velocityY(i, j) = wenoXPlus(i, j);
				//velocityY(i, j) = (streamFunction(grid.iEnd - 1, j) - streamFunction(grid.iStart + 1, j))*grid.oneOver2dx;
			}
			else if (i == grid.iEnd)
			{
				velocityY(i, j) = wenoXMinus(i, j);
				//velocityY(i, j) = (streamFunction(grid.iEnd - 1, j) - streamFunction(grid.iStart + 1, j))*grid.oneOver2dx;
			}
			else
			{
				if (abs(wenoXMinus(i, j))<abs(wenoXPlus(i, j)))
				{
					velocityY(i, j) = wenoXMinus(i, j);
				}
				else
				{
					velocityY(i, j) = wenoXPlus(i, j);
				}
			}
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
	return cflCondition*(grid.dx / max(maxVel1, maxVel2) + grid.dy / max(maxVel1, maxVel2));
}



inline double VortexSheet::DeltaFt(const double & ip)
{
	double eps = 3 * min(grid.dx, grid.dy);
	if (abs(ip)<eps)
	{
		return (1 + cos(PI*ip / eps)) / (2 * eps);
	}
	else
	{
		return 0.0;
	}
}
