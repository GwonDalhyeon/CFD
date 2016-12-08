#pragma once

#include "AdvectionMethod2D.h"
#include "EulerianFluidSolver.h"

class MovingInterface
{
public:
	MovingInterface(EulerianFluidSolver2D & ipFluid);
	~MovingInterface();

	int ExamNum;

	EulerianFluidSolver2D& Fluid;

	Grid2D& grid = Fluid.gridP;

	FD& U = Fluid.U; // x velocity
	FD& V = Fluid.V; // y velocity

	FD Surfactant;
	VectorND<double> tempSur;

	double initialSurfactant;


	FD SurfaceTension;
	double gamma0;
	double SurfactantMax;
	double IdealGasConstant;
	double AbsTemperature;

	LS levelSet;

	Array2D<double> term;
	Array2D<double> termOld;

	VectorND<double> vectorB;
	Array2D<double>A;
	// CG solver 1
	CSR<double> Acsr;


	double dt;
	double totalT;

	double cflCondition;

	int maxIteration;
	int writeOutputIteration;

	inline void InitialCondition(const int& example);

	////////////////////////////////////////////
	////   Surfactant Diffusion Solver      ////
	////////////////////////////////////////////
	inline void SurfactantDiffusionSolver(const int& example);

	inline void SurfactantDiffusion(const int& iter);
	inline void OneStepSemiImplicit();
	inline void TwoStepSemiImplicit();

	inline void SurfactantNormalTerm(FD& ipField, LS& ipLevelSet, Array2D<double>& term);
	inline double ExactSurfactant(const double& x, const double& y, const double& time);
	inline void GenerateLinearSystem1(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystem1(VectorND<double>& vectorB, const double & scaling);
	inline void GenerateLinearSystem2(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystem2(VectorND<double>& vectorB, const double & scaling);

	//////////////////////////////////////////////
	////   Local Surfactant Diffusion Solver  ////
	//////////////////////////////////////////////
	inline void LSurfactantDiffusionSolver(const int& example);

	///////////////////////////////////////////
	//// Eulerian Moving Interface Solver  ////
	///////////////////////////////////////////
	inline void EulerianMovingInterfaceSolver(const int& example);

	inline void SurfactantTube2Extrapolation();

	inline void LSurfactantDiffusion(const int& iter);
	inline void LOneStepSemiImplicit();
	inline void LTwoStepSemiImplicit();
	inline void LSurfactantNormalTerm(FD& ipField, LS& ipLevelSet, Array2D<double>& term);
	inline void LGenerateLinearSystem1(Array2D<double>& matrixA, const double & scaling);
	inline void LGenerateLinearSystem1(VectorND<double>& vectorB, const double & scaling);
	inline void LGenerateLinearSystem2(Array2D<double>& matrixA, const double & scaling);
	inline void LGenerateLinearSystem2(VectorND<double>& vectorB, const double & scaling);

	// Compute Surface Tension
	inline void DimlessNonlinearLangmu1rEOS();
	inline void DimlessNonlinearLangmu1rEOS(const int & tubeRange);
	inline void DimlessLinearLangmu1rEOS();
	inline void DimlessLinearLangmu1rEOS(const int & tubeRange);
	inline void LinearLangmu1rEOS();

	// Conserve Surfactant
	inline double IntegralSurfactant();
	inline void ConserveSurfactantFactorBeta();

private:

};

inline MovingInterface::MovingInterface(EulerianFluidSolver2D & ipFluid)
	:Fluid(ipFluid)
{
}

MovingInterface::~MovingInterface()
{
}

inline void MovingInterface::InitialCondition(const int & example)
{
	ExamNum = example;
	if (example == 1)
	{
		cout << "*************************" << endl;
		cout << "       An Eulerian Formulation  " << endl;
		cout << " for Solving PDE along Moving Interface" << endl;
		cout << "          --JJ Xu, HK Zhao-- " << endl;
		cout << "               Example 1 " << endl;
		cout << "*************************" << endl;

		int gridSize = 40;
		grid = Grid2D(-2, 2, gridSize + 1, -2, 2, gridSize + 1);

		U = FD(grid);
		V = FD(grid);

		levelSet = LS(grid);
		double radius = 1.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}

		Surfactant = FD(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);

			}
		}

		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		int innerIStart = grid.iStart + 1;
		int innerIEnd = grid.iEnd - 1;
		int innerJStart = grid.jStart + 1;
		int innerJEnd = grid.jEnd - 1;
		int innerIRes = grid.iRes - 2;
		int innerJRes = grid.jRes - 2;
		term = Array2D<double>(innerIStart, innerIRes, innerJStart, innerJRes);
		termOld = Array2D<double>(innerIStart, innerIRes, innerJStart, innerJRes);
		tempSur = VectorND<double>((grid.iRes - 2)*(grid.jRes - 2));
		vectorB = VectorND<double>((grid.iRes - 2)*(grid.jRes - 2));
		A = Array2D<double>(1, (grid.iRes - 2)*(grid.jRes - 2), 1, (grid.iRes - 2)*(grid.jRes - 2));

		totalT = 0;
		dt = grid.dx / 4.0;
		maxIteration = 80;
	}

	if (example == 2)
	{
		cout << "*************************" << endl;
		cout << "       An Eulerian Formulation  " << endl;
		cout << " for Solving PDE along Moving Interface" << endl;
		cout << "          --JJ Xu, HK Zhao-- " << endl;
		cout << "               Example 2 " << endl;
		cout << "*************************" << endl;

		int gridSize = 40;
		grid = Grid2D(-3, 5, gridSize + 1, -3, 3, gridSize * 3. / 4. + 1);

		U = FD(grid);
		V = FD(grid);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				U(i, j) = 1;
				V(i, j) = 0;
			}
		}

		levelSet = LS(grid);
		double radius = 2.0;
		double x0, y0;
#pragma omp parallel for private(x0, y0)
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				x0 = grid(i, j).x;
				y0 = grid(i, j).y;
				levelSet(i, j) = sqrt(x0*x0 + y0*y0) - radius;
			}
		}

		Surfactant = FD(grid);
		int i, j;
#pragma omp parallel for private(x0, y0)
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				x0 = grid(i, j).x;
				y0 = grid(i, j).y;
				Surfactant(i, j) = ExactSurfactant(x0, y0, totalT);
			}
		}

		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		int innerIStart = grid.iStart + 1;
		int innerIEnd = grid.iEnd - 1;
		int innerJStart = grid.jStart + 1;
		int innerJEnd = grid.jEnd - 1;
		int innerIRes = grid.iRes - 2;
		int innerJRes = grid.jRes - 2;
		term = Array2D<double>(innerIStart, innerIRes, innerJStart, innerJRes);
		termOld = Array2D<double>(innerIStart, innerIRes, innerJStart, innerJRes);
		tempSur = VectorND<double>((grid.iRes - 2)*(grid.jRes - 2));
		vectorB = VectorND<double>((grid.iRes - 2)*(grid.jRes - 2));
		A = Array2D<double>(1, (grid.iRes - 2)*(grid.jRes - 2), 1, (grid.iRes - 2)*(grid.jRes - 2));
		cflCondition = 1.0 / 4.0;
		dt = AdvectionMethod2D<double>::AdaptiveTimeStep(U, cflCondition);
		maxIteration = 80;
		totalT = 0;
	}

	if (example == 3)
	{
		cout << "*************************" << endl;
		cout << "       An Eulerian Formulation  " << endl;
		cout << " for Solving PDE along Moving Interface" << endl;
		cout << "          --JJ Xu, HK Zhao-- " << endl;
		cout << "               Example 1-1 " << endl;
		cout << "*************************" << endl;

		int gridSize = 80;
		grid = Grid2D(-2, 2, gridSize + 1, -2, 2, gridSize + 1);

		U = FD(grid);
		V = FD(grid);

		levelSet = LS(grid);
		double radius = 1.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}
		levelSet.InitialTube();

		Surfactant = FD(grid);
		int i, j;
#pragma omp parallel for private(i, j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			levelSet.TubeIndex(k, i, j);
			if (levelSet.tube(i, j) <= 1)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}
		//MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
		//grid.Variable();
		//Surfactant.Variable("Surfactant");
		//levelSet.tube.Variable("Tube");
		//MATLAB.Command("subplot(1,2,1),surf(X,Y,Surfactant), hold on, contour(X,Y,Tube),hold off;");
		int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.1*min(grid.dx, grid.dy)));
		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);
		//Surfactant.Variable("Surfactant");
		//MATLAB.Command("subplot(1,2,2),surf(X,Y,Surfactant), hold on, contour(X,Y,Tube),hold off;");
		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);
		dt = grid.dx / 4.0;
		maxIteration = 160;
		totalT = 0;

	}

	if (example == 4)
	{
		cout << "*************************" << endl;
		cout << "       An Eulerian Formulation  " << endl;
		cout << " for Solving PDE along Moving Interface" << endl;
		cout << "          --JJ Xu, HK Zhao-- " << endl;
		cout << "               Example 4 " << endl;
		cout << "*************************" << endl;

		int gridSize = 250;
		grid = Grid2D(-2, 8, gridSize + 1, -2, 2, gridSize * 4. / 10. + 1);

		U = FD(grid);
		V = FD(grid);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				U(i, j) = 1;
				V(i, j) = 0;
			}
		}

		levelSet = LS(grid, 3 * grid.dx);
		double radius = 1.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}
		levelSet.InitialTube();

		Surfactant = FD(grid);
		int i, j;
#pragma omp parallel for private(i, j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			levelSet.TubeIndex(k, i, j);
			if (levelSet.tube(i, j) <= 1)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}

		int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.1*min(grid.dx, grid.dy)));
		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);
		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);
		cflCondition = 1.0 / 4.0;
		dt = AdvectionMethod2D<double>::AdaptiveTimeStep(U, cflCondition);
		maxIteration = 800;
		totalT = 0;
	}

	if (example == 5)
	{
		cout << "*************************" << endl;
		cout << "       An Eulerian Formulation  " << endl;
		cout << " for Solving PDE along Moving Interface" << endl;
		cout << "          --JJ Xu, HK Zhao-- " << endl;
		cout << "               Example 5 " << endl;
		cout << "*************************" << endl;

		int gridSize = 150;
		grid = Grid2D(-3, 3, gridSize + 1, -3, 3, gridSize + 1);

		U = FD(grid);
		V = FD(grid);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (grid(i, j).y >= 0)
				{
					U(i, j) = grid(i, j).y*grid(i, j).y;
				}
				else
				{
					U(i, j) = -grid(i, j).y*grid(i, j).y;
				}
				V(i, j) = 0;
			}
		}

		levelSet = LS(grid, 3 * grid.dx);
		double radius = 1.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}
		levelSet.InitialTube();

		Surfactant = FD(grid);
		int i, j;
#pragma omp parallel for private(i, j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			levelSet.TubeIndex(k, i, j);
			if (levelSet.tube(i, j) <= 1)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}

		int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.1*min(grid.dx, grid.dy)));
		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);
		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);
		cflCondition = 1.0 / 4.0;
		dt = cflCondition*min(grid.dx, grid.dy);
		maxIteration = 800;
		totalT = 0;
	}

	if (example == 6)
	{
		cout << "*************************" << endl;
		cout << "       An Eulerian Formulation  " << endl;
		cout << " for Solving PDE along Moving Interface" << endl;
		cout << "          --JJ Xu, HK Zhao-- " << endl;
		cout << "               Example 6 " << endl;
		cout << "*************************" << endl;

		int gridSize = 200;
		grid = Grid2D(-2, 6, gridSize + 1, -2, 2, gridSize * 1. / 2. + 1);

		U = FD(grid);
		V = FD(grid);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				U(i, j) = (grid(i, j).y + 2)*(grid(i, j).y + 2) / 3.0;
				V(i, j) = 0;
			}
		}

		levelSet = LS(grid, 3 * grid.dx);
		double radius = 1.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}
		levelSet.InitialTube();

		Surfactant = FD(grid);
		int i, j;
#pragma omp parallel for private(i, j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			levelSet.TubeIndex(k, i, j);
			if (levelSet.tube(i, j) <= 1)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}

		int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.1*min(grid.dx, grid.dy)));
		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);
		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);
		cflCondition = 1.0 / 4.0;
		dt = cflCondition*min(grid.dx, grid.dy);
		maxIteration = 800;
		totalT = 0;
	}

	if (example == 7)
	{
		/////////////////////////////////////////////////////////
		/////  A level-set continuum method
		/////  for two-phase flows with insoluble surfactant
		/////  --JJ Xu, Y Yang, J Lowengrub--
		/////  Example 1
		/////////////////////////////////////////////////////////
		levelSet = LS(grid, 3 * grid.dx);
		double radius = 1.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}
		levelSet.InitialTube();

		Surfactant = FD(grid);
		int i, j;
#pragma omp parallel for private(i, j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			levelSet.TubeIndex(k, i, j);
			if (levelSet.tube(i, j) <= 2)
			{
				Surfactant(i, j) = 1;
			}
		}

		SurfaceTension = FD(grid);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (levelSet.tube(i, j) <= 2)
				{
					SurfaceTension(i, j) = 1 + Fluid.El*log(1 - Fluid.Xi* Surfactant(i, j));
				}
				else
				{
					SurfaceTension(i, j) = 1;
				}
			}
		}

		initialSurfactant = IntegralSurfactant();

		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);
		cflCondition = 1.0 / 4.0;
		dt = cflCondition*min(grid.dx, grid.dy);
		totalT = 0;
	}

	if (example == 8)
	{
		/////////////////////////////////////////////////////////
		/////  Simulations of surfactant effects
		/////  on the dynamics of coalescing drops and bubbles
		/////  --DW Martin and F Blanchette--
		/////  Example 1
		/////////////////////////////////////////////////////////

		levelSet.InitialTube();

		Surfactant = FD(grid);
		int i, j;
#pragma omp parallel for private(i, j)
		for (int k = 1; k <= levelSet.numTube; k++)
		{
			levelSet.TubeIndex(k, i, j);
			if (levelSet.tube(i, j) <= 2)
			{
				Surfactant(i, j) = 1;
			}
		}

		SurfaceTension = FD(grid);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (levelSet.tube(i, j) <= 2)
				{
					SurfaceTension(i, j) = 1 + Fluid.El*log(1 - Fluid.Xi* Surfactant(i, j));
				}
				else
				{
					SurfaceTension(i, j) = 1;
				}
			}
		}

		initialSurfactant = IntegralSurfactant();

		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);
		totalT = 0;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////
////
////       Surfactant Diffusion Solver
////
/////////////////////////////////////////////////////////////////////////////////////////
inline void MovingInterface::SurfactantDiffusionSolver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	const char* cmd;


	InitialCondition(example);

	grid.Variable();

	//// Initial Surfactant
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	Surfactant.Variable("Surfactant");
	MATLAB.Command("surf(X,Y,Surfactant);");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(totalT) + string(")]), ;");
	cmd = str.c_str();
	MATLAB.Command(cmd);

	FD exact(grid);

	for (int i = 1; i <= maxIteration; i++)
	{
		cout << "********************************" << endl;
		cout << "       Iteration " << to_string(i) << " : Start" << endl;

		totalT += dt;
		SurfactantDiffusion(i);


		cout << "       Iteration " << to_string(i) << " : End" << endl;
		cout << "********************************" << endl;

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				exact(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}

		MATLAB.Variable("i", i);
		MATLAB.Variable("totalT", totalT);
		Surfactant.Variable("Surfactant");
		exact.Variable("SurfactantExact");
		MATLAB.Command("surf(X,Y,Surfactant);");
		//MATLAB.Command("subplot(1,2,2),surf(X,Y,SurfactantExact);set(gca,'fontsize',20)");
		MATLAB.Command("loss=sum(sum((Surfactant-SurfactantExact).^2.))*(Y(2)-Y(1))^2;");
		str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', L2 error  :',num2str(loss)]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
	}

}

inline void MovingInterface::SurfactantDiffusion(const int & iter)
{
	if (iter == 1)
	{
		GenerateLinearSystem1(A, 1);
		Acsr = CSR<double>(A);
	}
	if (iter == 2)
	{
		GenerateLinearSystem2(A, 1);
		Acsr = CSR<double>(A);
	}
	if (iter == 1)
	{
		OneStepSemiImplicit();
	}
	if (iter >= 2)
	{
		TwoStepSemiImplicit();
	}
}

inline void MovingInterface::OneStepSemiImplicit()
{
	//// Linear Equation
	GenerateLinearSystem1(vectorB, 1);


	CGSolver::Solver(Acsr, vectorB, tempSur);


	int index;
#pragma omp parallel for private(index)
	for (int i = Surfactant.iStartI; i <= Surfactant.iEndI; i++)
	{
		for (int j = Surfactant.jStartI; j <= Surfactant.jEndI; j++)
		{
			index = (i - (grid.iStart + 1)) + (j - (grid.jStart + 1))*(grid.iRes - 2);
			Surfactant(i, j) = tempSur(index);
		}
	}

	//// Dirichlet Boundary Condition from the Exact Solution

#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		Surfactant(i, grid.jStart) = ExactSurfactant(grid(i, grid.jStart).x, grid(i, grid.jStart).y, totalT);
		Surfactant(i, grid.jEnd) = ExactSurfactant(grid(i, grid.jEnd).x, grid(i, grid.jEnd).y, totalT);
	}
#pragma omp parallel for
	for (int j = grid.jStart; j <= grid.jEnd; j++)
	{
		Surfactant(grid.iStart, j) = ExactSurfactant(grid(grid.iStart, j).x, grid(grid.iStart, j).y, totalT);
		Surfactant(grid.iEnd, j) = ExactSurfactant(grid(grid.iEnd, j).x, grid(grid.iEnd, j).y, totalT);
	}

	AdvectionMethod2D<double>::LSPropagatingTVDRK3(levelSet, U, V, dt);
}

inline void MovingInterface::TwoStepSemiImplicit()
{
	//// Linear Equation
	GenerateLinearSystem2(vectorB, 1);
	CGSolver::Solver(Acsr, vectorB, tempSur);

	int index;
#pragma omp parallel for private(index)
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			index = (i - (grid.iStart + 1)) + (j - (grid.jStart + 1))*(grid.iRes - 2);
			Surfactant(i, j) = tempSur(index);
		}
	}

	//// Dirichlet Boundary Condition from the Exact Solution
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		Surfactant(i, grid.jStart) = ExactSurfactant(grid(i, grid.jStart).x, grid(i, grid.jStart).y, totalT);
		Surfactant(i, grid.jEnd) = ExactSurfactant(grid(i, grid.jEnd).x, grid(i, grid.jEnd).y, totalT);
	}
#pragma omp parallel for
	for (int j = grid.jStart; j <= grid.jEnd; j++)
	{
		Surfactant(grid.iStart, j) = ExactSurfactant(grid(grid.iStart, j).x, grid(grid.iStart, j).y, totalT);
		Surfactant(grid.iEnd, j) = ExactSurfactant(grid(grid.iEnd, j).x, grid(grid.iEnd, j).y, totalT);
	}

	AdvectionMethod2D<double>::LSPropagatingTVDRK3(levelSet, U, V, dt);
}

inline void MovingInterface::SurfactantNormalTerm(FD& ipField, LS& ipLevelSet, Array2D<double>& term)
{
	levelSet.ComputeMeanCurvature();
	levelSet.ComputeNormal();
	Surfactant.Gradient();

	Array2D<Vector2D<double>>& gradientF = ipField.gradient;
	VT normal;
	Array2D<double> Hessian(2, 2);

	Array2D<double>& wenoDxMinus = levelSet.phi.dfdxM;
	Array2D<double>& wenoDxPlus = levelSet.phi.dfdxP;
	Array2D<double>& wenoDyMinus = levelSet.phi.dfdyM;
	Array2D<double>& wenoDyPlus = levelSet.phi.dfdyP;
	AdvectionMethod2D<double>::WENO5thDerivation(ipField, wenoDxMinus, wenoDxPlus, wenoDyMinus, wenoDyPlus);

	U.Gradient();
	V.Gradient();
	Array2D<Vector2D<double>>& gradientU = U.gradient;
	Array2D<Vector2D<double>>& gradientV = V.gradient;

	double curvatureThreshold = 3.0;
	double curvature;
#pragma omp parallel for private(normal, Hessian, curvature)
	for (int i = term.iStart; i <= term.iEnd; i++)
	{
		for (int j = term.jStart; j <= term.jEnd; j++)
		{
			normal = ipLevelSet.normal(i, j);
			//Normal(i, j) = normal;
			Hessian = Surfactant.Hessian(i, j);
			curvature = -ipLevelSet.meanCurvature(i, j); //// LEVELSET.MEANCURVATURE has a nagative sign. so mutiple -1.
			term(i, j) = 0;
			if (abs(curvature) < curvatureThreshold)
			{
				term(i, j) += -curvature*DotProduct(normal, gradientF(i, j));
			}
			else
			{
				term(i, j) += -AdvectionMethod2D<double>::sign(curvature)*curvatureThreshold*DotProduct(normal, gradientF(i, j));
			}
			term(i, j) += -normal(0)*(Hessian(0, 0)*normal(0) + Hessian(0, 1)*normal(1));
			term(i, j) += -normal(1)*(Hessian(1, 0)*normal(0) + Hessian(1, 1)*normal(1));

			//// Upwind WENO
			term(i, j) += -(AdvectionMethod2D<double>::Plus(U(i, j))*wenoDxMinus(i, j) + AdvectionMethod2D<double>::Minus(U(i, j))*wenoDxPlus(i, j));
			term(i, j) += -(AdvectionMethod2D<double>::Plus(V(i, j))*wenoDyMinus(i, j) + AdvectionMethod2D<double>::Minus(V(i, j))*wenoDyPlus(i, j));

			////
			term(i, j) += normal(0)*(gradientU(i, j).x*normal(0) + gradientU(i, j).y*normal(1))*ipField(i, j);
			term(i, j) += normal(1)*(gradientV(i, j).x*normal(0) + gradientV(i, j).y*normal(1))*ipField(i, j);
		}
	}
	//term.Variable("addedTerm");
	//ArrayVec2DVariable("normal", Normal.dataArray);
	//ArrayVec2DVariable("gradient", gradient.dataArray);
}

inline double MovingInterface::ExactSurfactant(const double & x, const double & y, const double & time)
{
	if (ExamNum == 1 || ExamNum == 3)
	{
		return exp(-time / (x*x + y*y + DBL_EPSILON))*(y / sqrt(x*x + y*y + DBL_EPSILON)) + 2;
	}
	if (ExamNum >= 2)
	{
		return exp(-time / ((x - totalT)*(x - totalT) + y*y + DBL_EPSILON))*(y / sqrt((x - totalT)*(x - totalT) + y*y + DBL_EPSILON)) + 2;
	}
	return 0.0;
}

inline void MovingInterface::GenerateLinearSystem1(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int innerIStart = grid.iStart + 1;
	int innerIEnd = grid.iEnd - 1;
	int innerJStart = grid.jStart + 1;
	int innerJEnd = grid.jEnd - 1;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

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
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
			}
			else if (j > innerJStart && j < innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;

				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;

				}
				else
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
			}
			else if (j == innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
			}
		}
	}
}

inline void MovingInterface::GenerateLinearSystem1(VectorND<double>& vectorB, const double & scaling)
{
	int innerIStart = grid.iStart + 1;
	int innerIEnd = grid.iEnd - 1;
	int innerJStart = grid.jStart + 1;
	int innerJEnd = grid.jEnd - 1;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

	SurfactantNormalTerm(Surfactant, levelSet, term);
	int index;
#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;
			vectorB(index) = Surfactant(i, j) + dt*term(i, j);

			if (i == innerIStart)
			{
				vectorB(index) += dt*grid.oneOverdx2*ExactSurfactant(grid(i - 1, j).x, grid(i - 1, j).y, totalT);
			}
			if (i == innerIEnd)
			{
				vectorB(index) += dt*grid.oneOverdx2*ExactSurfactant(grid(i + 1, j).x, grid(i + 1, j).y, totalT);
			}
			if (j == innerJStart)
			{
				vectorB(index) += dt*grid.oneOverdy2*ExactSurfactant(grid(i, j - 1).x, grid(i, j - 1).y, totalT);

			}
			if (j == innerJEnd)
			{
				vectorB(index) += dt*grid.oneOverdy2*ExactSurfactant(grid(i, j + 1).x, grid(i, j + 1).y, totalT);
			}
			vectorB(index) *= scaling;
		}
	}
}

inline void MovingInterface::GenerateLinearSystem2(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int innerIStart = grid.iStart + 1;
	int innerIEnd = grid.iEnd - 1;
	int innerJStart = grid.jStart + 1;
	int innerJEnd = grid.jEnd - 1;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

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
			//// Boundary condition.
			if (j == innerJStart)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1 / 2.0 * dt*grid.oneOverdy2;
				}
			}
			else if (j > innerJStart && j < innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;

				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;

				}
				else
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
			}
			else if (j == innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
			}
		}
	}
}

inline void MovingInterface::GenerateLinearSystem2(VectorND<double>& vectorB, const double & scaling)
{
	int innerIStart = grid.iStart + 1;
	int innerIEnd = grid.iEnd - 1;
	int innerJStart = grid.jStart + 1;
	int innerJEnd = grid.jEnd - 1;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

	termOld = term;
	SurfactantNormalTerm(Surfactant, levelSet, term);
	//SurfactantNormalTerm(SurfactantOld, levelSetOld, termOld);

	int index;
#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;
			vectorB(index) = Surfactant(i, j) + dt / 2.0*(Surfactant.dxxPhi(i, j) + Surfactant.dyyPhi(i, j))
				+ 3.0 / 2.0*dt*term(i, j) - 1.0 / 2.0*dt*termOld(i, j);

			if (i == innerIStart)
			{
				vectorB(index) += dt / 2.0*grid.oneOverdx2*ExactSurfactant(grid(i - 1, j).x, grid(i - 1, j).y, totalT);
			}
			if (i == innerIEnd)
			{
				vectorB(index) += dt / 2.0*grid.oneOverdx2*ExactSurfactant(grid(i + 1, j).x, grid(i + 1, j).y, totalT);
			}
			if (j == innerJStart)
			{
				vectorB(index) += dt / 2.0*grid.oneOverdy2*ExactSurfactant(grid(i, j - 1).x, grid(i, j - 1).y, totalT);

			}
			if (j == innerJEnd)
			{
				vectorB(index) += dt / 2.0*grid.oneOverdy2*ExactSurfactant(grid(i, j + 1).x, grid(i, j + 1).y, totalT);
			}
			vectorB(index) *= scaling;
		}
	}
}


/////////////////////////////////////////////////////////////////////////////////////////
////
////       "Local" Surfactant Diffusion Solver
////
/////////////////////////////////////////////////////////////////////////////////////////
inline void MovingInterface::LSurfactantDiffusionSolver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	const char* cmd;


	InitialCondition(example);

	grid.Variable();
	FD exact(grid);

	levelSet.phi.Variable("phi0");
	levelSet.tube.Variable("Tube");
	Surfactant.Variable("Surfactant0");

	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	MATLAB.Command("SurTube1 = Surfactant0.*(Tube<=1);");
	MATLAB.Command("surf(X,Y,SurTube1),set(gca,'fontsize',20)");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);
	int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.2*min(grid.dx, grid.dy)));;
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << "*******************************************************************" << endl;
		cout << "       Iteration " << to_string(i) << " : Start" << endl;
		totalT += dt;
		cout << "Diffusion Start" << endl;
		LSurfactantDiffusion(i);
		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);

		cout << "Diffusion End" << endl;


#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				exact(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}

		Surfactant.Variable("Surfactant");
		exact.Variable("exact");
		MATLAB.Command("SurTube1 = Surfactant.*(Tube<=1);");
		MATLAB.Command("surf(X,Y,SurTube1), hold on, contour(X,Y,Tube),hold off,set(gca,'fontsize',20);");
		MATLAB.Command("loss = sum(sum(SurTube1-exact.*(Tube==1)).^2)*(Y(2)-Y(1))*(Y(2)-Y(1));");
		MATLAB.Variable("i", i);
		MATLAB.Variable("totalT", totalT);
		str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', L2 error  :',num2str(loss)]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
	}
}


/////////////////////////////////////////////////////////////////////////////////////////
////
////       Interface Moving with Local Surfactant Diffusion
////
/////////////////////////////////////////////////////////////////////////////////////////
inline void MovingInterface::EulerianMovingInterfaceSolver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	const char* cmd;


	InitialCondition(example);

	grid.Variable();
	levelSet.phi.Variable("phi0");
	levelSet.tube.Variable("Tube");
	Surfactant.Variable("Surfactant");

	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	//MATLAB.Command("subplot(2,1,1),surf(X,Y,Surfactant),subplot(2,1,2),contour(X,Y,phi0,[0 0],'b'),grid on,axis equal");
	MATLAB.Command("SurTube1 = Surfactant.*(Tube<=1);");
	MATLAB.Command("surf(X,Y,SurTube1),grid on,axis equal,set(gca,'fontsize',20)");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);

	FD exact(grid);
	int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy)) * 2;
	int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.2*min(grid.dx, grid.dy)));
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << "*******************************************************************" << endl;
		cout << "       Iteration " << to_string(i) << " : Start" << endl;
		totalT += dt;
		cout << "Diffusion Start" << endl;
		LSurfactantDiffusion(i);
		cout << "Diffusion End" << endl;

		Surfactant.Variable("Surfactant");
		levelSet.tube.Variable("Tube");
		MATLAB.Command("SurTube1 = Surfactant.*(Tube<=1);");
		MATLAB.Command("surf(X,Y,SurTube1),axis equal,axis([X(1) X(end) Y(1) Y(end)]), hold on, contour(X,Y,Tube,'r'),hold off;set(gca,'fontsize',20)");

		if (ExamNum == 4)
		{
#pragma omp parallel for
			for (int i = grid.iStart; i <= grid.iEnd; i++)
			{
				for (int j = grid.jStart; j <= grid.jEnd; j++)
				{
					exact(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
				}
			}

			exact.Variable("exact");
			MATLAB.Command("loss = sum(sum(SurTube1-exact.*(Tube==1)).^2)*(Y(2)-Y(1))*(Y(2)-Y(1));");
			MATLAB.Variable("i", i);
			MATLAB.Variable("totalT", totalT);
			str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', L2 error  :',num2str(loss)]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
		}
		else
		{
			if (i == 1)
			{
#pragma omp parallel for
				for (int i = grid.iStart; i <= grid.iEnd; i++)
				{
					for (int j = grid.jStart; j <= grid.jEnd; j++)
					{
						exact(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, 0);
					}
				}
				exact.Variable("initial");
				MATLAB.Command("initialF = sum(sum(initial.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");
			}
			MATLAB.Command("SurIntegral = sum(sum(SurTube1))*(Y(2)-Y(1))*(Y(2)-Y(1));");
			MATLAB.Variable("i", i);
			MATLAB.Variable("totalT", totalT);
			MATLAB.Command("loss=(SurIntegral-initialF)/initialF*100;");
			str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', loss(%)  :',num2str(loss)]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
		}


		AdvectionMethod2D<double>::LLSPropagatingTVDRK3(levelSet, U, V, dt);
		AdvectionMethod2D<double>::LLSReinitializationTVDRK3(levelSet, dt, reinitialIter);
		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);
		levelSet.UpdateInterface();
		levelSet.UpdateLLS();

		cout << "       Iteration " << to_string(i) << " : End" << endl;

		//levelSet.phi.Variable("phi");
		//MATLAB.Command("subplot(2,1,2),contour(X,Y,phi0,[0 0],'b'),grid on,axis equal,hold on,contour(X,Y,Tube),hold off;");
		//str = string("title(['iteration : ', num2str(") + to_string(i) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
		//cmd = str.c_str();
		//MATLAB.Command(cmd);
	}

}

inline void MovingInterface::SurfactantTube2Extrapolation()
{
	int i, j;
	double temp;
#pragma omp parallel for private(i, j, temp)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) == 2)
		{
			temp = Surfactant(i, j);
			Surfactant(i, j) = 2 * Surfactant(i, j) - Surfactant.dataArrayOld(i, j);
			Surfactant.dataArrayOld(i, j) = temp;
		}
	}
}

inline void MovingInterface::LSurfactantDiffusion(const int & iter)
{
	A = Array2D<double>(1, levelSet.numTube1, 1, levelSet.numTube1);

	if (iter == 1)
	{
		LGenerateLinearSystem1(A, 1);

		Acsr = CSR<double>(A);

		A.Delete();

		LOneStepSemiImplicit();
	}

	if (iter >= 2)
	{
		SurfactantTube2Extrapolation();

		LGenerateLinearSystem2(A, 1);

		Acsr = CSR<double>(A);

		A.Delete();

		LTwoStepSemiImplicit();
	}

}

inline void MovingInterface::LOneStepSemiImplicit()
{
	//// Linear Equation
	vectorB = VectorND<double>(levelSet.numTube1);
	LGenerateLinearSystem1(vectorB, 1);

	tempSur = VectorND<double>(levelSet.numTube1);

	CGSolver::Solver(Acsr, vectorB, tempSur);


	int i, j;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;
		if (levelSet.tube(i, j) == 3)
		{
			Surfactant(i, j) = 0;
		}
		else
		{
			Surfactant(i, j) = tempSur(k - 1);
		}
	}
}

inline void MovingInterface::LTwoStepSemiImplicit()
{
	//// Linear Equation
	termOld = term;
	vectorB = VectorND<double>(levelSet.numTube1);
	LGenerateLinearSystem2(vectorB, 1);

	tempSur = VectorND<double>(levelSet.numTube1);

	CGSolver::Solver(Acsr, vectorB, tempSur);

	int i, j;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;
		Surfactant(i, j) = tempSur(k - 1);
	}

#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		i = levelSet.tubeIndex(k).i;
		j = levelSet.tubeIndex(k).j;
		if (levelSet.tube(i, j) == 3)
		{
			Surfactant(i, j) = 0;
		}
	}
}

inline void MovingInterface::LSurfactantNormalTerm(FD & ipField, LS & ipLevelSet, Array2D<double>& term)
{
	levelSet.LComputeMeanCurvature(1);
	levelSet.LComputeUnitNormal(1);
	Surfactant.Gradient();

	Array2D<VT>& gradientF = ipField.gradient;
	VT normal;
	Array2D<double> Hessian(2, 2);

	Array2D<double>& wenoDxMinus = levelSet.phi.dfdxM;
	Array2D<double>& wenoDxPlus = levelSet.phi.dfdxP;
	Array2D<double>& wenoDyMinus = levelSet.phi.dfdyM;
	Array2D<double>& wenoDyPlus = levelSet.phi.dfdyP;
	AdvectionMethod2D<double>::LLSWENO5thDerivation(ipLevelSet, ipField, wenoDxMinus, wenoDxPlus, wenoDyMinus, wenoDyPlus);

	U.Gradient();
	V.Gradient();
	Array2D<VT>& gradientU = U.gradient;
	Array2D<VT>& gradientV = V.gradient;

	double curvatureThreshold = 3.0;
	double curvature;
	int i, j;
#pragma omp parallel for private(i,j, normal, Hessian, curvature)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) == 1)
		{
			normal = ipLevelSet.unitNormal(i, j);
			//Normal(i, j) = normal;
			Hessian = Surfactant.Hessian(i, j);
			curvature = -ipLevelSet.meanCurvature(i, j); //// LEVELSET.MEANCURVATURE has a nagative sign. so mutiple -1.
			term(i, j) = 0;
			if (abs(curvature) < curvatureThreshold)
			{
				term(i, j) += -curvature*DotProduct(normal, gradientF(i, j));
			}
			else
			{
				term(i, j) += -AdvectionMethod2D<double>::sign(curvature)*curvatureThreshold*DotProduct(normal, gradientF(i, j));
			}
			term(i, j) += -normal(0)*(Hessian(0, 0)*normal(0) + Hessian(0, 1)*normal(1));
			term(i, j) += -normal(1)*(Hessian(1, 0)*normal(0) + Hessian(1, 1)*normal(1));

			//// Upwind WENO
			term(i, j) += -(AdvectionMethod2D<double>::Plus(U(i, j))*wenoDxMinus(i, j) + AdvectionMethod2D<double>::Minus(U(i, j))*wenoDxPlus(i, j));
			term(i, j) += -(AdvectionMethod2D<double>::Plus(V(i, j))*wenoDyMinus(i, j) + AdvectionMethod2D<double>::Minus(V(i, j))*wenoDyPlus(i, j));

			////
			term(i, j) += normal(0)*(gradientU(i, j).x*normal(0) + gradientU(i, j).y*normal(1))*ipField(i, j);
			term(i, j) += normal(1)*(gradientV(i, j).x*normal(0) + gradientV(i, j).y*normal(1))*ipField(i, j);
		}
	}
}

inline void MovingInterface::LGenerateLinearSystem1(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	VI leftIndex, rightIndex, bottomIndex, topIndex;
	int i, j, l, m, n;
#pragma omp parallel for private(i, j, l, m, n, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;

		matrixA(k, k) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);

		leftIndex = VI(i - 1, j);
		l = levelSet.tubeIJ2K(leftIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 1)
		{
			l = levelSet.tube1(m, n);
			matrixA(k, l) = scaling * -1 * dt*grid.oneOverdx2;
		}

		rightIndex = VI(i + 1, j);
		l = levelSet.tubeIJ2K(rightIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 1)
		{
			l = levelSet.tube1(m, n);
			matrixA(k, l) = scaling * -1 * dt*grid.oneOverdx2;
		}

		bottomIndex = VI(i, j - 1);
		l = levelSet.tubeIJ2K(bottomIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 1)
		{
			l = levelSet.tube1(m, n);
			matrixA(l, k) = scaling * -1 * dt*grid.oneOverdy2;
		}

		topIndex = VI(i, j + 1);
		l = levelSet.tubeIJ2K(topIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 1)
		{
			l = levelSet.tube1(m, n);
			matrixA(l, k) = scaling * -1 * dt*grid.oneOverdy2;
		}
	}
}

inline void MovingInterface::LGenerateLinearSystem1(VectorND<double>& vectorB, const double & scaling)
{
	LSurfactantNormalTerm(Surfactant, levelSet, term);
	VI leftIndex, rightIndex, bottomIndex, topIndex;
	int i, j, l, m, n;
#pragma omp parallel for private(i, j, l, m, n, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;
		vectorB(k - 1) = Surfactant(i, j) + dt*term(i, j);

		leftIndex = VI(i - 1, j);
		l = levelSet.tubeIJ2K(leftIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 2)
		{
			vectorB(k - 1) += dt*grid.oneOverdx2*Surfactant(m, n);
		}

		rightIndex = VI(i + 1, j);
		l = levelSet.tubeIJ2K(rightIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 2)
		{
			vectorB(k - 1) += dt*grid.oneOverdx2*Surfactant(m, n);
		}

		bottomIndex = VI(i, j - 1);
		l = levelSet.tubeIJ2K(bottomIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 2)
		{
			vectorB(k - 1) += dt*grid.oneOverdy2*Surfactant(m, n);
		}

		topIndex = VI(i, j + 1);
		l = levelSet.tubeIJ2K(topIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 2)
		{
			vectorB(k - 1) += dt*grid.oneOverdy2*Surfactant(m, n);
		}
		vectorB(k - 1) *= scaling;
	}
}

inline void MovingInterface::LGenerateLinearSystem2(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	VI leftIndex, rightIndex, bottomIndex, topIndex;
	int i, j, l, m, n;
#pragma omp parallel for private(i, j, l, m, n, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;

		matrixA(k, k) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);

		leftIndex = VI(i - 1, j);
		l = levelSet.tubeIJ2K(leftIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 1)
		{
			l = levelSet.tube1(m, n);
			matrixA(k, l) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
		}

		rightIndex = VI(i + 1, j);
		l = levelSet.tubeIJ2K(rightIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 1)
		{
			l = levelSet.tube1(m, n);
			matrixA(k, l) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
		}

		bottomIndex = VI(i, j - 1);
		l = levelSet.tubeIJ2K(bottomIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 1)
		{
			l = levelSet.tube1(m, n);
			matrixA(l, k) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
		}

		topIndex = VI(i, j + 1);
		l = levelSet.tubeIJ2K(topIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 1)
		{
			l = levelSet.tube1(m, n);
			matrixA(l, k) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
		}
	}
}

inline void MovingInterface::LGenerateLinearSystem2(VectorND<double>& vectorB, const double & scaling)
{
	LSurfactantNormalTerm(Surfactant, levelSet, term);
	VI leftIndex, rightIndex, bottomIndex, topIndex;
	int i, j, l, m, n;
#pragma omp parallel for private(i, j, l, m, n, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;
		vectorB(k - 1) = Surfactant(i, j) + dt / 2.0*(Surfactant.dxxPhi(i, j) + Surfactant.dyyPhi(i, j))
			+ 3.0 / 2.0*dt*term(i, j) - 1.0 / 2.0*dt*termOld(i, j);

		leftIndex = VI(i - 1, j);
		l = levelSet.tubeIJ2K(leftIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 2)
		{
			vectorB(k - 1) += dt / 2.0*grid.oneOverdx2*Surfactant(m, n);
		}

		rightIndex = VI(i + 1, j);
		l = levelSet.tubeIJ2K(rightIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 2)
		{
			vectorB(k - 1) += dt / 2.0*grid.oneOverdx2*Surfactant(m, n);
		}

		bottomIndex = VI(i, j - 1);
		l = levelSet.tubeIJ2K(bottomIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 2)
		{
			vectorB(k - 1) += dt / 2.0*grid.oneOverdy2*Surfactant(m, n);
		}

		topIndex = VI(i, j + 1);
		l = levelSet.tubeIJ2K(topIndex);
		levelSet.TubeIndex(l, m, n);
		if (levelSet.tube(m, n) == 2)
		{
			vectorB(k - 1) += dt / 2.0*grid.oneOverdy2*Surfactant(m, n);
		}
		vectorB(k - 1) *= scaling;
	}
}

inline void MovingInterface::DimlessNonlinearLangmu1rEOS()
{
	int i, j;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) <= 1)
		{
			SurfaceTension(i, j) = 1 + Fluid.El*log(1 - Fluid.Xi* Surfactant(i, j));
		}
		else
		{
			SurfaceTension(i, j) = 1;
		}

	}
}

inline void MovingInterface::DimlessNonlinearLangmu1rEOS(const int & tubeRange)
{
	int i, j;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) <= tubeRange)
		{
			SurfaceTension(i, j) = 1 + Fluid.El*log(1 - Fluid.Xi* Surfactant(i, j));
		}
		else
		{
			SurfaceTension(i, j) = 1;
		}

	}
}

inline void MovingInterface::DimlessLinearLangmu1rEOS()
{
	int i, j;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) <= 1)
		{
			SurfaceTension(i, j) = 1 - Fluid.El* Fluid.Xi* Surfactant(i, j);
		}
		else
		{
			SurfaceTension(i, j) = 1;
		}

	}
}

inline void MovingInterface::DimlessLinearLangmu1rEOS(const int & tubeRange)
{
	int i, j;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) <= tubeRange)
		{
			SurfaceTension(i, j) = 1 - Fluid.El* Fluid.Xi* Surfactant(i, j);
		}
		else
		{
			SurfaceTension(i, j) = 1;
		}

	}
}

inline void MovingInterface::LinearLangmu1rEOS()
{
	int i, j;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) <= 1)
		{
			SurfaceTension(i, j) = gamma0 - IdealGasConstant*AbsTemperature*Surfactant(i, j);
		}
		else
		{
			SurfaceTension(i, j) = gamma0;
		}

	}
}

inline double MovingInterface::IntegralSurfactant()
{
	double  tempIntegral = 0;
	int i, j;
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) <= 1)
		{
			tempIntegral += Surfactant(i, j)*grid.dx*grid.dy;
		}
	}
	return tempIntegral;
}

inline void MovingInterface::ConserveSurfactantFactorBeta()
{
	double currentSurfactant = IntegralSurfactant();
	int i, j;
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) <= 1)
		{
			Surfactant(i, j) = Surfactant(i, j)*initialSurfactant / currentSurfactant;
		}
	}
}

