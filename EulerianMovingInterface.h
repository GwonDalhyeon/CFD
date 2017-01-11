#pragma once

#include "AdvectionMethod2D.h"
//#include "EulerianFluidSolver.h"
#include "FluidSolver2D.h"


class MovingInterface
{
public:
	MovingInterface(FluidSolver2D & ipFluid);
	~MovingInterface();

	int ExamNum;

	FluidSolver2D& Fluid;

	Grid2D& grid = Fluid.grid;

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
	// CG solver 1
	CSR<double> Acsr;


	double& reynoldNum = Fluid.reynoldNum;
	double& cflCondition = Fluid.cflCondition;
	double& Ca = Fluid.Ca;
	double& Xi = Fluid.Xi;
	double& El = Fluid.El;
	double& Pe = Fluid.Pe;
	double& dt = Fluid.dt;
	double totalT;

	int& maxIteration = Fluid.maxIteration;
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
	inline void GenerateLinearSystem1(CSR<double>& ipCSR);
	inline void GenerateLinearSystem1(VectorND<double>& vectorB);
	inline void GenerateLinearSystem2(CSR<double>& ipCSR);
	inline void GenerateLinearSystem2(VectorND<double>& vectorB);

	inline void PlotSurfactant();

	//////////////////////////////////////////////
	////   Local Surfactant Diffusion Solver  ////
	//////////////////////////////////////////////
	inline void LSurfactantDiffusionSolver(const int& example);

	inline void PlotLocalSurfactant();

	///////////////////////////////////////////
	//// Eulerian Moving Interface Solver  ////
	///////////////////////////////////////////
	inline void EulerianMovingInterfaceSolver(const int& example);

	inline void SurfactantTube2Extrapolation();

	inline void LSurfactantDiffusion(const int& iter);
	inline void LCountNonZero();
	inline void LOneStepSemiImplicit();
	inline void LTwoStepSemiImplicit();
	inline void LSurfactantNormalTerm(FD& ipField, LS& ipLevelSet, Array2D<double>& term);
	inline void LGenerateLinearSystem1(CSR<double>& ipCSR);
	inline void LGenerateLinearSystem1(VectorND<double>& vectorB);
	inline void LGenerateLinearSystem2(CSR<double>& ipCSR);
	inline void LGenerateLinearSystem2(VectorND<double>& vectorB);

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

inline MovingInterface::MovingInterface(FluidSolver2D & ipFluid)
	:Fluid(ipFluid)
{
}

MovingInterface::~MovingInterface()
{
}

inline void MovingInterface::InitialCondition(const int & example)
{
	AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

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
		int tempBC = 0;
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			for (int i = grid.iStart; i <= grid.iEnd; i++)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, 0);
				if (i == grid.iStart || i == grid.iEnd || j == grid.jStart || j == grid.jEnd)
				{
					Surfactant.BC(i, j) = BC_DIR;
					continue;
				}
				Surfactant.BC(i, j) = tempBC++;
			}
		}

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);

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
		int tempBC = 0;
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			for (int i = grid.iStart; i <= grid.iEnd; i++)
			{
				x0 = grid(i, j).x;
				y0 = grid(i, j).y;
				Surfactant(i, j) = ExactSurfactant(x0, y0, 0);
				if (i == grid.iStart || i == grid.iEnd || j == grid.jStart || j == grid.jEnd)
				{
					Surfactant.BC(i, j) = BC_DIR;
					continue;
				}
				Surfactant.BC(i, j) = tempBC++;
			}
		}
		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);

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
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius - 100*DBL_EPSILON;
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
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
				//Surfactant(i, j) = sqrt((grid(i, j).x)*(grid(i, j).x) + (grid(i, j).y)*(grid(i, j).y)); ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}

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
		initialSurfactant = IntegralSurfactant();

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);
		cflCondition = 1.0 / 4.0;
		dt = AdvectionMethod2D<double>::AdaptiveTimeStep(U, cflCondition);
		maxIteration = 800;
		totalT = 0;

		Pe = 10;
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
			if (levelSet.tube(i, j) <= 2)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}

		initialSurfactant = IntegralSurfactant();

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
			if (levelSet.tube(i, j) <= 2)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}

		initialSurfactant = IntegralSurfactant();

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);
		cflCondition = 1.0 / 4.0;
		dt = cflCondition*min(grid.dx, grid.dy);
		double finalT = 2;
		maxIteration = ceil(finalT / dt);
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
		double x, y;
#pragma omp parallel for private(x, y)
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				x = grid(i, j).x;
				y = grid(i, j).y;
				levelSet(i, j) = sqrt(x*x + y*y) - radius;
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
		double El = Fluid.El;
		double Xi = Fluid.Xi;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (levelSet.tube(i, j) <= 2)
				{
					SurfaceTension(i, j) = 1 + El*log(1 - Xi* Surfactant(i, j));
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
		//cflCondition = 1.0 / 4.0;
		//dt = cflCondition*min(grid.dx, grid.dy);
		//totalT = 0;
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

		term = Array2D<double>(grid);
		termOld = Array2D<double>(grid);
		totalT = 0;
	}

	AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

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


	InitialCondition(example);

	grid.Variable();

	//// Initial Surfactant
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	PlotSurfactant();

	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(totalT) + string(")]), ;");
	MATLAB.Command(str.c_str());

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

		PlotSurfactant();

		MATLAB.Variable("i", i);
		MATLAB.Variable("totalT", totalT);
		exact.Variable("SurfactantExact");
		//MATLAB.Command("subplot(1,2,2),surf(X,Y,SurfactantExact);set(gca,'fontsize',20)");
		MATLAB.Command("loss=sum(sum((Surfactant-SurfactantExact).^2.))*(Y(2)-Y(1))^2;");
		str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', L2 error  :',num2str(loss)]);");
		MATLAB.Command(str.c_str());
	}

}

inline void MovingInterface::SurfactantDiffusion(const int & iter)
{
	if (iter == 1)
	{
		GenerateLinearSystem1(Acsr);
	}
	if (iter == 2)
	{
		GenerateLinearSystem2(Acsr);
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
	GenerateLinearSystem1(vectorB);


	CGSolver::Solver(Acsr, vectorB, tempSur);

	Fluid.VectorToGrid(tempSur, Surfactant);

	AdvectionMethod2D<double>::LSPropagatingTVDRK3(levelSet, U, V, dt);
}

inline void MovingInterface::TwoStepSemiImplicit()
{
	//// Linear Equation
	GenerateLinearSystem2(vectorB);

	CGSolver::Solver(Acsr, vectorB, tempSur);

	Fluid.VectorToGrid(tempSur, Surfactant);

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
	double curvature, vel;
	VT velG;
#pragma omp parallel for private(normal, Hessian, curvature, vel, velG)
	for (int i = term.iStart; i <= term.iEnd; i++)
	{
		for (int j = term.jStart; j <= term.jEnd; j++)
		{
			if (Surfactant.BC(i,j) < 0)
			{
				continue;
			}
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
			vel = 0.5 * (U(i, j) + U(i + 1, j));
			term(i, j) += -(AdvectionMethod2D<double>::Plus(vel)*wenoDxMinus(i, j) + AdvectionMethod2D<double>::Minus(vel)*wenoDxPlus(i, j));
			vel = 0.5 * (V(i, j) + V(i, j + 1));
			term(i, j) += -(AdvectionMethod2D<double>::Plus(vel)*wenoDyMinus(i, j) + AdvectionMethod2D<double>::Minus(vel)*wenoDyPlus(i, j));

			////
			velG = 0.5 * (gradientU(i, j) + gradientU(i + 1, j));
			term(i, j) += normal(0)*(velG.x*normal(0) + velG.y*normal(1))*ipField(i, j);
			velG = 0.5 * (gradientV(i, j) + gradientV(i, j + 1));
			term(i, j) += normal(1)*(velG.x*normal(0) + velG.y*normal(1))*ipField(i, j);
		}
	}
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

inline void MovingInterface::GenerateLinearSystem1(CSR<double>& ipCSR)
{
	Surfactant.CountNonZero();
	vectorB = VectorND<double>(Surfactant.num_all_full_cells);
	tempSur = VectorND<double>(Surfactant.num_all_full_cells);
	Acsr = CSR<double>(Surfactant.num_all_full_cells, Surfactant.nnz);


	int iStart = Surfactant.iStart, iEnd = Surfactant.iEnd, jStart = Surfactant.jStart, jEnd = Surfactant.jEnd;

	double dx = Surfactant.dx, dy = Surfactant.dy;
	double dx2 = dx*dx, dy2 = dy*dy, dxdy = dx*dy;
	double oneOverdx = 1 / dx, oneOverdx2 = 1 / dx2;
	double oneOverdy = 1 / dy, oneOverdy2 = 1 / dy2;

	Array2D<int>& BC = Surfactant.BC;

	double coefIJ;

	for (int j = jStart; j <= jEnd; j++)
	{
		for (int i = iStart; i <= iEnd; i++)
		{
			coefIJ = 0;

			if (BC(i, j) < 0)
			{
				continue;
			}

			//// If neighbor is full cell
			if (i>iStart)
			{
				if (BC(i - 1, j) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i - 1, j), -oneOverdx2);
				}
			}
			if (i<iEnd)
			{
				if (BC(i + 1, j) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i + 1, j), -oneOverdx2);
				}
			}
			if (j>jStart)
			{
				if (BC(i, j - 1) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i, j - 1), -oneOverdy2);
				}
			}
			if (j<jEnd)
			{
				if (BC(i, j + 1) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i, j + 1), -oneOverdy2);
				}
			}

			//// Dirichlet Boundary Condition
			if (i>iStart)
			{
				if (BC(i - 1, j) == BC_DIR)	coefIJ++;
			}
			if (i<iEnd)
			{
				if (BC(i + 1, j) == BC_DIR)	coefIJ++;
			}
			if (j>jStart)
			{
				if (BC(i, j - 1) == BC_DIR)	coefIJ++;
			}
			if (j<jEnd)
			{
				if (BC(i, j + 1) == BC_DIR)	coefIJ++;
			}

			if ((i == Surfactant.iStartI) && (j == Surfactant.jStartI))
			{
				//if (BC(i - 1, j) == BC_NEUM)
				//{
				//	coefIJ ++;
				//}
				//if (BC(i + 1, j) == BC_NEUM)
				//{
				//	coefIJ ++;
				//}
				//if (BC(i, j - 1) == BC_NEUM)
				//{
				//	coefIJ ++;
				//}
				//if (BC(i, j + 1) == BC_NEUM)
				//{
				//	coefIJ ++;
				//}
			}
			else
			{
				if (i>iStart)
				{
					if (BC(i - 1, j) == BC_NEUM) coefIJ += 0;
				}
				if (i<iEnd)
				{
					if (BC(i + 1, j) == BC_NEUM) coefIJ += 0;
				}
				if (j>jStart)
				{
					if (BC(i, j - 1) == BC_NEUM) coefIJ += 0;
				}
				if (j<jEnd)
				{
					if (BC(i, j + 1) == BC_NEUM) coefIJ += 0;
				}
			}

			if (coefIJ == 0)
			{
				coefIJ = 1;
			}

			ipCSR.AssignValue(BC(i, j), BC(i, j), oneOverdx2*coefIJ + 1./dt);

		}
	}

}

inline void MovingInterface::GenerateLinearSystem1(VectorND<double>& vectorB)
{
	Array2D<int>& BC = Surfactant.BC;
	SurfactantNormalTerm(Surfactant, levelSet, term);
	double oneOverdt = 1 / dt;
#pragma omp parallel for
	for (int j = Surfactant.jStart; j <= Surfactant.jEnd; j++)
	{
		for (int i = Surfactant.iStart; i <= Surfactant.iEnd; i++)
		{
			if (BC(i, j) < 0)
			{
				continue;
			}
			vectorB(BC(i, j)) = Surfactant(i, j)*oneOverdt + term(i, j);

			if (i > Surfactant.iStart)
			{
				if (BC(i - 1, j) == BC_DIR)
				{
					Surfactant(i - 1, j) = ExactSurfactant(grid(i - 1, j).x, grid(i - 1, j).y, totalT);
					vectorB(BC(i, j)) += Surfactant(i - 1, j)*Surfactant.oneOverdx2;
				}
			}
			if (i < Surfactant.iEnd)
			{
				if (BC(i + 1, j) == BC_DIR)
				{
					Surfactant(i + 1, j) = ExactSurfactant(grid(i + 1, j).x, grid(i + 1, j).y, totalT);
					vectorB(BC(i, j)) += Surfactant(i + 1, j)*Surfactant.oneOverdx2;
				}
			}
			if (j > Surfactant.jStart)
			{
				if (BC(i, j - 1) == BC_DIR)
				{
					Surfactant(i, j - 1)= ExactSurfactant(grid(i, j - 1).x, grid(i, j - 1).y, totalT);
					vectorB(BC(i, j)) += Surfactant(i, j - 1)*Surfactant.oneOverdy2;
				}
			}
			if (j < Surfactant.jEnd)
			{
				if (BC(i, j + 1) == BC_DIR)
				{
					Surfactant(i, j + 1)= ExactSurfactant(grid(i, j + 1).x, grid(i, j + 1).y, totalT);
					vectorB(BC(i, j)) += Surfactant(i, j + 1)*Surfactant.oneOverdy2;
				}
			}
		}
	}
	
	Surfactant(Surfactant.iStart, Surfactant.jStart) = ExactSurfactant(grid(grid.iStart, grid.jStart).x, grid(grid.iStart, grid.jStart).y, totalT);
	Surfactant(Surfactant.iStart, Surfactant.jEnd) = ExactSurfactant(grid(grid.iStart, grid.jEnd).x, grid(grid.iStart, grid.jEnd).y, totalT);
	Surfactant(Surfactant.iEnd, Surfactant.jStart) = ExactSurfactant(grid(grid.iEnd, grid.jStart).x, grid(grid.iEnd, grid.jStart).y, totalT);
	Surfactant(Surfactant.iEnd, Surfactant.jEnd) = ExactSurfactant(grid(grid.iEnd, grid.jEnd).x, grid(grid.iEnd, grid.jEnd).y, totalT);
}

inline void MovingInterface::GenerateLinearSystem2(CSR<double>& ipCSR)
{
	int iStart = Surfactant.iStart, iEnd = Surfactant.iEnd, jStart = Surfactant.jStart, jEnd = Surfactant.jEnd;

	double dx = Surfactant.dx, dy = Surfactant.dy;
	double dx2 = dx*dx, dy2 = dy*dy, dxdy = dx*dy;
	double oneOverdx = 1 / dx, oneOverdx2 = 1 / dx2;
	double oneOverdy = 1 / dy, oneOverdy2 = 1 / dy2;
	double oneOverdt = 1 / dt;
	Array2D<int>& BC = Surfactant.BC;

	double coefIJ;

	for (int j = jStart; j <= jEnd; j++)
	{
		for (int i = iStart; i <= iEnd; i++)
		{
			coefIJ = 0;

			if (BC(i, j) < 0)
			{
				continue;
			}

			//// If neighbor is full cell
			if (i>iStart)
			{
				if (BC(i - 1, j) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i - 1, j), -oneOverdx2 / 2);
				}
			}
			if (i<iEnd)
			{
				if (BC(i + 1, j) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i + 1, j), -oneOverdx2 / 2);
				}
			}
			if (j>jStart)
			{
				if (BC(i, j - 1) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i, j - 1), -oneOverdy2 / 2);
				}
			}
			if (j<jEnd)
			{
				if (BC(i, j + 1) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i, j + 1), -oneOverdy2 / 2);
				}
			}

			//// Dirichlet Boundary Condition
			if (i>iStart)
			{
				if (BC(i - 1, j) == BC_DIR)	coefIJ++;
			}
			if (i<iEnd)
			{
				if (BC(i + 1, j) == BC_DIR)	coefIJ++;
			}
			if (j>jStart)
			{
				if (BC(i, j - 1) == BC_DIR)	coefIJ++;
			}
			if (j<jEnd)
			{
				if (BC(i, j + 1) == BC_DIR)	coefIJ++;
			}

			if ((i == Surfactant.iStartI) && (j == Surfactant.jStartI))
			{
				//if (BC(i - 1, j) == BC_NEUM)
				//{
				//	coefIJ ++;
				//}
				//if (BC(i + 1, j) == BC_NEUM)
				//{
				//	coefIJ ++;
				//}
				//if (BC(i, j - 1) == BC_NEUM)
				//{
				//	coefIJ ++;
				//}
				//if (BC(i, j + 1) == BC_NEUM)
				//{
				//	coefIJ ++;
				//}
			}
			else
			{
				if (i>iStart)
				{
					if (BC(i - 1, j) == BC_NEUM) coefIJ += 0;
				}
				if (i<iEnd)
				{
					if (BC(i + 1, j) == BC_NEUM) coefIJ += 0;
				}
				if (j>jStart)
				{
					if (BC(i, j - 1) == BC_NEUM) coefIJ += 0;
				}
				if (j<jEnd)
				{
					if (BC(i, j + 1) == BC_NEUM) coefIJ += 0;
				}
			}

			if (coefIJ == 0)
			{
				coefIJ = 1;
			}

			ipCSR.AssignValue(BC(i, j), BC(i, j), oneOverdx2*coefIJ / 2 + oneOverdt);
		}
	}
}

inline void MovingInterface::GenerateLinearSystem2(VectorND<double>& vectorB)
{
	Array2D<int>& BC = Surfactant.BC;

	termOld = term;

	SurfactantNormalTerm(Surfactant, levelSet, term);
	double oneOverdt = 1 / dt;
#pragma omp parallel for
	for (int j = Surfactant.jStart; j <= Surfactant.jEnd; j++)
	{
		for (int i = Surfactant.iStart; i <= Surfactant.iEnd; i++)
		{
			if (BC(i, j) < 0)
			{
				continue;
			}
			vectorB(BC(i, j)) = Surfactant(i, j)*oneOverdt + 0.5 *(Surfactant.dxxPhi(i, j) + Surfactant.dyyPhi(i, j))
				+ 1.5*term(i, j) - 0.5*termOld(i, j);
			if (i > Surfactant.jStart)
			{
				if (BC(i - 1, j) == BC_DIR)
				{
					Surfactant(i - 1, j) = ExactSurfactant(grid(i - 1, j).x, grid(i - 1, j).y, totalT);
					vectorB(BC(i, j)) += Surfactant(i - 1, j)*Surfactant.oneOverdx2 * 0.5;
				}
			}
			if (i < Surfactant.iEnd)
			{
				if (BC(i + 1, j) == BC_DIR)
				{
					Surfactant(i + 1, j) = ExactSurfactant(grid(i + 1, j).x, grid(i + 1, j).y, totalT);
					vectorB(BC(i, j)) += Surfactant(i + 1, j)*Surfactant.oneOverdx2 * 0.5;
				}
			}
			if (j > Surfactant.jStart)
			{
				if (BC(i, j - 1) == BC_DIR)
				{
					Surfactant(i, j - 1) = ExactSurfactant(grid(i, j - 1).x, grid(i, j - 1).y, totalT);
					vectorB(BC(i, j)) += Surfactant(i, j - 1)*Surfactant.oneOverdy2 * 0.5;
				}
			}
			if (j < Surfactant.jEnd)
			{
				if (BC(i, j + 1) == BC_DIR)
				{
					Surfactant(i, j + 1) = ExactSurfactant(grid(i, j + 1).x, grid(i, j + 1).y, totalT);
					vectorB(BC(i, j)) += Surfactant(i, j + 1)*Surfactant.oneOverdy2 * 0.5;
				}
			}
		}
	}

	Surfactant(Surfactant.iStart, Surfactant.jStart) = ExactSurfactant(grid(grid.iStart, grid.jStart).x, grid(grid.iStart, grid.jStart).y, totalT);
	Surfactant(Surfactant.iStart, Surfactant.jEnd) = ExactSurfactant(grid(grid.iStart, grid.jEnd).x, grid(grid.iStart, grid.jEnd).y, totalT);
	Surfactant(Surfactant.iEnd, Surfactant.jStart) = ExactSurfactant(grid(grid.iEnd, grid.jStart).x, grid(grid.iEnd, grid.jStart).y, totalT);
	Surfactant(Surfactant.iEnd, Surfactant.jEnd) = ExactSurfactant(grid(grid.iEnd, grid.jEnd).x, grid(grid.iEnd, grid.jEnd).y, totalT);
}

inline void MovingInterface::PlotSurfactant()
{
	Surfactant.Variable("Surfactant");
	MATLAB.Command("surf(X,Y,Surfactant);");
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
	Surfactant.Variable("Surfactant0");
	levelSet.tube.Variable("tube");
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	PlotLocalSurfactant();
	MATLAB.Command("IntSur0=sum(sum(Surfactant.*(Tube<=1)));");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);
	int extensionIter = (int)ceil(levelSet.gamma2 / (0.5*min(grid.dx, grid.dy)));;
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << "*******************************************************************" << endl;
		cout << "       Iteration " << to_string(i) << " : Start" << endl;
		totalT += dt;
		cout << "------------- Diffusion Start ------------- " << endl;

		LSurfactantDiffusion(i);

		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);

		cout << "------------ Diffusion End ------------- " << endl;

		PlotLocalSurfactant();
		MATLAB.Command("IntSur=sum(sum(Surfactant.*(Tube<=1)));");
		MATLAB.Command("loss = (IntSur0-IntSur)/IntSur0*100");
		MATLAB.Variable("i", i);
		MATLAB.Variable("totalT", totalT);
		str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', Loss(%)  :',num2str(loss)]);");
		MATLAB.Command(str.c_str());
	}
}

inline void MovingInterface::PlotLocalSurfactant()
{
	Surfactant.Variable("Surfactant");
	MATLAB.Command("SurTube1=Surfactant.*(Tube<=1);");
	levelSet.tube.Variable("Tube");
	bool isLookDown = true;
	if (isLookDown)
	{
		MATLAB.Command("surf(X,Y,SurTube1), hold on, contour(X,Y,Tube),h=colorbar,h.Limits=[0 max(max(SurTube1))],hold off,set(gca,'fontsize',20);axis([X(1) X(end) Y(1) Y(end)]),axis equal");
	}
	else
	{
		MATLAB.Command("surf(X,Y,SurTube1), hold on, contour(X,Y,Tube),h=colorbar,h.Limits=[0 max(max(SurTube1))],hold off,set(gca,'fontsize',20);");
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
	PlotLocalSurfactant();
	MATLAB.Command("IntSur0=sum(sum(Surfactant.*(Tube<=1)));");

	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);

	FD exact(grid);
	int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy)) * 2;
	int extensionIter = (int)ceil(levelSet.gamma2 / (0.5*min(grid.dx, grid.dy)));
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << "*******************************************************************" << endl;
		cout << "       Iteration " << to_string(i) << " : Start" << endl;
		totalT += dt;
		cout << "------------- Diffusion Start ------------- " << endl;

		LSurfactantDiffusion(i);

		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);

		cout << "------------ Diffusion End ------------- " << endl;

		//////////////////////////
		//// Moving Level Set ////
		//////////////////////////
		AdvectionMethod2D<double>::LLSPropagatingTVDRK3MACGrid(levelSet, U, V, dt);
		AdvectionMethod2D<double>::LLSReinitializationTVDRK3(levelSet, dt, reinitialIter);

		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);
		levelSet.UpdateInterface();
		levelSet.UpdateLLS();

		//// Surfactant Conservation Method
		ConserveSurfactantFactorBeta();

		PlotLocalSurfactant();
		MATLAB.Command("IntSur=sum(sum(Surfactant.*(Tube<=1)));");
		MATLAB.Command("loss = (IntSur0-IntSur)/IntSur0*100");
		MATLAB.Variable("i", i);
		MATLAB.Variable("totalT", totalT);
		MATLAB.Command("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', Loss(%)  :',num2str(loss)]);");

		cout << "       Iteration " << to_string(i) << " : End" << endl;
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
	if (iter == 1)
	{
		LGenerateLinearSystem1(Acsr);

		LOneStepSemiImplicit();
	}

	if (iter >= 2)
	{
		SurfactantTube2Extrapolation();
		
		LGenerateLinearSystem2(Acsr);

		LTwoStepSemiImplicit();
	}

}

inline void MovingInterface::LCountNonZero()
{
	int & num_all_full_cells = Surfactant.num_all_full_cells;
	int & nnz = Surfactant.nnz;
	int & iStart = Surfactant.iStart;
	int & iEnd = Surfactant.iEnd;
	int & jStart = Surfactant.jStart;
	int & jEnd = Surfactant.jEnd;

	VI leftIndex, rightIndex, bottomIndex, topIndex;
	int i, j, l, m, n, k;

	num_all_full_cells = 0;
	nnz = 0;

	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;

		num_all_full_cells++;
		nnz++;

		if (i>iStart)
		{
			leftIndex = VI(i - 1, j);
			l = levelSet.tubeIJ2K(leftIndex);
			levelSet.TubeIndex(l, m, n);
			if (levelSet.tube(m, n) == 1) nnz++;
		}
		
		if (i<iEnd)
		{
			rightIndex = VI(i + 1, j);
			l = levelSet.tubeIJ2K(rightIndex);
			levelSet.TubeIndex(l, m, n);
			if (levelSet.tube(m, n) == 1) nnz++;
		}
		
		if (j>jStart)
		{
			bottomIndex = VI(i, j - 1);
			l = levelSet.tubeIJ2K(bottomIndex);
			levelSet.TubeIndex(l, m, n);
			if (levelSet.tube(m, n) == 1) nnz++;
		}
		
		if (j<jEnd)
		{
			topIndex = VI(i, j + 1);
			l = levelSet.tubeIJ2K(topIndex);
			levelSet.TubeIndex(l, m, n);
			if (levelSet.tube(m, n) == 1) nnz++;
		}
		
	}
}

inline void MovingInterface::LOneStepSemiImplicit()
{
	//// Linear Equation
	LGenerateLinearSystem1(vectorB);
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

inline void MovingInterface::LTwoStepSemiImplicit()
{
	//// Linear Equation
	termOld = term;
	LGenerateLinearSystem2(vectorB);
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
	double velX, velY;
	int i, j;
#pragma omp parallel for private(i,j, normal, Hessian, curvature, velX, velY)
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
			velX = 0.5*(U(i + 1, j) + U(i, j));
			velY = 0.5*(V(i, j + 1) + V(i, j));
			//// Upwind WENO
			term(i, j) += -(AdvectionMethod2D<double>::Plus(velX)*wenoDxMinus(i, j) + AdvectionMethod2D<double>::Minus(velX)*wenoDxPlus(i, j));
			term(i, j) += -(AdvectionMethod2D<double>::Plus(velY)*wenoDyMinus(i, j) + AdvectionMethod2D<double>::Minus(velY)*wenoDyPlus(i, j));

			////
			term(i, j) += normal(0)*(gradientU(i, j).x*normal(0) + gradientU(i, j).y*normal(1))*ipField(i, j);
			term(i, j) += normal(1)*(gradientV(i, j).x*normal(0) + gradientV(i, j).y*normal(1))*ipField(i, j);
		}
	}
}

inline void MovingInterface::LGenerateLinearSystem1(CSR<double>& ipCSR)
{
	LCountNonZero();
	vectorB = VectorND<double>(Surfactant.num_all_full_cells);
	tempSur = VectorND<double>(Surfactant.num_all_full_cells);
	Acsr = CSR<double>(Surfactant.num_all_full_cells, Surfactant.nnz);

	int iStart = Surfactant.iStart, iEnd = Surfactant.iEnd, jStart = Surfactant.jStart, jEnd = Surfactant.jEnd;

	Array2D<int>& tube = levelSet.tube;
	Array2D<int>& tube1 = levelSet.tube1;

	double dx = Surfactant.dx, dy = Surfactant.dy;
	double dx2 = dx*dx, dy2 = dy*dy, dxdy = dx*dy;
	double oneOverdx = 1 / dx, oneOverdx2 = 1 / dx2;
	double oneOverdy = 1 / dy, oneOverdy2 = 1 / dy2;
	double oneOverdt = 1 / dt;
	double coefIJ;
	double oneOverPe = 1 / Pe;

	VI leftIndex, rightIndex, bottomIndex, topIndex;
	int i, j, l, m, n, k1, l1;
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		coefIJ = 0;
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;
		k1 = k - 1;

		if (i>iStart)
		{
			leftIndex = VI(i - 1, j);
			l = levelSet.tubeIJ2K(leftIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 1)
			{
				l = tube1(m, n);
				l1 = l - 1;
				ipCSR.AssignValue(k1, l1, -oneOverdx2*oneOverPe);
			}
			if (tube(m,n)==1 || tube(m,n)==2)
			{
				coefIJ++;
			}
		}
		
		if (i<iEnd)
		{
			rightIndex = VI(i + 1, j);
			l = levelSet.tubeIJ2K(rightIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 1)
			{
				l = tube1(m, n);
				l1 = l - 1;
				ipCSR.AssignValue(k1, l1, -oneOverdx2*oneOverPe);
			}
			if (tube(m, n) == 1 || tube(m, n) == 2)
			{
				coefIJ++;
			}
		}
		
		if (j>jStart)
		{
			bottomIndex = VI(i, j - 1);
			l = levelSet.tubeIJ2K(bottomIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 1)
			{
				l = tube1(m, n);
				l1 = l - 1;
				ipCSR.AssignValue(k1, l1, -oneOverdy2*oneOverPe);
			}
			if (tube(m, n) == 1 || tube(m, n) == 2)
			{
				coefIJ++;
			}
		}
		
		if (j<jEnd)
		{
			topIndex = VI(i, j + 1);
			l = levelSet.tubeIJ2K(topIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 1)
			{
				l = tube1(m, n);
				l1 = l - 1;
				ipCSR.AssignValue(k1, l1, -oneOverdy2*oneOverPe);
			}
			if (tube(m, n) == 1 || tube(m, n) == 2)
			{
				coefIJ++;
			}
		}

		if (coefIJ==0)
		{
			coefIJ = 1;
		}
		ipCSR.AssignValue(k1, k1, oneOverdx2 * coefIJ*oneOverPe + oneOverdt);
	}
}

inline void MovingInterface::LGenerateLinearSystem1(VectorND<double>& vectorB)
{
	int iStart = Surfactant.iStart, iEnd = Surfactant.iEnd, jStart = Surfactant.jStart, jEnd = Surfactant.jEnd;
	double oneOverdt = 1 / dt;
	double oneOverdx2 = grid.oneOverdx2;
	double oneOverdy2 = grid.oneOverdy2;

	Array2D<int>& tube = levelSet.tube;

	double* bVal(vectorB.values);

	LSurfactantNormalTerm(Surfactant, levelSet, term);
	VI leftIndex, rightIndex, bottomIndex, topIndex;
	double oneOverPe = 1 / Pe;
	int i, j, l, m, n;
#pragma omp parallel for private(i, j, l, m, n, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;
		bVal[k - 1] = (Surfactant(i, j)*oneOverdt + term(i, j)*oneOverPe);

		if (i>iStart)
		{
			leftIndex = VI(i - 1, j);
			l = levelSet.tubeIJ2K(leftIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 2 || (tube(m, n) == 1 && i == iStart + 1))
			{
				bVal[k - 1] += oneOverdx2*Surfactant(m, n) * oneOverPe;
			}
		}
		
		if (i<iEnd)
		{
			rightIndex = VI(i + 1, j);
			l = levelSet.tubeIJ2K(rightIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 2 || (tube(m, n) == 1 && i == iEnd - 1))
			{
				bVal[k - 1] += oneOverdx2*Surfactant(m, n) * oneOverPe;
			}
		}
		
		if (j>jStart)
		{
			bottomIndex = VI(i, j - 1);
			l = levelSet.tubeIJ2K(bottomIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 2 || (tube(m, n) == 1 && j == jStart + 1))
			{
				bVal[k - 1] += oneOverdy2*Surfactant(m, n) * oneOverPe;
			}
		}
		
		if (j<jEnd)
		{
			topIndex = VI(i, j + 1);
			l = levelSet.tubeIJ2K(topIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 2 || (tube(m, n) == 1 && j == jEnd - 1))
			{
				bVal[k - 1] += oneOverdy2*Surfactant(m, n) * oneOverPe;
			}
		}
	}
}

inline void MovingInterface::LGenerateLinearSystem2(CSR<double>& ipCSR)
{
	LCountNonZero();
	vectorB = VectorND<double>(Surfactant.num_all_full_cells);
	tempSur = VectorND<double>(Surfactant.num_all_full_cells);
	Acsr = CSR<double>(Surfactant.num_all_full_cells, Surfactant.nnz);

	int iStart = Surfactant.iStart, iEnd = Surfactant.iEnd, jStart = Surfactant.jStart, jEnd = Surfactant.jEnd;

	Array2D<int>& tube = levelSet.tube;
	Array2D<int>& tube1 = levelSet.tube1;

	double dx = Surfactant.dx, dy = Surfactant.dy;
	double dx2 = dx*dx, dy2 = dy*dy, dxdy = dx*dy;
	double oneOverdx = 1 / dx, oneOverdx2 = 1 / dx2;
	double oneOverdy = 1 / dy, oneOverdy2 = 1 / dy2;
	double oneOverdt = 1 / dt;
	double coefIJ;
	double oneOverPe = 1 / Pe;

	VI leftIndex, rightIndex, bottomIndex, topIndex;
	int i, j, l, m, n, k1, l1;
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		coefIJ = 0;
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;
		k1 = k - 1;

		if (i>iStart)
		{
			leftIndex = VI(i - 1, j);
			l = levelSet.tubeIJ2K(leftIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 1)
			{
				l = tube1(m, n);
				l1 = l - 1;
				ipCSR.AssignValue(k1, l1, -oneOverdx2 * oneOverPe * 0.5);
			}
			if (tube(m, n) == 1 || tube(m, n) == 2)
			{
				coefIJ ++;
			}
		}
		
		if (i<iEnd)
		{
			rightIndex = VI(i + 1, j);
			l = levelSet.tubeIJ2K(rightIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 1)
			{
				l = tube1(m, n);
				l1 = l - 1;
				ipCSR.AssignValue(k1, l1, -oneOverdx2 * oneOverPe * 0.5);
			}
			if (tube(m, n) == 1 || tube(m, n) == 2)
			{
				coefIJ++;
			}
		}

		if (j>jStart)
		{
			bottomIndex = VI(i, j - 1);
			l = levelSet.tubeIJ2K(bottomIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 1)
			{
				l = tube1(m, n);
				l1 = l - 1;
				ipCSR.AssignValue(k1, l1, -oneOverdy2 * oneOverPe * 0.5);
			}
			if (tube(m, n) == 1 || tube(m, n) == 2)
			{
				coefIJ++;
			}
		}

		if (j<jEnd)
		{
			topIndex = VI(i, j + 1);
			l = levelSet.tubeIJ2K(topIndex);
			levelSet.TubeIndex(l, m, n);
			if (tube(m, n) == 1)
			{
				l = tube1(m, n);
				l1 = l - 1;
				ipCSR.AssignValue(k1, l1, -oneOverdy2 * oneOverPe * 0.5);
			}
			if (tube(m, n) == 1 || tube(m, n) == 2)
			{
				coefIJ++;
			}
		}

		if (coefIJ == 0)
		{
			coefIJ = 1;
		}
		ipCSR.AssignValue(k1, k1, oneOverdx2 * coefIJ * 0.5 * oneOverPe + oneOverdt);

	}
}

inline void MovingInterface::LGenerateLinearSystem2(VectorND<double>& vectorB)
{
	int iStart = Surfactant.iStart, iEnd = Surfactant.iEnd, jStart = Surfactant.jStart, jEnd = Surfactant.jEnd;
	LSurfactantNormalTerm(Surfactant, levelSet, term);
	VI leftIndex, rightIndex, bottomIndex, topIndex;
	double oneOverdt = 1 / dt;
	double oneOverdx2 = grid.oneOverdx2;
	double oneOverdy2 = grid.oneOverdy2;
	
	double* bVal(vectorB.values);
	double oneOverPe = 1 / Pe;
	int i, j, l, m, n;
#pragma omp parallel for private(i, j, l, m, n, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int k = 1; k <= levelSet.numTube1; k++)
	{
		i = levelSet.tube1Index(k).i;
		j = levelSet.tube1Index(k).j;
		bVal[k - 1] = Surfactant(i, j) * oneOverdt + 0.5 * oneOverPe *((Surfactant.dxxPhi(i, j) + Surfactant.dyyPhi(i, j))
			+ 3 * term(i, j) - 1 * termOld(i, j));

		if (i>iStart)
		{
			leftIndex = VI(i - 1, j);
			l = levelSet.tubeIJ2K(leftIndex);
			levelSet.TubeIndex(l, m, n);
			if (levelSet.tube(m, n) == 2 || (levelSet.tube(m, n) == 1 && i == iStart + 1))
			{
				bVal[k - 1] += oneOverdx2*Surfactant(m, n) * 0.5 * oneOverPe;
			}
		}

		if (i<iEnd)
		{
			rightIndex = VI(i + 1, j);
			l = levelSet.tubeIJ2K(rightIndex);
			levelSet.TubeIndex(l, m, n);
			if (levelSet.tube(m, n) == 2 || (levelSet.tube(m, n) == 1 && i == iEnd - 1))
			{
				bVal[k - 1] += oneOverdx2*Surfactant(m, n) * 0.5 * oneOverPe;
			}
		}

		if (j>jStart)
		{
			bottomIndex = VI(i, j - 1);
			l = levelSet.tubeIJ2K(bottomIndex);
			levelSet.TubeIndex(l, m, n);
			if (levelSet.tube(m, n) == 2 || (levelSet.tube(m, n) == 1 && j == jStart + 1))
			{
				bVal[k - 1] += oneOverdy2*Surfactant(m, n) * 0.5 * oneOverPe;
			}
		}

		if (j<jEnd)
		{
			topIndex = VI(i, j + 1);
			l = levelSet.tubeIJ2K(topIndex);
			levelSet.TubeIndex(l, m, n);
			if (levelSet.tube(m, n) == 2 || (levelSet.tube(m, n) == 1 && j == jEnd - 1))
			{
				bVal[k - 1] += oneOverdy2*Surfactant(m, n) * 0.5 * oneOverPe;
			}
		}
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

