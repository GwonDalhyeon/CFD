#include "PoissonSolver.h"

#include "LevelSetAdvectionProblem.h"
#include "LevelSetReinitializationProblem.h"
#include "LocalLevelSetProb.h"

#include "FluidSolver2D.h"
#include "EulerianFluidSolver.h"
#include "VortexSheet.h"
#include "EulerianMovingInterface.h"
#include "Incompressible2PhasewithSurfactantProb.h"
#include "CoalescingDropProblem.h"

#include "SurfaceReconstruction.h"
#include "BregmaMethodSolver.h"
#include "PointIntegralMethod.h"

#include "EulerEquation1D.h"

#include "mgmres.h"

void main()
{
	//	int omp_get_max_threads(void);
	//	void omp_set_num_threads(int);
	//	cout << omp_get_max_threads() << endl;
	//
	//	omp_set_num_threads(omp_get_max_threads());
	//
	//	cout << omp_get_thread_num() << endl;
	cout << "DALHYEON's Fluid Solver" << endl;
	cout << endl;
	MATLAB.Command("clc; clear all; close all;");
	MATLAB.Command("workspace");


	bool selection = false;
	int problem;
	int example;

	if (selection)
	{
		cout << "  1 : Not Yet" << endl;
		//cout << " 1 : POINT INTEGRAL METHOD " << endl;
		cout << "  2 : VORTEX SHEET" << endl;
		cout << "  3 : LEVEL SET ADVECTION" << endl;
		cout << "  4 : NAVIER STOKES EQUATION" << endl;
		cout << "  5 : Not Yet" << endl;
		//cout << " 5 : SURFACE RECONSTRUCTION USING SPLIT BREGMAN" << endl;
		cout << "  6 : SURFACE RECONSTRUCTION" << endl;
		cout << "  7 : LEVEL SET REINITIALIZATION" << endl;
		cout << "  8 : POISSON EQUATION, BOUNDARY CAPTURING" << endl;
		cout << "  9 : Not Yet" << endl;
		cout << " 10 : Not Yet" << endl;
		cout << " 11 : Not Yet" << endl;
		cout << " 12 : Eulerian Moving Interface. (By JJ Xu, HK Zhao)" << endl;
		cout << " 13 : Local Level Set Problem Solver" << endl;
		cout << " 14 : Surfactant + NS equation" << endl;
		cout << " 15 : Surfactant + Bubble Coalescing" << endl;
		cout << endl;
		cout << "Choose Problem Number! : ";
		cin >> problem;
		if (problem == 1)
		{
			cout << endl;
			cout << " It does not work. To be Continue... " << endl;
			system("pause");
			return;
		}
		if (problem == 2)
		{
			cout << endl;
			cout << " 1 : Vortex Sheet in 2D " << endl;
			cout << " 2 : Vortex Sheet Dipole" << endl;

		}
		if (problem == 3)
		{
			cout << endl;
			cout << " 1 : Level set advection Test - Rotating a circle" << endl;
			cout << " 2 : Level set advection Test - Spiral" << endl;
			cout << " 3 : Level set advection Test - Seven-point Star" << endl;
			cout << " 4 : Level set advection Test - Unit square" << endl;
			cout << " 5 : Level set advection Test - Unit sircle" << endl;
			cout << " 6 : Level set advection Test - Diamond" << endl;
			cout << " 7 : Level set advection Test - Concave" << endl;

		}
		if (problem == 4)
		{
			cout << endl;
			cout << " 1 : Cavity Flow" << endl;
			cout << " 2 : Tube" << endl;
			cout << " 3 : No !!!!" << endl;
			cout << " 4 : Not Yet..." << endl;
		}
		if (problem == 5)
		{
			cout << endl;
			cout << " It does not work. To be Continue... " << endl;
			system("pause");
			return;
		}
		if (problem == 6)
		{
			cout << endl;
			cout << " 1 : Surface reconstruction : One circles" << endl;
			cout << " 2 : Surface reconstruction : One circles" << endl;
			cout << " 3 : Surface reconstruction : Two circles" << endl;
			cout << " 4 : Surface reconstruction : Two circles with outlier, Initial Level Set is far from real shape" << endl;
			cout << " 5 : Surface reconstruction : Two circles with outlier, Good Initial Level Set" << endl;
			cout << " 6 : Surface reconstruction : Two circles with outlier" << endl;
		}
		if (problem == 7)
		{
			cout << endl;
			cout << " 1 : A circle with center at the origen and radius 1" << endl;
			cout << " 2 : two circles intersect each other" << endl;
			cout << " 3 : Not Yet..." << endl;
			cout << " 4 : Discontinue Level Set" << endl;
		}
		if (problem == 8)
		{
			cout << endl;
			cout << " 2 : Unique example." << endl;
		}
		if (problem == 9)
		{
			cout << endl;
			cout << " It does not work. To be Continue... " << endl;
			system("pause");
			return;
		}
		if (problem == 10)
		{
			cout << endl;
			cout << " It does not work. To be Continue... " << endl;
			system("pause");
			return;
		}
		if (problem == 11)
		{
			cout << endl;
			cout << " It does not work. To be Continue... " << endl;
			system("pause");
			return;
		}
		if (problem == 12)
		{
			cout << endl;
			cout << " 1 : Surfactant Diffusion" << endl;
			cout << " 2 : Surfactant Diffusion with Velocity" << endl;
			cout << " 3 : Locally Surfactant Diffusion" << endl;
			cout << " 4 : Locally Surfactant Diffusion with Velocity 1" << endl;
			cout << " 5 : Locally Surfactant Diffusion with Velocity 2" << endl;
			cout << " 6 : Locally Surfactant Diffusion with Velocity 3" << endl;
		}
		if (problem == 13)
		{
			cout << endl;
			cout << " 1 : Local Local Level set advection Test - Rotating a circle" << endl;
			cout << " 2 : Local Level set advection Test - Unit circle" << endl;
			cout << " 3 : Quantity Extension using Local LS - Unit circle" << endl;
		}
		if (problem == 14)
		{
			cout << endl;
			cout << " 1 : Uniform Surfactant + NS equation" << endl;
		}
		if (problem == 15)
		{
			cout << endl;
			cout << " 1 : Surfactant + Bubble Coalescing " << endl;
		}

		cout << endl;
		cout << "Choose Example Number! : ";
		cin >> example;
	}
	else
	{
		problem = 14;
		example = 1;
	}

	if (problem == 0) //// GMRES test :: example 1,2,3,4
	{
		switch (example)
		{
		case(1):
			GMRESLib::test01();
			break;
		case(2):
			GMRESLib::test02();
			break;
		case(3):
			GMRESLib::test03();
			break;
		case(4):
			GMRESLib::test04();
			break;
		}
	}
	else if (problem == 1) // POINT INTEGRAL METHOD : EXAMPLE 1
	{
		PointIntegralMethod<double> PIM;
		PIM.PointIntegralMethodnSolver(example);
	}
	else if (problem == 2) // VORTEX SHEET : example 1,2
	{
		VortexSheet vortex;
		vortex.VortexSolver(example);
	}
	else if (problem == 3) // LEVEL SET ADVECTION : EXAMPLE 1,2,3,4,5,6,7
	{
		LevelSetAdvection levelSet;
		levelSet.AdvectionSolver(example);
	}
	else if (problem == 4) // NAVIER STOKES EQUATION SOLVER : EXAMPLE 1 ??
	{
		// Single Phase
		if (example == 1 || example == 2)
		{
			FluidSolver2D Fluid;
			Fluid.Solver(example);
		}
		// Multi Phase
		if (example == 4 || example == 5 || example == 6)
		{
			FluidSolver2D Fluid;
			Fluid.Solver(example);
		}
	}
	else if (problem == 5) // SURFACE RECONSTRUCTION USING SPLIT BREGMAN : EXAMPLE 1 ~ 6 ??
	{
		SurfaceReconst<double> surface;
		surface.SurfaceReconstructionSplitBregman(example);
	}
	else if (problem == 6) // SURFACE RECONSTRUCTION : EXAMPLE 1 ~ 6
	{
		SurfaceReconst<double> surface;
		surface.SurfaceReconstructionSolver(example);
	}
	else if (problem == 7) // LEVEL SET REINITIALIZATION : EXAMPLE 
	{
		Reinitialzation reinitial;
		reinitial.ReinitializationSolver(example);
	}
	else if (problem == 8) // POISSON EQUATION, BOUNDARY CAPTURING : EXAMPLE 2
	{
		PoissonSolver testPoisson2d;
		testPoisson2d.SolvePoissonJumpCondi(example);
	}
	else if (problem == 9)
	{
		//GridInfo testGrid1d(0.0, 1.0, 101);
		//LaplaceEquationSolver testLaplace = LaplaceEquationSolver(testGrid1d);
		//testLaplace.solveLaplaceEquationJumpCondi(example);
	}
	else if (problem == 10)
	{
		//GridInfo testGrid1d(0.0, 1.0, 101);
		//PoissonEquationSolver testPoisson1d = PoissonEquationSolver(testGrid1d);
		//testPoisson1d.solvePoissonEquationJumpCondi(example);
	}
	else if (problem == 11)
	{
		//GridInfo testGrid2d(0.0, 1.0, 101, 0.0, 1.0, 101);
		//PoissonEquationSolver testPoisson2d = PoissonEquationSolver(testGrid2d);
		//testPoisson2d.solvePoissonEquationJumpCondi(example);
	}
	else if (problem == 12)  // Eulerian Moving Interface.(JJ Xu, HK Zhao) : Example 1
	{
		FluidSolver2D Fluid;
		MovingInterface Eulerian(Fluid);
		if (example <= 2) // Whole domain
		{
			Eulerian.SurfactantDiffusionSolver(example);
		}
		else if (example == 3) // Only Near Interface
		{
			Eulerian.LSurfactantDiffusionSolver(example);
		}
		else if (example == 4 || example == 5 || example == 6)
		{
			Eulerian.EulerianMovingInterfaceSolver(example);
		}
	}
	else if (problem == 13) // Local Level Set Problem Solver : Example 1,2
	{
		LocalLevelSetAdvection LLS;
		if (example == 1 || example == 2)
		{
			LLS.AdvectionSolver(example);
		}
		else if (example == 3)
		{
			LLS.QuantityExtensionSolver(example);
		}
	}
	else if (problem == 14)
	{
		FluidSolver2D Fluid;
		MovingInterface Surfactant(Fluid);
		InsolubleSurfactant ContinuumMethod(Fluid, Surfactant);
		ContinuumMethod.ContinuumMethodWithSurfactantSolver(example);
	}
	else if (problem == 15)
	{
		FluidSolver2D Fluid;
		MovingInterface Surfactant(Fluid);
		Coalescence Coalescing(Fluid, Surfactant);
		Coalescing.DropCoalescenceSolver(example);
	}
	//system("pause");
	return;
}
