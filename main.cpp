#include "PoissonSolver.h"

#include "LevelSetAdvectionProblem.h"
#include "LevelSetReinitializationProblem.h"
#include "LocalLevelSetProb.h"

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

	MATLAB.Command("clc; clear all; close all;");
	MATLAB.Command("workspace");

	int problem = 15;
	int example = 1;
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
		EulerianFluidSolver2D Euler;
		if (example == 1 || example == 4)
		{
			Euler.FluidSolver(example);
		}
		else if (example ==2)
		{
			Euler.NSSolver2ndOrder(example);
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
		EulerianFluidSolver2D Fluid;
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
		EulerianFluidSolver2D Fluid;
		MovingInterface Surfactant(Fluid);
		InsolubleSurfactant ContinuumMethod(Fluid, Surfactant);
		ContinuumMethod.ContinuumMethodWithSurfactantSolver(example);
		
	}
	else if (problem == 15)
	{
		EulerianFluidSolver2D Fluid;
		MovingInterface Surfactant(Fluid);
		CoalescingDrop Coalescing(Fluid, Surfactant);

	}
	system("pause");
}
