#include "PoissonSolver.h"

#include "AdvectionMethod2D.h"
#include "SurfaceReconstruction.h"
#include "LevelSetReinitializationProblem.h"

#include "LevelSetAdvectionProblem.h"
#include "EulerianFluidSolver.h"

#include "VortexSheet.h"
#include "ToMATLAB.h"
void main()
{
	MATLAB.Command("clear all;close all;");
	MATLAB.Command("workspace");

//	int omp_get_max_threads(void);
//	void omp_set_num_threads(int);
//	cout << omp_get_max_threads() << endl;
//
//	omp_set_num_threads(omp_get_max_threads());
//
//	cout << omp_get_thread_num() << endl;


	//testVN.Variable("testVN");
	
	//VortexSheet vortex;
	//vortex.VortexSolver(1, 0.1);

	LevelSetAdvection levelSet;
	levelSet.advectionSolver(4, false, false, true, 0.1);

	//EulerianFluidSolver2D Euler;
	//Euler.FluidSolver(1, 1);

	//LevelSetAdvection levelSet;
	//levelSet.SurfaceReconstructionSplitBregman(7, false, false, true);

	//SurfaceReconst<double> surface;
	//surface.surfaceReconstructionSolver(4);
	
	//Reinitialzation reinitial;
	//reinitial. reinitializationSolver(4);
	
	//PoissonSolver testPoisson2d;
	//testPoisson2d.solvePoissonJumpCondi(2, testGrid2d);
	
	//GridInfo testGrid1d(0.0, 1.0, 101);
	//LaplaceEquationSolver testLaplace = LaplaceEquationSolver(testGrid1d);
	//testLaplace.solveLaplaceEquationJumpCondi(1);

	//GridInfo testGrid1d(0.0, 1.0, 101);
	//PoissonEquationSolver testPoisson1d = PoissonEquationSolver(testGrid1d);
	//testPoisson1d.solvePoissonEquationJumpCondi(1);

	//GridInfo testGrid2d(0.0, 1.0, 101, 0.0, 1.0, 101);
	//PoissonEquationSolver testPoisson2d = PoissonEquationSolver(testGrid2d);
	//testPoisson2d.solvePoissonEquationJumpCondi(2);
}