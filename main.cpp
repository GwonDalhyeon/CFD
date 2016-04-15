#include "PoissonSolver.h"

#include "AdvectionMethod2D.h"
#include "SurfaceReconstruction.h"
#include "LevelSetReinitializationProblem.h"

#include "PointIntegralMethod.h"

#include "LevelSetAdvectionProblem.h"
#include "EulerianFluidSolver.h"

#include "VortexSheet.h"
#include "ToMATLAB.h"


void main()
{
	MATLAB.Command("clc; clear all; close all;");
	MATLAB.Command("workspace");
	
	int numP = 15;
	VectorND<Vector2D<double>>pointCloud(1, numP);
	Vector2D<double> tempVector;
	srand(time(NULL));
	//for (int i = 0; i <= 4; i++)
	//{
	//	pointCloud(3 * i + 1) = Vector2D<double>(i, 0);
	//	pointCloud(3 * i + 2) = Vector2D<double>(i, 1);
	//	pointCloud(3 * i + 3) = Vector2D<double>(i, 2);

	//}

	for (int i = 1; i <= numP; i++)
	{
		tempVector = Vector2D<double>(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX)) - 0.5;
		pointCloud(i) = (tempVector);
		for (int j = i; j >= 2; j--)
		{
			if (pointCloud(j).x < pointCloud(j-1).x)
			{
				tempVector = pointCloud(j);
				pointCloud(j) = pointCloud(j - 1);
				pointCloud(j - 1) = tempVector;
			}
			else
			{
				continue;
			}
		}
		
	}

	//VectorND<Polygon2D> DelTri = DelaunayTriangulate(pointCloud);
	DelaunayTriangulate(pointCloud);
	
//	int omp_get_max_threads(void);
//	void omp_set_num_threads(int);
//	cout << omp_get_max_threads() << endl;
//
//	omp_set_num_threads(omp_get_max_threads());
//
//	cout << omp_get_thread_num() << endl;

	//PointIntegralMethod<double> PIM;
	//PIM.PointIntegralMethodnSolver(1);

	//VortexSheet vortex;
	//vortex.VortexSolver(1);

	//LevelSetAdvection levelSet;
	//levelSet.advectionSolver(4, 0.1);

	//EulerianFluidSolver2D Euler;
	//Euler.FluidSolver(1, 1);


	//SurfaceReconst<double> surface;
	//surface.SurfaceReconstructionSplitBregman(1);


	//SurfaceReconst<double> surface;
	//surface.SurfaceReconstructionSolver(1);
	
	//Reinitialzation reinitial;
	//reinitial.ReinitializationSolver(1);
	
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

	system("pause");
}