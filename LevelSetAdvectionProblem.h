#pragma once


//double AdvectionMethod2D<double>::alpha;
//#ifndef LevelSetAdvectionProblem_H
//#define LevelSetAdvectionProblem_H
#include "CommonDef.h"
#include "Grid2D.h"
#include "Field2D.h"
#include "LevelSet2D.h"
#include "AdvectionMethod2D.h"

class LevelSetAdvection
{
public:
	Grid2D grid;
	LevelSet2D levelSet;

	double cflCondition;
	double dt;

	bool isPropagation, isReinitialization, isSurfaceReconstruction;

	// Level set propagation velocity.
	bool isVelocity;
	bool needReinitial;
	Field2D<double> velocityX;
	Field2D<double> velocityY;

	int reinitialIter;
	int propaMaxIteration;
	int propaWriteIter;


	// Zero level set point.
	int givenPointNum;
	VectorND<Vector2D<double>> givenPoint;


	// Reinitialization Test Variable
	LevelSet2D exactLevelSet;

	int reinitialMaxIteration;
	int reinitialWriteIter;


	// Surface reconstruction Variable
	Field2D<double> distance;
	Field2D<double> reconstructionVelocity;

	int reconstMaxIteration;
	int reconstWriteIter;

	double LpNorm;
	double distanceThreshold;
	double curvatureThreshold;




	LevelSetAdvection();
	~LevelSetAdvection();




	void initialCondition(const int & example, const bool & propa, const bool & reinitial, const bool & surfReconst);
	void advectionSolver(const int & example, const bool & propa, const bool & reinitial, const bool & surfReconst, const double & cfl);




	// Propagation.
	void propaInitialCondition(const int& example);

	// Reinitialization
	void reinitializationInitialCondition(const int& example);

	// Surface reconstruction
	void surfReconstInitialCondition(const int& example);

	void computeVelocity();
	double computeIntegralTerm();
	void levelSetPropagatingTVDRK3();
	bool stoppingCriterion();


	// Distance functions.
	double distance2Data(const int&i, const int& j);
	void exactDistance();
	void sweepingDistance();

	// Adaptive time step functions.
	double adaptiveTimeStep();
	double adaptiveTimeStep(const Field2D<double>& velocity1);
	double adaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2);

	void outputResult(const int& iter);


private:

};

//#endif // !LevelSetAdvectionProblem_H




LevelSetAdvection::LevelSetAdvection()
{
}

LevelSetAdvection::~LevelSetAdvection()
{
}

inline void LevelSetAdvection::initialCondition(const int& example, const bool & propa, const bool & reinitial, const bool & surfReconst)
{
	isPropagation = propa;
	isReinitialization = reinitial;
	isSurfaceReconstruction = surfReconst;

	if (isPropagation)
	{
		propaInitialCondition(example);
	}
	else if (isReinitialization)
	{
		reinitializationInitialCondition(example);
	}
	else if (isSurfaceReconstruction)
	{
		surfReconstInitialCondition(example);
	}
	else
	{
		cout << "No example." << endl;
	}
}


inline void LevelSetAdvection::advectionSolver(const int & example, const bool & propa, const bool & reinitial, const bool & surfReconst, const double & cfl)
{
	cflCondition = cfl;
	initialCondition(example, propa, reinitial, surfReconst);

	// Reconstruction
	if (isSurfaceReconstruction)
	{
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy); // delta function parameter.
		outputResult(0);

		int i = 0;
		while (!stoppingCriterion() && i < reconstMaxIteration)
		{
			i++;
			cout << "Surface reconstruction : " << i << endl;

			computeVelocity();

			dt = adaptiveTimeStep(reconstructionVelocity);
			levelSetPropagatingTVDRK3();

			if (needReinitial)
			{
				dt = adaptiveTimeStep();
				for (int j = 0; j < reinitialIter; j++)
				{
					cout << "Reinitialization : " << i << "-" << j + 1 << endl;
					AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
				}
			}

			if (i %reconstWriteIter == 0)
			{
				outputResult(i);
			}
			cout << endl;
		}
		outputResult(i);
	}

	// Reinitialization
	if (isReinitialization)
	{
		outputResult(0);
		for (int i = 1; i <= reinitialMaxIteration; i++)
		{

			cout << "Reinitialization : " << i << endl;
			dt = adaptiveTimeStep();
			AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
			if (i%reinitialWriteIter == 0)
			{
				outputResult(i);
			}
		}
	}


	// Propagation.
	if (isPropagation)
	{
		outputResult(0);

		for (int i = 1; i <= propaMaxIteration; i++)
		{
			cout << "Level set advection : " << i << endl;

			if (isVelocity)
			{
				dt = adaptiveTimeStep(velocityX, velocityY);
				AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, velocityX, velocityY, dt);
			}
			else
			{
				dt = adaptiveTimeStep();
				AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, dt);
			}


			if (needReinitial)
			{
				for (int j = 0; j < reinitialIter; j++)
				{
					cout << "Reinitialization : " << i << "-" << j + 1 << endl;
					dt = adaptiveTimeStep();
					AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
				}
			}

			if (i%propaWriteIter == 0)
			{
				outputResult(i);
			}
		}
	}
}


inline void LevelSetAdvection::propaInitialCondition(const int & example)
{
	if (example == 1)
	{
		isVelocity = true;
		needReinitial = false;

		cout << "Level set advection Test - Rotating a circle" << endl;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);

		levelSet = LevelSet2D(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt((grid(i, j).x - 0.5)*(grid(i, j).x - 0.5) + grid(i, j).y*grid(i, j).y) - 0.25;
			}
		}

		velocityX = Field2D<double>(grid);
		velocityY = Field2D<double>(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				velocityX(i, j) = -grid(i, j).y;
				velocityY(i, j) = grid(i, j).x;
			}
		}

		isVelocity = true;

		//dt = grid.dx*grid.dy;
		propaMaxIteration = 10000;
		propaWriteIter = 100;
	}
	else if (example == 2)
	{
		cout << "Level set advection Test - Spiral" << endl;

		isVelocity = false;
		needReinitial = false;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);


		levelSet = LevelSet2D(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt((grid(i, j).x - 0.5)*(grid(i, j).x - 0.5) + grid(i, j).y*grid(i, j).y) - 0.25;
			}
		}


		//dt = grid.dx*grid.dy;
		propaMaxIteration = 1000;
		propaWriteIter = 10;
	}
	else if (example == 3)
	{
		cout << "Level set advection Test - Seven-point Star" << endl;

		isVelocity = true;
		needReinitial = false;

		grid = Grid2D(-0.25, 0.25, 101, -0.25, 0.25, 101);

		givenPointNum = 100;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum);

		double s;
#pragma omp parallel for private (s)
		for (int i = 0; i < givenPointNum; i++)
		{
			s = double(i) / double(givenPointNum);
			givenPoint(i) = (0.1 + 0.065*sin(7 * 2 * PI*s))*Vector2D<double>(cos(2 * PI*s), sin(2 * PI*s));
		}

		levelSet = LevelSet2D(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{

				levelSet(i, j) = distance2Data(i, j);
			}
		}

		isVelocity = false;

		//dt = grid.dx*grid.dy;
		propaMaxIteration = 1000;
		propaWriteIter = 10;
	}
}

inline void LevelSetAdvection::reinitializationInitialCondition(const int & example)
{
	grid = Grid2D(-2, 2, 101, -2, 2, 101);

	//dt = grid.dx*grid.dy;
	reinitialMaxIteration = 1000;
	reinitialWriteIter = 10;

	exactLevelSet = LevelSet2D(grid);
	levelSet = LevelSet2D(grid);

	isVelocity = false;
	needReinitial = false;

	double a = 0.7;
	double r = 1.0;

	if (example == 1) //// Exmple1. A circle with center at the origen and radius 1.
	{
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.jEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = (grid(i, j).magnitude() - 1.0)*((grid(i, j) - 1.0).magnitude() + 0.1);
				exactLevelSet(i, j) = grid(i, j).magnitude() - 1.0;
			}
		}
	}
	else if (example == 2) //// Exmple2. Two circles of radius r are placed at (+-a,0)  and a sqruar on the plane. Let 0<a<r, sh that the two circles intersect each other.
	{

		double temp1, temp2, temp3;
#pragma omp parallel for private(temp1,temp2,temp3)
		for (int i = grid.iStart; i <= grid.jEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{

				temp1 = (a - grid(i, j).x) / sqrt((a - grid(i, j).x)*(a - grid(i, j).x) + grid(i, j).y * grid(i, j).y);
				temp2 = (a + grid(i, j).x) / sqrt((a + grid(i, j).x)*(a + grid(i, j).x) + grid(i, j).y * grid(i, j).y);
				if (temp1 >= a / r && temp2 >= a / r)
				{
					temp3 = min(grid(i, j).x * grid(i, j).x + (grid(i, j).y + sqrt(r*r - a*a))*(grid(i, j).y + sqrt(r*r - a*a)), grid(i, j).x * grid(i, j).x + (grid(i, j).y - sqrt(r*r - a*a))*(grid(i, j).y - sqrt(r*r - a*a)));
					temp3 = sqrt(temp3);
					levelSet(i, j) = temp3 * ((grid(i, j) - 1).magnitude2() + 0.1);
					exactLevelSet(i, j) = temp3;
				}
				else
				{
					temp3 = min(sqrt((grid(i, j).x + a)*(grid(i, j).x + a) + grid(i, j).y * grid(i, j).y), sqrt((grid(i, j).x - a)*(grid(i, j).x - a) + grid(i, j).y * grid(i, j).y)) - r;
					levelSet(i, j) = temp3 * ((grid(i, j) - 1).magnitude2() + 0.1);
					exactLevelSet(i, j) = temp3;
				}
			}
		}

	}
	else if (example == 3) ////Exmple3.Two circles of radius r are placed at(+-a, 0) on the plane.Let 0<a<r, sh that the two circles intersect each other.
	{
		//double temp1, temp2, temp3;
		//for (int i = grid.iStart; i <= grid.jEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		temp1 = (a - grid(i, j).x) / sqrt((a - grid(i, j).x)*(a - grid(i, j).x) + grid(i, j).y * grid(i, j).y);
		//		temp2 = (a + grid(i, j).x) / sqrt((a + grid(i, j).x)*(a + grid(i, j).x) + grid(i, j).y * grid(i, j).y);
		//		if (temp1 >= a / r && temp2 >= a / r)
		//		{
		//			temp3 = min(grid(i, j).x * grid(i, j).x + (grid(i, j).y + sqrt(r*r - a*a))*(grid(i, j).y + sqrt(r*r - a*a)), grid(i, j).x * grid(i, j).x + (grid(i, j).y - sqrt(r*r - a*a))*(grid(i, j).y - sqrt(r*r - a*a)));
		//			levelSet(i, j) = sqrt(temp3) * ((grid(i, j) - 1).magnitude2() + 0.1);
		//			exactLevelSet(i, j) = sqrt(temp3);
		//		}
		//		else
		//		{
		//			temp3 = min(sqrt((grid(i, j).x + a)*(grid(i, j).x + a) + grid(i, j).y * grid(i, j).y) - r, sqrt((grid(i, j).x - a)*(grid(i, j).x - a) + grid(i, j).y * grid(i, j).y) - r);
		//			levelSet(i, j) = temp3 * ((grid(i, j) - 1).magnitude2() + 0.1);
		//			exactLevelSet(i, j) = temp3;
		//		}
		//	}
		//}
	}
	else if (example == 4)
	{
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.jEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (grid(i, j).magnitude() < 1)
				{
					levelSet(i, j) = -1.0;

				}
				else
				{
					levelSet(i, j) = 0.5;

				}
				exactLevelSet(i, j) = grid(i, j).magnitude() - 1.0;
			}
		}
	}
	else if (example == 5)
	{
		//for (int i = grid.iStart; i <= grid.jEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		if (abs(grid(i,j).x) < 0.5 || abs(grid(i, j).y)<0.5)
		//		{
		//			levelSet(i,j) = -1.0;
		//			exactLevelSet(i,j) = ;
		//		}
		//		else
		//		{
		//			levelSet(i, j) = 1.0;
		//		}
		//	}
		//}
	}
}

inline void LevelSetAdvection::surfReconstInitialCondition(const int & example)
{
	isVelocity = false;
	needReinitial = true;
	reinitialIter = 20;

	if (example == 1)
	{
		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		//distance.dataArray = 100;
		reconstructionVelocity = Field2D<double>(grid);
		givenPoint = VectorND<Vector2D<double>>(18);
		//dt = grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 500;
		reconstWriteIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		for (int i = 0; i < 5; i++)
		{
			givenPoint(i) = Vector2D<double>(0.3 + 0.1*double(i), 0.3) + grid.dx / 2;
			givenPoint(13 + i) = Vector2D<double>(0.3 + 0.1*double(i), 0.7) + grid.dx / 2;
		}
		for (int i = 0; i < 4; i++)
		{
			givenPoint(5 + i) = Vector2D<double>(0.3, 0.3 + 0.1*double(i)) + grid.dx / 2;
			givenPoint(9 + i) = Vector2D<double>(0.7, 0.3 + 0.1*double(i)) + grid.dx / 2;
		}
		sweepingDistance();

		// level set initialzation
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (abs(grid(i, j)(0) - 0.5)>abs(grid(i, j)(1) - 0.5))
				{
					levelSet(i, j) = abs(grid(i, j)(0) - 0.5) - 0.25;
				}
				else
				{
					levelSet(i, j) = abs(grid(i, j)(1) - 0.5) - 0.25;
				}
			}
		}
	}
	else if (example == 2)
	{

		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		//distance.dataArray = 100;
		reconstructionVelocity = Field2D<double>(grid);
		givenPointNum = 100;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum);
		//dt = 5 * grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 500;
		reconstWriteIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		for (int i = 0; i < givenPointNum; i++)
		{
			givenPoint(i) = 0.25*Vector2D<double>(cos(2 * PI*i / givenPointNum), sin(2 * PI*i / givenPointNum)) + 0.5 + grid.dx / 2;
		}

		//int i0, i1, j0, j1;
		//for (int k = 0; k < givenPointNum; k++)
		//{
		//	i0 = floor((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j0 = floor((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);
		//	i1 = ceil((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j1 = ceil((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);

		//	distance2Data(i0, j0);
		//	distance2Data(i0, j1);
		//	distance2Data(i1, j0);
		//	distance2Data(i1, j1);
		//}

		//distance2Data(); // distance initialization
		exactDistance();

		// level set initialization
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if ((grid(i, j) - 0.5 - grid.dx / 2).magnitude() < 0.25)
				{
					levelSet.phi(i, j) = -distance(i, j) - grid.dx * 3;
				}
				else
				{
					levelSet.phi(i, j) = distance(i, j) - grid.dx * 3;
				}

			}
		}

		//for (int i = grid.iStart; i <= grid.iEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		levelSet.phi(i, j) = -(abs((grid(i, j) - 0.5).x) + abs((grid(i, j) - 0.5).y) - 0.28*sqrt(2));
		//	}
		//}

	}
	else if (example == 3)
	{
		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		//distance.dataArray = 100;
		reconstructionVelocity = Field2D<double>(grid);
		givenPointNum = 100;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum);
		//dt = 5 * grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 500;
		reconstWriteIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		for (int i = 0; i < givenPointNum; i++)
		{
			givenPoint(i) = 0.25*Vector2D<double>(cos(2 * PI*i / givenPointNum), sin(2 * PI*i / givenPointNum)) + 0.5 + grid.dx / 2;
		}

		//int i0, i1, j0, j1;
		//for (int k = 0; k < givenPointNum; k++)
		//{
		//	i0 = floor((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j0 = floor((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);
		//	i1 = ceil((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j1 = ceil((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);

		//	distance2Data(i0, j0);
		//	distance2Data(i0, j1);
		//	distance2Data(i1, j0);
		//	distance2Data(i1, j1);
		//}

		//distance2Data(); // distance initialization
		exactDistance();

		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet.phi(i, j) = -(grid(i, j) - 0.5).magnitude() + 0.28;
			}
		}

		//for (int i = grid.iStart; i <= grid.iEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		levelSet.phi(i, j) = -(abs((grid(i, j) - 0.5).x) + abs((grid(i, j) - 0.5).y) - 0.28*sqrt(2));
		//	}
		//}
	}
	else if (example == 4)
	{
		cout << "Surface reconstruction : Two circles." << endl;

		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		//distance.dataArray = 100;
		reconstructionVelocity = Field2D<double>(grid);
		givenPointNum = 400;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum);
		//dt = 5 * grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 5000;
		reconstWriteIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		Vector2D<double> point1(0.6, 0.4);
		Vector2D<double> point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		//int i0, i1, j0, j1;
		//for (int k = 0; k < givenPointNum; k++)
		//{
		//	i0 = floor((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j0 = floor((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);
		//	i1 = ceil((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j1 = ceil((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);

		//	distance2Data(i0, j0);
		//	distance2Data(i0, j1);
		//	distance2Data(i1, j0);
		//	distance2Data(i1, j1);
		//}

		//distance2Data(); // distance initialization
		exactDistance();

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet.phi(i, j) = -(grid(i, j) - 0.5).magnitude() + 0.28;
			}
		}

		//for (int i = grid.iStart; i <= grid.iEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		levelSet.phi(i, j) = -(abs((grid(i, j) - 0.5).x) + abs((grid(i, j) - 0.5).y) - 0.28*sqrt(2));
		//	}
		//}
	}
	else if (example == 5)
	{
		cout << "Surface reconstruction : Two circles with outlier." << endl;
		cout << "Initial Level Set is far from real shape." << endl;

		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		reconstructionVelocity = Field2D<double>(grid);
		givenPointNum = 400;
		int outlier = 20;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum + outlier);
		reconstMaxIteration = 2000;
		reconstWriteIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		Vector2D<double> point1(0.6, 0.4);
		Vector2D<double> point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		srand(time(NULL));
		for (int i = 0; i < outlier; i++)
		{
			Vector2D<double> tempVector(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
			if (((tempVector-0.5)/3).magnitude()<0.28)
			{
				givenPoint(i + givenPointNum) = (tempVector - 0.5) / 3 + 0.5;
			}
			else
			{
				i--;
			}
			
		}
		givenPointNum = givenPointNum + outlier;

		exactDistance();

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet.phi(i, j) = -(grid(i, j) - 0.5).magnitude() + 0.28;
			}
		}
	}
	else if (example == 6)
	{
		cout << "Surface reconstruction : Two circles with outlier." << endl;
		cout << "Good Initial Level Set." << endl;

		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		reconstructionVelocity = Field2D<double>(grid);
		givenPointNum = 400;
		int outlier = 20;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum + outlier);
		reconstMaxIteration = 2000;
		reconstWriteIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		Vector2D<double> point1(0.6, 0.4);
		Vector2D<double> point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		srand(time(NULL));
		for (int i = 0; i < outlier; i++)
		{
			Vector2D<double> tempVector(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
			if (((tempVector - 0.5) / 3).magnitude()<0.28)
			{
				givenPoint(i + givenPointNum) = (tempVector - 0.5) / 3 + 0.5;
			}
			else
			{
				i--;
			}

		}
		givenPointNum = givenPointNum + outlier;

		exactDistance();

		double initialEpsilon = 0.02;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if ((grid(i, j) - point1 - grid.dx / 2).magnitude()<0.1 || (grid(i, j) - point2 - grid.dx / 2).magnitude()<0.1)
				{
					levelSet(i, j) = distance(i, j) + initialEpsilon;
				}
				else
				{
					levelSet(i, j) = -distance(i, j) + initialEpsilon;
				}
			}
		}
	}
}

inline void LevelSetAdvection::computeVelocity()
{
	double integralTerm = computeIntegralTerm();
	double lastTerm;
	Vector2D <double> gradPhi;

	double tempDist;
	double tempCurvature;

	levelSet.computeMeanCurvature();

	//ofstream solutionFile3;
	//solutionFile3.open("D:\\Data/curvature.txt", ios::binary);
	//for (int i = grid.iStart; i <= grid.iEnd; i++)
	//{
	//	for (int j = grid.jStart; j <= grid.jEnd; j++)
	//	{
	//		solutionFile3 << i << " " << j << " " << grid(i, j) << " " << levelSet.meanCurvature(i, j) <<" " << levelSet.dxxPhi(i, j)*levelSet.dyPhi(i, j)*levelSet.dyPhi(i, j) << " "<< -2.0*levelSet.dxyPhi(i, j)*levelSet.dxPhi(i, j)*levelSet.dyPhi(i, j) <<" " << levelSet.dyyPhi(i, j)*levelSet.dxPhi(i, j)*levelSet.dxPhi(i, j)<<" " << levelSet.dxPhi(i, j) << " " << levelSet.dyPhi(i, j) << " " << levelSet.dxyPhi(i, j)<<" " << levelSet.dxxPhi(i, j) << " " << levelSet.dyyPhi(i, j) << endl;
	//		cout <<i<<" "<<j<<" "<< levelSet.meanCurvature(i, j) << endl;
	//		////cout << -(levelSet.dxxPhi(i, j)*levelSet.dyPhi(i, j)*levelSet.dyPhi(i, j) - 2.0*levelSet.dxyPhi(i, j)*levelSet.dxPhi(i, j)*levelSet.dyPhi(i, j) + levelSet.dyyPhi(i, j)*levelSet.dxPhi(i, j)*levelSet.dxPhi(i, j)) << endl;
	//		cout << "dxx :" << levelSet.dxxPhi(i, j) << endl;
	//		cout << "dy  :" << levelSet.dyPhi(i, j) << endl;
	//		cout << "term:" << levelSet.dxxPhi(i, j)*levelSet.dyPhi(i, j)*levelSet.dyPhi(i, j) << endl;
	//		//cout << "dxy :" << levelSet.dxyPhi(i, j) << endl;
	//		////cout << "term:" << -2.0*levelSet.dxyPhi(i, j)*levelSet.dxPhi(i, j)*levelSet.dyPhi(i, j) << endl;
	//		cout << "dyy :" << levelSet.dyyPhi(i, j) << endl;
	//		cout << "dx  :" << levelSet.dxPhi(i, j) << endl;
	//		cout << "term:" << levelSet.dyyPhi(i, j)*levelSet.dxPhi(i, j)*levelSet.dxPhi(i, j) << endl;
	//		////cout << endl;	
	//		//cout << pow(levelSet.dxPhi(i, j)*levelSet.dxPhi(i, j) + levelSet.dyPhi(i, j)*levelSet.dyPhi(i, j) + DBL_EPSILON, 3.0 / 2.0) << endl;
	//	}
	//}
	//solutionFile3.close();

	//ofstream solutionFile2;
	//solutionFile2.open("D:\\Data/dotproduct.txt", ios::binary);

	//ofstream solutionFile4;
	//solutionFile4.open("D:\\Data/lastterm.txt", ios::binary);
	#pragma omp parallel for private(lastTerm, tempDist,tempCurvature, gradPhi)
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			reconstructionVelocity(i, j) = 0;
			gradPhi = levelSet.gradient(i, j);
			if (gradPhi.magnitude() < 3 * DBL_EPSILON)
			{
				lastTerm = 0;
			}
			else
			{
				//lastTerm = dotProduct(distance.gradient(i, j), gradPhi / (gradPhi.magnitude() + DBL_EPSILON)) + 1.0 / LpNorm*distance(i, j)*levelSet.meanCurvature(i, j);
				tempDist = min(distance(i, j), distanceThreshold);
				tempCurvature = AdvectionMethod2D<double>::sign(levelSet.meanCurvature(i, j))* min(abs(levelSet.meanCurvature(i, j)), curvatureThreshold);
				lastTerm = dotProduct(distance.gradient(i, j), gradPhi / (gradPhi.magnitude() + DBL_EPSILON)) + 1.0 / LpNorm*tempDist*tempCurvature;
			}
			//solutionFile2 << i << " " << j << " " << grid(i, j) << " " << dotProduct(distance.gradient(i, j), gradPhi / (gradPhi.magnitude() + DBL_EPSILON)) << endl;
			//solutionFile4 << i << " " << j << " " << grid(i, j) << " " << lastTerm << endl;

			//reconstructionVelocity(i, j) = gradPhi.magnitude()*integralTerm*pow(distance(i, j), LpNorm - 1)*lastTerm;
			//reconstructionVelocity(i, j) = gradPhi.magnitude()*integralTerm*pow(min(distance(i, j), distanceThreshold), LpNorm - 1)*lastTerm;
			reconstructionVelocity(i, j) = gradPhi.magnitude()*integralTerm*lastTerm;

			//if (i < 1 && j < 1)
			//{
			//	cout << i << " " << j << endl;
			//	cout << "inte : " << integralTerm << endl;
			//	cout << "dist : " << distance(i, j) << endl;
			//	cout << "disG : " << distance.gradient(i, j) << endl;
			//	cout << "phiG : " << gradPhi << endl;
			//	cout << "phiG : " << gradPhi.magnitude() << endl;
			//	cout << "dotP : " << dotProduct(distance.gradient(i, j), gradPhi / (gradPhi.magnitude() + DBL_EPSILON)) << endl;
			//	cout << "curv : " << levelSet.meanCurvature(i, j) << endl;
			//	cout << "last : " << lastTerm << endl;
			//	cout << "velo : " << velocity(i, j) << endl;
			//	cout << endl;
			//}
		}
	}
	//solutionFile2.close();
	//solutionFile4.close();

}

inline double LevelSetAdvection::computeIntegralTerm()
{
	double sum = 0.0;
	
#pragma omp parallel for reduction(+:sum)
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			sum += pow(distance(i, j), LpNorm)*AdvectionMethod2D<double>::deltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude()*grid.dx*grid.dy;
			//cout << i << " " << j << endl;
			//cout << "distance : " << pow(distance(i, j), LpNorm) << endl;
			//cout << "deltaft  : " << AdvectionMethod2D<double>::deltaFt(levelSet(i, j)) << endl;
			//cout << "gradient : " << levelSet.gradient(i, j) << endl;
			//cout << "magnitude: " << levelSet.gradient(i, j).magnitude() << endl;
			//cout << sum << endl;
			//cout << endl;
			//cout << i << " " << j << " " << sum << endl;
		}
	}
	return pow(sum, 1.0 / LpNorm - 1.0);
}

inline void LevelSetAdvection::levelSetPropagatingTVDRK3()
{
	LevelSet2D originLevelSet = levelSet;
	LevelSet2D tempLevelSet(originLevelSet.grid);

	Field2D<double> k1(levelSet.grid);
	Field2D<double> k2(levelSet.grid);
	Field2D<double> k3(levelSet.grid);

	Field2D<double> wenoXMinus(levelSet.grid);
	Field2D<double> wenoXPlus(levelSet.grid);
	Field2D<double> wenoYMinus(levelSet.grid);
	Field2D<double> wenoYPlus(levelSet.grid);

	AdvectionMethod2D<double>::WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = dt*reconstructionVelocity(i, j);
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}

	AdvectionMethod2D<double>::WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = dt*reconstructionVelocity(i, j);
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}

	AdvectionMethod2D<double>::WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = dt*reconstructionVelocity(i, j);
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*levelSet(i, j) + 2.0 / 3.0*k3(i, j);
		}
	}
}

inline bool LevelSetAdvection::stoppingCriterion()
{
	bool criterion = true;

	for (int i = 0; i < givenPointNum; i++)
	{
		if (abs(levelSet.interpolation(givenPoint(i)))>grid.dx + grid.dy)
		{
			criterion = false;
			return criterion;
		}
	}

	return criterion;
}


inline double LevelSetAdvection::distance2Data(const int & i, const int & j)
{
	double distance = 100;
	double tempDist;

	for (int k = 0; k < givenPointNum; k++)
	{
		tempDist = (givenPoint(k) - grid(i, j)).magnitude();
		if (distance > tempDist)
		{
			distance = tempDist;
		}
	}

	return distance;
}

inline void LevelSetAdvection::exactDistance()
{
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			distance(i, j) = distance2Data(i, j);
		}
	}
}

inline void LevelSetAdvection::sweepingDistance()
{
	double xMin, yMin;
	double h = min(grid.dx, grid.dy);

	for (int i = distance.iStart; i <= distance.iEnd; i++)
	{
		for (int j = distance.jStart; j <= distance.jEnd; j++)
		{
			xMin = min(distance(max(i - 1, distance.iStart), j), distance(min(i + 1, distance.iEnd), j));
			yMin = min(distance(i, max(j - 1, distance.jStart)), distance(i, min(j + 1, distance.jEnd)));
			//#pragma omp parallel for
			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				distance(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}


	for (int i = distance.iStart; i <= distance.iEnd; i++)
	{
		for (int j = distance.jEnd; j >= distance.jStart; j--)
		{
			xMin = min(distance(max(i - 1, distance.iStart), j), distance(min(i + 1, distance.iEnd), j));
			yMin = min(distance(i, max(j - 1, distance.jStart)), distance(i, min(j + 1, distance.jEnd)));

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				distance(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}



	for (int i = distance.iEnd; i >= distance.iStart; i--)
	{
		for (int j = distance.jStart; j <= distance.jEnd; j++)
		{
			xMin = min(distance(max(i - 1, distance.iStart), j), distance(min(i + 1, distance.iEnd), j));
			yMin = min(distance(i, max(j - 1, distance.jStart)), distance(i, min(j + 1, distance.jEnd)));

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				distance(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}


	for (int i = distance.iEnd; i >= distance.iStart; i--)
	{
		for (int j = distance.jEnd; j >= distance.jStart; j--)
		{
			xMin = min(distance(max(i - 1, distance.iStart), j), distance(min(i + 1, distance.iEnd), j));
			yMin = min(distance(i, max(j - 1, distance.jStart)), distance(i, min(j + 1, distance.jEnd)));

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				distance(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}
}

inline double LevelSetAdvection::adaptiveTimeStep()
{
	return cflCondition*max(grid.dx, grid.dy);
}

inline double LevelSetAdvection::adaptiveTimeStep(const Field2D<double>& velocity1)
{
	double maxVel1 = 0;

	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			if (abs(velocity1(i, j)) > maxVel1)
			{
				maxVel1 = abs(velocity1(i, j));
			}
		}
	}
	return cflCondition*max(grid.dx, grid.dy) / maxVel1;
}

inline double LevelSetAdvection::adaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2)
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

inline void LevelSetAdvection::outputResult(const int & iter)
{
	cout << "Write results" << endl;

	ofstream solutionFile1;
	solutionFile1.open("D:\\Data/phi" + to_string(iter) + ".txt", ios::binary);
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			solutionFile1 << i << " " << j << " " << grid(i, j) << " " << levelSet(i, j) << endl;
		}
	}
	solutionFile1.close();


	if (isSurfaceReconstruction)
	{
		ofstream solutionFile2;
		solutionFile2.open("D:\\Data/velocity" + to_string(iter) + ".txt", ios::binary);
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				solutionFile2 << i << " " << j << " " << grid(i, j) << " " << reconstructionVelocity(i, j) << endl;
			}
		}
		solutionFile2.close();
	}


	if (iter == 0)
	{
		if (isPropagation)
		{
			ofstream solutionFile2;
			solutionFile2.open("D:\\Data/velocity.txt", ios::binary);
			for (int i = grid.iStart; i <= grid.iEnd; i++)
			{
				for (int j = grid.jStart; j <= grid.jEnd; j++)
				{
					solutionFile2 << i << " " << j << " " << grid(i, j) << " " << velocityX(i, j) << " " << velocityY(i, j) << endl;
				}
			}
			solutionFile2.close();
		}

		if (isSurfaceReconstruction)
		{
			ofstream solutionFile3;
			solutionFile3.open("D:\\Data/distance.txt", ios::binary);
			for (int i = grid.iStart; i <= grid.iEnd; i++)
			{
				for (int j = grid.jStart; j <= grid.jEnd; j++)
				{
					solutionFile3 << i << " " << j << " " << grid(i, j) << " " << distance(i, j) << endl;
				}
			}
			solutionFile3.close();


			ofstream solutionFile4;
			solutionFile4.open("D:\\Data/pointData.txt", ios::binary);
			for (int i = 0; i <= givenPoint.iEnd; i++)
			{
				solutionFile4 << givenPoint(i) << endl;
			}
			solutionFile4.close();
		}


	}
}
