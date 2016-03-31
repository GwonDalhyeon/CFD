#pragma once
#include "CommonDef.h"
#include "Grid2D.h"
#include "Field2D.h"
#include "CombineStructure.h"
#include "LevelSet2D.h"
#include "AdvectionMethod2D.h"

class LevelSetAdvection
{
public:
	Grid2D grid;
	LevelSet2D levelSet;

	double cflCondition;
	double dt;

	// Level set propagation velocity.
	bool isVelocity;
	bool needReinitial;
	Field2D<double> velocityX;
	Field2D<double> velocityY;

	int reinitialIter;
	int maxIteration;
	int writeIter;


	// Zero level set point.
	int givenPointNum;
	VectorND<Vector2D<double>> givenPoint;
	Field2D<double> distance;


	LevelSetAdvection();
	~LevelSetAdvection();


	void InitialCondition(const int & example);
	void AdvectionSolver(const int & example);

	// Adaptive time step functions.
	double AdaptiveTimeStep();
	double AdaptiveTimeStep(const Field2D<double>& velocity1);
	double AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2);
private:

};

//#endif // !LevelSetAdvectionProblem_H




LevelSetAdvection::LevelSetAdvection()
{
}

LevelSetAdvection::~LevelSetAdvection()
{
}

inline void LevelSetAdvection::InitialCondition(const int & example)
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
		maxIteration = 1000;
		writeIter = 100;
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
		maxIteration = 1000;
		writeIter = 10;
	}
	else if (example == 3)
	{
		cout << "Level set advection Test - Seven-point Star" << endl;
		
		isVelocity = false;
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

				//levelSet(i, j) = distance2Data(i, j);
			}
		}


		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 10;
	}
}


inline void LevelSetAdvection::AdvectionSolver(const int & example)
{
	bool writeFile = false;
	string str;
	const char*cmd;

	cflCondition;
	InitialCondition(example);

	grid.Variable();
	levelSet.phi.Variable("phi0");

	
	if (writeFile)
	{
		str = "phi0";
		levelSet.phi.WriteFile(str);
		str = "velocityX";
		velocityX.WriteFile(str);
		str = "velocityY";
		velocityY.WriteFile(str);
	}
	
	

	for (int i = 1; i <= maxIteration; i++)
	{
		cout << "Level set advection : " << i << endl;

		if (isVelocity)
		{
			dt = AdaptiveTimeStep(velocityX, velocityY);
			AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, velocityX, velocityY, dt);
		}
		else
		{
			dt = AdaptiveTimeStep();
			AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, dt);
		}


		if (needReinitial)
		{
			for (int j = 0; j < reinitialIter; j++)
			{
				cout << "Reinitialization : " << i << "-" << j + 1 << endl;
				dt = AdaptiveTimeStep();
				AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
			}
		}

		if (i%writeIter == 0 && writeFile)
		{
			str = "phi" + to_string(i);
			levelSet.phi.WriteFile(str);
		}
	}
	
}


inline double LevelSetAdvection::AdaptiveTimeStep()
{
	return cflCondition*max(grid.dx, grid.dy);
}

inline double LevelSetAdvection::AdaptiveTimeStep(const Field2D<double>& velocity1)
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

inline double LevelSetAdvection::AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2)
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