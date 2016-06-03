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

	double Distance2Data(const int& i, const int& j);
	void ExactDistance();

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
		cflCondition = 0.5;

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

		cflCondition = 0.5;

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

		cflCondition = 0.5;

		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 10;
	}
	else if (example ==4)
	{
		cout << "Level set advection Test - Unit square" << endl;

		isVelocity = false;
		needReinitial = false;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);


		levelSet = LevelSet2D(grid);
//#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = max(abs(grid(i,j).x), abs(grid(i, j).y)) - 0.5;
			}
		}

		cflCondition = 0.05;

		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 10;
	}
	else if (example == 5)
	{
		cout << "Level set advection Test - Unit sircle" << endl;

		isVelocity = false;
		needReinitial = false;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);


		levelSet = LevelSet2D(grid);
		//#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i,j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - 0.5;
			}
		}

		cflCondition = 0.05;

		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 10;
	}
	else if (example == 6)
	{
		cout << "Level set advection Test - Diamond" << endl;

		isVelocity = false;
		needReinitial = false;

		grid = Grid2D(-1.5, 1.5, 301, -1.5, 1.5, 301);


		levelSet = LevelSet2D(grid);
		//#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = (4 * abs(grid(i, j).x) + abs(grid(i, j).y)) - 1;
			}
		}

		cflCondition = 0.05;

		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 10;
	}
	else if (example == 7)
	{
		cout << "Level set advection Test - Concave" << endl;

		isVelocity = false;
		needReinitial = false;

		grid = Grid2D(-1.5, 1.5, 301, -1.5, 1.5, 301);

		int edgePointNum = 200;
		givenPointNum = edgePointNum * 5;
		givenPoint = VectorND<Vector2D<double>>(1, givenPointNum);

		for (int i = givenPoint.iStart; i <= edgePointNum; i++)
		{
			givenPoint(i).x = 0.25 / edgePointNum * i;
			givenPoint(i).y = (1 - 4 * abs(givenPoint(i).x));
			givenPoint(i + edgePointNum).x = 0.25 / edgePointNum * i;
			givenPoint(i + edgePointNum).y = -(1 - 4 * abs(givenPoint(i + edgePointNum).x));
			givenPoint(i + edgePointNum * 2).x = 0.25 / edgePointNum * i - 0.25;
			givenPoint(i + edgePointNum * 2).y = (1 - 4 * abs(givenPoint(i + edgePointNum * 2).x));
		}
		for (int i = 1; i <= edgePointNum; i++)
		{
			givenPoint(i + edgePointNum * 3).x = 0.25 / edgePointNum * i - 0.25;
			givenPoint(i + edgePointNum * 3).y = 0;
		}
		for (int i = 1; i <= edgePointNum; i++)
		{
			givenPoint(i + edgePointNum * 4).x = 0;
			givenPoint(i + edgePointNum * 4).y = -1 / double(edgePointNum) * i;

		}

		distance = Field2D<double>(grid);
		ExactDistance();
		levelSet = LevelSet2D(grid);
		//#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if ((4 * abs(grid(i, j).x) + abs(grid(i, j).y)) - 1 >= 0)
				{
					levelSet(i, j) = distance(i, j);
				}
				else if (grid(i,j).x<=0 && grid(i,j).y<=0)
				{
					levelSet(i, j) = distance(i, j);
				}
				else
				{
					levelSet(i, j) = -distance(i, j);
				}
			}
		}

		cflCondition = 0.05;

		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 10;
	}
}


inline void LevelSetAdvection::AdvectionSolver(const int & example)
{
	bool writeFile = true;
	string str;
	const char*cmd;

	InitialCondition(example);

	grid.Variable();
	levelSet.phi.Variable("phi0");

	
	if (writeFile)
	{
		str = "phi0";
		levelSet.phi.WriteFile(str);
		
		if (isVelocity)
		{
			str = "velocityX";
			velocityX.WriteFile(str);
			str = "velocityY";
			velocityY.WriteFile(str);
		}
	}
	
	
	//MATLAB.Command("figure('units','normalized','outerposition',[0 0 1/2 1])");
	//MATLAB.Command("v = VideoWriter('newfile.avi');");
	//MATLAB.Command("open(v)");
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << "Level set advection : " << i << endl;

		if (isVelocity)
		{
			dt = AdaptiveTimeStep(velocityX, velocityY);
			AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, velocityX, velocityY, dt);
			levelSet.phi.Variable("phi");

			MATLAB.Command("contour(X, Y, phi0, [0 0],'b');");
			MATLAB.Command("hold on");
			MATLAB.Command("grid on");
			MATLAB.Command("axis([-1 1 -1 1]);axis equal;");
			MATLAB.Command("contour(X, Y, phi, [0 0],'r');");
			MATLAB.Command("hold off");
		}
		else
		{
			dt = AdaptiveTimeStep();
			AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, dt);
			//MATLAB.Variable("i", i);
			//levelSet.phi.Variable("phi");
			//MATLAB.Command("axis([-1 1 -1 1]);axis equal;");
			//MATLAB.Command("contour(X, Y, phi0, [0 0],'b');");
			//MATLAB.Command("hold on");
			//MATLAB.Command("grid on");
			//MATLAB.Command("axis([-1 1 -1 1]);axis equal;");
			////MATLAB.Command("hold on");
			//MATLAB.Command("contour(X, Y, phi, [0 0],'r');");
			//MATLAB.Command("hold off");
			////MATLAB.Command("F=getframe;");
			////MATLAB.Command("writeVideo(v,F)");
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
	//MATLAB.Command("close(v)");
	
}

inline double LevelSetAdvection::Distance2Data(const int & i, const int & j)
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

inline void LevelSetAdvection::ExactDistance()
{
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			distance(i, j) = Distance2Data(i, j);
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