#pragma once

#include "AdvectionMethod2D.h"

class LocalLevelSetAdvection
{
public:
	Grid2D grid;
	LS LLS;
	double cflCondition;
	double dt;

	bool isVelocity;
	bool needReinitial;
	FD velocityX;
	FD velocityY;

	int reinitialIter;
	int maxIteration;
	int writeIter;

	// Zero level set point.
	int givenPointNum;
	VectorND<VT> givenPoint;
	FD distance;

	LocalLevelSetAdvection();
	~LocalLevelSetAdvection();

	void InitialCondition(const int & example);
	void AdvectionSolver(const int & example);

	double Distance2Data(const int& i, const int& j);
	void ExactDistance();

	// Adaptive time step functions.
	double AdaptiveTimeStep();
	double AdaptiveTimeStep(const FD& velocity1);
	double AdaptiveTimeStep(const FD& velocity1, const FD& velocity2);

private:

};

LocalLevelSetAdvection::LocalLevelSetAdvection()
{
}

LocalLevelSetAdvection::~LocalLevelSetAdvection()
{
}

inline void LocalLevelSetAdvection::InitialCondition(const int & example)
{
	if (example == 1)
	{
		isVelocity = true;
		needReinitial = false;

		cout << "Local Level set advection Test - Rotating a circle" << endl;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);

		LLS = LS(grid, 3*grid.dx);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				LLS(i, j) = sqrt((grid(i, j).x - 0.5)*(grid(i, j).x - 0.5) + grid(i, j).y*grid(i, j).y) - 0.25;
			}
		}
		LLS.UpdateInterface();

		velocityX = FD(grid);
		velocityY = FD(grid);
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


		LLS = LS(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				LLS(i, j) = sqrt((grid(i, j).x - 0.5)*(grid(i, j).x - 0.5) + grid(i, j).y*grid(i, j).y) - 0.25;
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
		givenPoint = VectorND<VT>(givenPointNum);

		double s;
#pragma omp parallel for private (s)
		for (int i = 0; i < givenPointNum; i++)
		{
			s = double(i) / double(givenPointNum);
			givenPoint(i) = (0.1 + 0.065*sin(7 * 2 * PI*s))*VT(cos(2 * PI*s), sin(2 * PI*s));
		}

		LLS = LS(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{

				//LLS(i, j) = distance2Data(i, j);
			}
		}

		cflCondition = 0.5;

		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 10;
	}
	else if (example == 4)
	{
		cout << "Level set advection Test - Unit square" << endl;

		isVelocity = false;
		needReinitial = false;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);


		LLS = LS(grid);
		//#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				LLS(i, j) = max(abs(grid(i, j).x), abs(grid(i, j).y)) - 0.5;
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


		LLS = LS(grid);
		//#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				LLS(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - 0.5;
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


		LLS = LS(grid);
		//#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				LLS(i, j) = (4 * abs(grid(i, j).x) + abs(grid(i, j).y)) - 1;
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
		givenPoint = VectorND<VT>(1, givenPointNum);

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

		distance = FD(grid);
		ExactDistance();
		LLS = LS(grid);
		//#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if ((4 * abs(grid(i, j).x) + abs(grid(i, j).y)) - 1 >= 0)
				{
					LLS(i, j) = distance(i, j);
				}
				else if (grid(i, j).x <= 0 && grid(i, j).y <= 0)
				{
					LLS(i, j) = distance(i, j);
				}
				else
				{
					LLS(i, j) = -distance(i, j);
				}
			}
		}

		cflCondition = 0.05;

		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 10;
	}
}


inline void LocalLevelSetAdvection::AdvectionSolver(const int & example)
{
	bool writeFile = false;
	string str;
	const char*cmd;

	InitialCondition(example);

	grid.Variable();
	LLS.phi.Variable("phi0");


	if (writeFile)
	{
		str = "phi0";
		LLS.phi.WriteFile(str);

		if (isVelocity)
		{
			str = "velocityX";
			velocityX.WriteFile(str);
			str = "velocityY";
			velocityY.WriteFile(str);
		}
	}

	//LLS.T1.Variable("T1");
	//LLS.T2.Variable("T2");
	//MATLAB.Command("subplot(1,2,1)");
	//MATLAB.Command("plot3(X,Y,T1,'o')");
	//MATLAB.Command("subplot(1,2,2)");
	//MATLAB.Command("plot3(X,Y,T2,'o')");

	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1/2 1])");
	//MATLAB.Command("v = VideoWriter('newfile.avi');");
	//MATLAB.Command("open(v)");
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << "Level set advection : " << i << endl;

		if (isVelocity)
		{
			dt = AdaptiveTimeStep(velocityX, velocityY);
			AdvectionMethod2D<double>::LLSPropagatingTVDRK3(LLS, velocityX, velocityY, dt);
			LLS.phi.Variable("phi");
			MATLAB.Command("subplot(1,2,1)");
			MATLAB.Command("contour(X, Y, phi0, [0 0],'b');hold on;grid on;axis([-1 1 -1 1]);axis equal;contour(X, Y, phi, [0 0],'r');");
			MATLAB.Command("hold off");
			MATLAB.Command("subplot(1,2,2)");
			MATLAB.Command("surf(X,Y,phi);");
		}
		else
		{
			dt = AdaptiveTimeStep();
			AdvectionMethod2D<double>::LLSPropagatingTVDRK3(LLS, dt);
			MATLAB.Variable("i", i);
			LLS.phi.Variable("phi");
			MATLAB.Command("contour(X, Y, phi0, [0 0],'b');hold on;grid on;contour(X, Y, phi, [0 0],'r');axis([-1 1 -1 1]);axis equal;hold off;");
			//MATLAB.Command("F=getframe;");
			//MATLAB.Command("writeVideo(v,F)");
		}


		if (needReinitial)
		{
			for (int j = 0; j < reinitialIter; j++)
			{
				cout << "Reinitialization : " << i << "-" << j + 1 << endl;
				dt = AdaptiveTimeStep();
				AdvectionMethod2D<double>::LLSReinitializationTVDRK3(LLS, dt);
			}
		}

		if (i%writeIter == 0 && writeFile)
		{
			str = "phi" + to_string(i);
			LLS.phi.WriteFile(str);
		}
	}
	//MATLAB.Command("close(v)");

}

inline double LocalLevelSetAdvection::Distance2Data(const int & i, const int & j)
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

inline void LocalLevelSetAdvection::ExactDistance()
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


inline double LocalLevelSetAdvection::AdaptiveTimeStep()
{
	return cflCondition*max(grid.dx, grid.dy);
}

inline double LocalLevelSetAdvection::AdaptiveTimeStep(const FD& velocity1)
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

inline double LocalLevelSetAdvection::AdaptiveTimeStep(const FD& velocity1, const FD& velocity2)
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