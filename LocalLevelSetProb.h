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
		LLS.InitialTube();

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
		LLS.InitialTube();
		cflCondition = 0.05;

		//dt = grid.dx*grid.dy;
		maxIteration = 1000;
		writeIter = 10;
	}
	else if (example == 3)
	{

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
	LLS.phi.Variable("phi");
	LLS.Tube.Variable("Tube");
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	MATLAB.Command("subplot(1,2,1)");
	MATLAB.Command("surf(X, Y, Tube);grid on;axis([-1 1 -1 1]);axis equal;");
	MATLAB.Command("subplot(1,2,2)");
	MATLAB.Command("surf(X,Y,phi);");

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

	//MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	//MATLAB.Command("v = VideoWriter('newfile.avi');");
	//MATLAB.Command("open(v)");
	clock_t before;
	double  result;
	
	double totalT = 0;
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << endl;
		cout << "********************************" << endl;
		cout << "Level set advection : " << i << endl;

		if (isVelocity)
		{
			dt = AdaptiveTimeStep(velocityX, velocityY);
			totalT += dt;
			before = clock();
			AdvectionMethod2D<double>::LLSPropagatingTVDRK3(LLS, velocityX, velocityY, dt);
			result = (double)(clock() - before) / CLOCKS_PER_SEC;
			cout << result << endl;
			LLS.phi.Variable("phi");
			MATLAB.Command("subplot(1,2,1)");
			MATLAB.Command("surf(X,Y,phi);");
			str = string("title(['iteration : ', num2str(") + to_string(i) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
			MATLAB.Command("subplot(1,2,2)");
			MATLAB.Command("contour(X, Y, phi0, [0 0],'b');hold on;grid on;contour(X, Y, phi,[0 0],'r');axis([-1 1 -1 1]);axis equal;hold off;");
		}
		else
		{
			dt = AdaptiveTimeStep();
			before = clock();
			AdvectionMethod2D<double>::LLSPropagatingTVDRK3(LLS, dt);
			result = (double)(clock() - before) / CLOCKS_PER_SEC;
			cout << result << endl;
			LLS.phi.Variable("phi");
			MATLAB.Command("subplot(1,2,1)");
			MATLAB.Command("surf(X,Y,phi);");
			str = string("title(['iteration : ', num2str(") + to_string(i) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
			MATLAB.Command("subplot(1,2,2)");
			MATLAB.Command("contour(X, Y, phi0, [0 0],'b');hold on;grid on;contour(X, Y, phi,[0 0],'r');axis([-1 1 -1 1]);axis equal;hold off;");
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