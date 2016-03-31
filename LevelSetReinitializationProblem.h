#pragma once
#include "CommonDef.h"
#include "Grid2D.h"
#include "Field2D.h"
#include "CombineStructure.h"
#include "LevelSet2D.h"
#include "AdvectionMethod2D.h"

class Reinitialzation
{
public:
	Grid2D grid;
	LevelSet2D exactLevelSet;
	LevelSet2D levelSet;

	double cflCondition;
	double dt;
	int maxIteration;
	int writeIter;

	Reinitialzation();
	~Reinitialzation();

	inline void InitialCondition(const int & example);
	inline void ReinitializationSolver(const int & example);

	double AdaptiveTimeStep();

	inline void OutputResult(const int & iter);
private:

};

Reinitialzation::Reinitialzation()
{
}

Reinitialzation::~Reinitialzation()
{
}

inline void Reinitialzation::InitialCondition(const int & example)
{
	grid = Grid2D(-2, 2, 101, -2, 2, 101);

	//dt = grid.dx*grid.dy;
	maxIteration = 120;
	writeIter = 10;

	exactLevelSet = LevelSet2D(grid);
	levelSet = LevelSet2D(grid);

	double a = 0.7;
	double r = 1.0;

	if (example == 1) 
	{
		cout << "*******************************************************" << endl;
		cout << "    A circle with center at the origen and radius 1" << endl;
		cout << "*******************************************************" << endl;


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
	else if (example == 2) 
	{
		cout << "*******************************************************" << endl;
		cout << "    Two circles of radius r are placed at (+-a,0)" << endl;
		cout << "         and a sqruar on the plane. Let 0<a<r,"<< endl;
		cout<<"        so that the two circles intersect each other." << endl;
		cout << "*******************************************************" << endl;

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
	else if (example == 3) 
	{

		cout << "*******************************************************" << endl;
		cout << "    Two circles of radius r are placed at(+-a, 0)" << endl;
		cout << "         on the plane.Let 0<a<r," << endl;
		cout << "        so that the two circles intersect each other." << endl;
		cout << "*******************************************************" << endl;

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

inline void Reinitialzation::ReinitializationSolver(const int & example)
{
	bool writeFile = false;
	string str;
	const char*cmd;

	cflCondition = 0.8;
	InitialCondition(example);

	if (writeFile)
	{
		str = "phi0";
		levelSet.phi.WriteFile(str);
	}

	grid.Variable();
	levelSet.phi.Variable("phi0");

	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	levelSet.phi.Variable("phi0");
	MATLAB.Command("subplot(1, 2, 1)");
	MATLAB.Command("surf(X,Y,phi0)");
	MATLAB.Command("subplot(1, 2, 2)");
	MATLAB.Command("contour(X, Y, phi0, [0 0],'b');");
	MATLAB.Command("grid on");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);

	//OutputResult(0);

	/////////////////////////
	////                /////
	////    Iteration   /////
	////                /////
	/////////////////////////
	for (int i = 1; i <= maxIteration; i++)
	{

		cout << "Reinitialization : " << i << endl;
		dt = AdaptiveTimeStep();
		AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
		levelSet.phi.Variable("phi");

		MATLAB.Command("subplot(1, 2, 1)");
		MATLAB.Command("surf(X,Y,phi)");
		MATLAB.Command("subplot(1, 2, 2)");
		MATLAB.Command("plot(0,0);");
		MATLAB.Command("contour(X, Y, phi, [0 0],'b');");
		MATLAB.Command("hold on");
		MATLAB.Command("contour(X, Y, phi);");
		MATLAB.Command("contour(X, Y, phi0,[0 0],'r');");
		MATLAB.Command("grid on");
		MATLAB.Command("hold off");
		str = string("title(['iteration : ', num2str(") + to_string(i) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);

		if (writeFile && i%writeIter == 0)
		{
			str = "phi" + to_string(i);
			levelSet.phi.WriteFile(str);
		}
	}
}

inline double Reinitialzation::AdaptiveTimeStep()
{
	return cflCondition*max(grid.dx, grid.dy);
}