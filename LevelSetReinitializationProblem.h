#pragma once
#include <iostream>

#include "AdvectionMethod2D.h"

class Reinitialzation
{
public:
	Grid2D grid;
	LevelSet2D exactLevelSet;
	LevelSet2D levelSet;
	double dt;
	int maxIteration;
	int writeIter;

	Reinitialzation();
	~Reinitialzation();

	void initialCondition(const int& example);

	void reinitializationSolver(const int& example);
	void outputResult(const int& iter);
	void outputResult(const int& iter, const string& method);
private:

};

Reinitialzation::Reinitialzation()
{
}

Reinitialzation::~Reinitialzation()
{
}

inline void Reinitialzation::initialCondition(const int& example)
{
	grid = Grid2D(-2, 2, 101, -2, 2, 101);

	dt = grid.dx*grid.dy;
	maxIteration = 1000;
	writeIter = 10;

	exactLevelSet = LevelSet2D(grid);
	levelSet = LevelSet2D(grid);

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

inline void Reinitialzation::reinitializationSolver(const int& example)
{
	initialCondition(example);


	outputResult(0, "TVDRK");
	for (int i = 1; i <= maxIteration; i++)
	{

		cout << "Reinitialization TVD R-K : " << i << endl;
		AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, 5 * dt);
		if (i%writeIter == 0)
		{
			outputResult(i, "TVDRK");
		}
	}

	//outputResult(0, "FE");
	//for (int i = 1; i <= maxIteration; i++)
	//{
	//	cout << "Reinitialization Forward Euler : " << i << endl;
	//	AdvectionMethod2D<double>::levelSetReinitializationFE(levelSet, 5 * dt);
	//	if (i%writeIter == 0)
	//	{
	//		outputResult(i);
	//		outputResult(i,"FE");
	//	}
	//}

	//outputResult(0, "GS");
	//for (int i = 1; i <= maxIteration; i++)
	//{
	//	cout << "Reinitialization Gauss Seidel : " << i << endl;
	//	AdvectionMethod2D<double>::levelSetReinitializationGS(levelSet, 5 * dt);
	//	if (i%writeIter == 0)
	//	{
	//		outputResult(i, "GS");
	//	}
	//}


}

inline void Reinitialzation::outputResult(const int & iter)
{
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

	if (iter == 0)
	{
		ofstream solutionFile2;
		solutionFile2.open("D:\\Data/phi.txt", ios::binary);
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				solutionFile2 << i << " " << j << " " << grid(i, j) << " " << exactLevelSet(i, j) << endl;
			}
		}
		solutionFile2.close();

	}

}

inline void Reinitialzation::outputResult(const int & iter, const string & method)
{
	ofstream solutionFile1;
	solutionFile1.open("D:\\Data/phi" + method + to_string(iter) + ".txt", ios::binary);
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			solutionFile1 << i << " " << j << " " << grid(i, j) << " " << levelSet(i, j) << endl;
		}
	}
	solutionFile1.close();

	if (iter == 0)
	{
		ofstream solutionFile2;
		solutionFile2.open("D:\\Data/phi.txt", ios::binary);
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				solutionFile2 << i << " " << j << " " << grid(i, j) << " " << exactLevelSet(i, j) << endl;
			}
		}
		solutionFile2.close();

	}
}
