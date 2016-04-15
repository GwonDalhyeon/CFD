#pragma once
#include "CommonDef.h"
#include "Grid2D.h"
#include "Field2D.h"
#include "CombineStructure.h"
#include "LevelSet2D.h"
#include "VoronoiDiagram2D.h"
#include "DelaunayTriangle.h"

template<class TT>
class PointIntegralMethod
{
public:
	Grid2D grid;
	
	int pointNum;
	VectorND<Vector2D<double>> pointCloud;
	
	int innerPointNum;
	VectorND<int> innerPointIndex;
	
	int bdryPointNum;
	VectorND<int> bdryPointIndex;

	int nbhdNum;
	Array2D<int> nNearstNbhdIndex;
	Array2D<double> nNearstNbhdDist;

	double t;
	double beta;



	PointIntegralMethod();
	~PointIntegralMethod();

	inline void InitialCondition(const int & example);

	inline void PointIntegralMethodnSolver(int example);


private:

};

template<class TT>
PointIntegralMethod<TT>::PointIntegralMethod()
{
}

template<class TT>
PointIntegralMethod<TT>::~PointIntegralMethod()
{
}

template<class TT>
inline void PointIntegralMethod<TT>::InitialCondition(const int & example)
{
	if (example == 1)
	{
		cout << "******************************************************" << endl;
		cout << "       Point Integral Method : One circles." << endl;
		cout << "******************************************************" << endl;

		grid = Grid2D(-1, 1, 2001, -1, 1, 2001);

		innerPointNum = 1200;
		innerPointIndex = VectorND<int>(1, innerPointNum);

		bdryPointNum = 100;
		bdryPointIndex = VectorND<int>(1, bdryPointNum);

		pointNum = innerPointNum + bdryPointNum;
		pointCloud = VectorND<Vector2D<double>>(1, pointNum);

		nbhdNum = 20;
		nNearstNbhdIndex = Array2D<int>(1, pointNum, 1, nbhdNum);
		nNearstNbhdDist = Array2D<double>(1, pointNum, 1, nbhdNum);

#pragma omp parallel for
		for (int i = 1; i <= bdryPointNum; i++)
		{
			pointCloud(i) = 0.5*Vector2D<double>(cos(2 * PI*i / bdryPointNum), sin(2 * PI*i / bdryPointNum));
			bdryPointIndex(i) = i;
		}

		Vector2D<double> tempVector;
		srand(time(NULL));
		for (int i = 1; i <= innerPointNum; i++)
		{
			tempVector= Vector2D<double>(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
			if (((tempVector - 0.5)).magnitude()<0.5-grid.dx/2)
			{
				pointCloud(i + bdryPointNum) = (tempVector - 0.5);
				innerPointIndex(i) = i + bdryPointNum;
			}
			else
			{
				i--;
			}
		}

		t;
		beta;
	}
	else if (example == 2)
	{

	}
}

template<class TT>
inline void PointIntegralMethod<TT>::PointIntegralMethodnSolver(int example)
{
	InitialCondition(example);

	grid.Variable();
	VecND2DVariable("pointData", pointCloud);
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	MATLAB.Command("plot(pointData(:,1), pointData(:,2),'ro');axis([-1 1 -1 1]);grid on");
	
	Field2D<double> section(grid);

	Voronoi<double> tempVoronoi;// (grid, pointCloud);
	tempVoronoi.Make(grid, pointCloud, section);

}
