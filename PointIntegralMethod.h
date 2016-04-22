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
		bdryPointNum = 300;

		pointNum = innerPointNum + bdryPointNum;
		innerPointIndex = VectorND<int>(1, pointNum);
		bdryPointIndex = VectorND<int>(1, pointNum);
		pointCloud = VectorND<Vector2D<double>>(1, pointNum);

		nbhdNum = 20;
		nNearstNbhdIndex = Array2D<int>(1, pointNum, 1, nbhdNum);
		nNearstNbhdDist = Array2D<double>(1, pointNum, 1, nbhdNum);

		Vector2D<double> tempVector;
		for (int i = 1; i <= bdryPointNum; i++)
		{
			pointCloud(i) = 0.5*Vector2D<double>(cos(2 * PI*i / bdryPointNum), sin(2 * PI*i / bdryPointNum));
			bdryPointIndex(i) = 1;

			for (int j = i; j >= 2; j--)
			{
				if (pointCloud(j).x < pointCloud(j - 1).x)
				{
					tempVector = pointCloud(j);
					pointCloud(j) = pointCloud(j - 1);
					pointCloud(j - 1) = tempVector;
				}
				else
				{
					break;
				}
			}
		}

		double bdryBand = 2 * PI / double(bdryPointNum) / 2;
		int temp;
		srand(time(NULL));
		for (int i = 1; i <= innerPointNum; i++)
		{
			tempVector = Vector2D<double>(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
			if (((tempVector - 0.5)).magnitude()<0.5-grid.dx/2- bdryBand)
			{
				pointCloud(i + bdryPointNum) = tempVector - 0.5;
				innerPointIndex(i + bdryPointNum) = 1;

				for (int j = i + bdryPointNum; j >= 2; j--)
				{
					if (pointCloud(j).x < pointCloud(j - 1).x)
					{
						tempVector = pointCloud(j);
						pointCloud(j) = pointCloud(j - 1);
						pointCloud(j - 1) = tempVector;

						temp = bdryPointIndex(j);
						bdryPointIndex(j) = bdryPointIndex(j - 1);
						bdryPointIndex(j - 1) = temp;
						
						temp = innerPointIndex(j);
						innerPointIndex(j) = innerPointIndex(j - 1);
						innerPointIndex(j - 1) = temp;
					}
					else
					{
						break;
					}
				}
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

	innerPointIndex.Variable("inner");
	bdryPointIndex.Variable("bdry");

	grid.Variable();
	VecND2DVariable("pointData", pointCloud);
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	MATLAB.Command("plot(pointData(:,1), pointData(:,2),'ro');axis([-1 1 -1 1]);grid on;axis equal");
	
	VectorND<Polygon2D> voronoi(pointCloud.iStart, pointCloud.iLength);

	VoronoiDiagram<double>::Make(grid, pointCloud, innerPointIndex, bdryPointIndex, voronoi);

	VectorND<double> area(pointCloud.iStart, pointCloud.iLength);

	double totalArea = 0;
#pragma omp parallel for reduction(+:totalArea)
	for (int i = pointCloud.iStart; i <= pointCloud.iLength; i++)
	{
		area(i) = voronoi(i).Area();
		totalArea += area(i);
	}

	MATLAB.Variable("area", totalArea);


}
