#pragma once
#include "DelaunayTriangle.h"

template<class TT>
class VoronoiDiagram
{
public:
	VoronoiDiagram();
	~VoronoiDiagram();


	static void Make(const Grid2D & grid, const VectorND<Vector2D<double>> & ipPoints, const VectorND<int> & innerIndex, const VectorND<int> & bdryIndex, VectorND<Polygon2D> & rVoronoi);

	static void FindNbhdPolygon(const VectorND<Polygon2D> & polygon, const VectorND<int> & innerIndex, const VectorND<int> & bdryIndex, VectorND<int> & nbhdNum, VectorND<VectorND<int>> & nbhdPoly);

	static void RearrangeNbhdPolygonClockwise(const VectorND<Vector2D<double>>& ipPoints, const VectorND<int> & innerIndex, const VectorND<int> & bdryIndex, VectorND<VectorND<int>>& nbhdPoly, VectorND<Vector2D<double>> & center);

	static void BdryAddedVoronoi(const VectorND<Vector2D<double>>& ipPoints, const VectorND<int> & bdryIndex, const VectorND<Polygon2D> & ipPolygon, const VectorND<VectorND<int>>& nbhdPoly, const VectorND<Vector2D<double>>& polyCenter, VectorND<Polygon2D> & rVoronoi);

	static double Angle(const Vector2D<double> & P1, const Vector2D<double> & P2);
private:

};

template<class TT>
inline VoronoiDiagram<TT>::VoronoiDiagram()
{
}

template<class TT>
inline VoronoiDiagram<TT>::~VoronoiDiagram()
{
}

template<class TT>
inline void VoronoiDiagram<TT>::Make(const Grid2D & grid, const VectorND<Vector2D<double>> & ipPoints, const VectorND<int> & innerIndex, const VectorND<int> & bdryIndex, VectorND<Polygon2D> & rVoronoi)
{
	VectorND<Polygon2D> Triangles;
	VectorND<Vector2D<double>> TriCenter;
	VectorND<double> radius;
	DelaunayTriangulization::DelaunayTriangulate(ipPoints, Triangles, TriCenter, radius);

	//for (int i = Triangles.iStart; i <= Triangles.iEnd; i++)
	//{
	//	TriCenter(i) = 0;
	//	TriCenter(i) += Triangles(i)(1) / 3;
	//	TriCenter(i) += Triangles(i)(2) / 3;
	//	TriCenter(i) += Triangles(i)(3) / 3;
	//}

	int maxNbhdTri = 20;
	VectorND<VectorND<int>> nbhdTri(ipPoints.iStart, ipPoints.iLength);
#pragma omp parallel for
	for (int i = ipPoints.iStart; i <= ipPoints.iEnd; i++)
	{
		nbhdTri(i) = VectorND<int>(1, maxNbhdTri);
	}

	VectorND<int> nbhdNum(ipPoints.iStart, ipPoints.iLength);
	FindNbhdPolygon(Triangles, innerIndex, bdryIndex, nbhdNum, nbhdTri);

	RearrangeNbhdPolygonClockwise(ipPoints, innerIndex, bdryIndex, nbhdTri, TriCenter);

	BdryAddedVoronoi(ipPoints, bdryIndex, Triangles, nbhdTri, TriCenter, rVoronoi);

}

template<class TT>
inline void VoronoiDiagram<TT>::FindNbhdPolygon(const VectorND<Polygon2D>& polygon, const VectorND<int> & innerIndex, const VectorND<int> & bdryIndex, VectorND<int>& nbhdNum, VectorND<VectorND<int>>& nbhdPoly)
{
	int vIdx;
	//#pragma omp parallel for private(vIdx) shared(nbhdNum)
	for (int i = polygon.iStart; i <= polygon.iEnd; i++)
	{
		//cout << polygon(i).Index << endl;
		for (int j = 1; j <= polygon(i).nGon; j++)
		{
			vIdx = polygon(i).Index(j);
			nbhdNum(vIdx)++;
			nbhdPoly(vIdx)(nbhdNum(vIdx)) = i;
		}
	}

	VectorND<int> tempV;
#pragma omp parallel for private(tempV)
	for (int i = nbhdPoly.iStart; i <= nbhdPoly.iEnd; i++)
	{
		tempV = nbhdPoly(i);
		nbhdPoly(i) = VectorND<int>(tempV.iStart, nbhdNum(i));
		for (int j = 1; j <= nbhdNum(i); j++)
		{
			nbhdPoly(i)(j) = tempV(j);
		}
	}
}

template<class TT>
inline void VoronoiDiagram<TT>::RearrangeNbhdPolygonClockwise(const VectorND<Vector2D<double>>& ipPoints, const VectorND<int> & innerIndex, const VectorND<int> & bdryIndex, VectorND<VectorND<int>>& nbhdPoly, VectorND<Vector2D<double>>& center)
{
	int maxNbhdTri = 20;
	Vector2D<double> currentPoint;
	Vector2D<double> P1, P2;
	VectorND<VectorND<int>> angleRank(ipPoints.iStart, ipPoints.iLength);
	VectorND<VectorND<double>> angle(ipPoints.iStart, ipPoints.iLength);
	VectorND<int> tempV;


#pragma omp parallel for private(currentPoint, P1, P2, tempV)
	for (int i = ipPoints.iStart; i <= ipPoints.iEnd; i++)
	{


		angle(i) = VectorND<double>(nbhdPoly(i).iStart, nbhdPoly(i).iLength);
		angleRank(i) = VectorND<int>(nbhdPoly(i).iStart, nbhdPoly(i).iLength);


		currentPoint = ipPoints(i);
		P1 = center(nbhdPoly(i)(1));

		for (int j = 1; j <= nbhdPoly(i).iLength; j++)
		{
			P2 = center(nbhdPoly(i)(j));
			angle(i)(j) = Angle(P2 - currentPoint, P1 - currentPoint);
		}

		for (int j = 1; j <= nbhdPoly(i).iLength; j++)
		{
			for (int k = 1; k <= nbhdPoly(i).iLength; k++)
			{
				if (angle(i)(k) <= angle(i)(j))
				{
					angleRank(i)(j)++;
				}
			}
		}

		tempV = nbhdPoly(i);
		for (int j = nbhdPoly(i).iStart; j <= nbhdPoly(i).iEnd; j++)
		{
			nbhdPoly(i)(angleRank(i)(j)) = tempV(j);
		}
	}

}

template<class TT>
inline void VoronoiDiagram<TT>::BdryAddedVoronoi(const VectorND<Vector2D<double>>& ipPoints, const VectorND<int>& bdryIndex, const VectorND<Polygon2D>& ipPolygon, const VectorND<VectorND<int>>& nbhdPoly, const VectorND<Vector2D<double>>& polyCenter, VectorND<Polygon2D> & rVoronoi)
{

	int maxNbhdTri = 20;
	Vector2D<double> currentPoint;
	Vector2D<double> P1, P2;
	int tempI;
	VectorND<VectorND<int>> angleRank(ipPoints.iStart, ipPoints.iLength);
	VectorND<VectorND<double>> angle(ipPoints.iStart, ipPoints.iLength);
	VectorND<int> tempV;
	int addedPointNum;
	VectorND<Vector2D<double>> addedPoint(1, 3);

	VectorND<Vector2D<double>> tempCenter(polyCenter.iStart, 2 * polyCenter.iLength);
	VectorND<VectorND<int>> tempNbhdPoly(ipPoints.iStart, ipPoints.iLength);

#pragma omp parallel for
	for (int i = polyCenter.iStart; i <= polyCenter.iEnd; i++)
	{
		tempCenter(i) = polyCenter(i);
	}

	int temp2;
	int tempPolyNum = ipPolygon.iEnd;
	//#pragma omp parallel for private(currentPoint, P1, P2, tempV)
	for (int i = ipPoints.iStart; i <= ipPoints.iEnd; i++)
	{
		if (bdryIndex(i))
		{
			addedPointNum = 2;
			temp2 = 0;
			angle(i) = VectorND<double>(nbhdPoly(i).iStart, nbhdPoly(i).iLength + addedPointNum);
			angleRank(i) = VectorND<int>(nbhdPoly(i).iStart, nbhdPoly(i).iLength + addedPointNum);

			tempNbhdPoly(i) = VectorND<int>(1, nbhdPoly(i).iLength + addedPointNum);

			//tempPolyNum++;
			//temp2++;
			//tempNbhdPoly(i)(nbhdPoly(i).iLength + temp2) = tempPolyNum;
			//tempCenter(tempPolyNum) = ipPoints(i);

			for (int j = nbhdPoly(i).iStart; j <= nbhdPoly(i).iLength; j++)
			{
				tempNbhdPoly(i)(j) = nbhdPoly(i)(j);
				for (int k = 1; k <= ipPolygon(nbhdPoly(i)(j)).nGon; k++)
				{
					tempI = ipPolygon(nbhdPoly(i)(j)).Index(k);
					if (bdryIndex(tempI) && tempI != i)
					{
						tempPolyNum++;
						temp2++;
						tempNbhdPoly(i)(nbhdPoly(i).iLength + temp2) = tempPolyNum;
						tempCenter(tempPolyNum) = ipPoints(i) / 2 + ipPoints(tempI) / 2;
					}
				}
			}

			currentPoint = ipPoints(i);
			P1 = tempCenter(nbhdPoly(i)(1));
			for (int j = nbhdPoly(i).iStart; j <= nbhdPoly(i).iEnd; j++)
			{
				P2 = tempCenter(tempNbhdPoly(i)(j));
				angle(i)(j) = Angle(P2 - currentPoint, P1 - currentPoint);
			}

			for (int j = nbhdPoly(i).iEnd + 1; j <= nbhdPoly(i).iEnd + 2; j++)
			{
				P2 = tempCenter(tempNbhdPoly(i)(j));
				angle(i)(j) = Angle(P2 - currentPoint, P1 - currentPoint);
			}

			for (int j = 1; j <= tempNbhdPoly(i).iLength; j++)
			{
				for (int k = 1; k <= tempNbhdPoly(i).iLength; k++)
				{
					if (angle(i)(k) <= angle(i)(j))
					{
						angleRank(i)(j)++;
					}
				}
			}
			tempV = tempNbhdPoly(i);
#pragma omp parallel for
			for (int j = tempNbhdPoly(i).iStart; j <= tempNbhdPoly(i).iEnd; j++)
			{
				tempNbhdPoly(i)(angleRank(i)(j)) = tempV(j);
			}

		}
		else
		{
			tempCenter(i) = polyCenter(i);
			tempNbhdPoly(i) = nbhdPoly(i);
		}
	}

#pragma omp parallel for
	for (int i = rVoronoi.iStart; i <= rVoronoi.iEnd; i++)
	{
		rVoronoi(i) = Polygon2D(tempNbhdPoly(i).iLength);
		for (int j = tempNbhdPoly(i).iStart; j <= tempNbhdPoly(i).iEnd; j++)
		{
			rVoronoi(i)(j) = tempCenter(tempNbhdPoly(i)(j));
		}
	}


	bool isShow = true;
	if (isShow)
	{
		string str;
		const char* cmd;

		MATLAB.Command("figure('units','normalized','outerposition',[0 0 1/2 1])");
		MATLAB.Command("plot(pointData(:,1), pointData(:,2),'ro');axis([-1/2 1/2 -1/2 1/2]);grid on;axis equal");
		MATLAB.Command("hold on");
		for (int i = rVoronoi.iStart; i <= rVoronoi.iEnd; i++)
		{
			//MATLAB.Variable("i", i);
			//MATLAB.Command("plot(pointData(i,1), pointData(i,2),'ro');;grid on");
			MATLAB.Command("axis equal");

			//rVoronoi(i).Plot("b");
			rVoronoi(i).Plot();

			str = string("title(['Voronoi : ', num2str(") + to_string(i) + string("),'/', num2str(") + to_string(rVoronoi.iLength) + string(")]);");
			cmd = str.c_str();
			MATLAB.Command(cmd);
		}
	}
}

template<class TT>
inline double VoronoiDiagram<TT>::Angle(const Vector2D<double>& P1, const  Vector2D<double>& P2)
{
	double alpha = atan2(P1.y, P1.x);
	double beta = atan2(P2.y, P2.x);

	return beta - alpha;
}



