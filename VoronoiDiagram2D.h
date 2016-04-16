#pragma once
#include "DelaunayTriangle.h"

template<class TT>
class VoronoiDiagram
{
public:
	VoronoiDiagram();
	~VoronoiDiagram();


	static void Make(const Grid2D & grid, const VectorND<Vector2D<TT>> & ipPoints, VectorND<Polygon2D> & rVoronoi);

	static void FindNbhdPolygon(const VectorND<Polygon2D> & polygon, VectorND<int> nbhdNum, VectorND<VectorND<int>> & nbhdPoly);
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
inline void VoronoiDiagram<TT>::Make(const Grid2D & grid, const VectorND<Vector2D<TT>> & ipPoints, VectorND<Polygon2D> & rVoronoi)
{
	VectorND<Polygon2D> Triangles;
	VectorND<Vector2D<double>> center;
	VectorND<double> radius;
	DelaunayTriangulization::DelaunayTriangulate(ipPoints, Triangles, center, radius);


	int maxNbhdTri = 20;
	VectorND<VectorND<int>> nbhdTri(ipPoints.iStart, ipPoints.iLength);
#pragma omp parallel for
	for (int i = ipPoints.iStart; i <= ipPoints.iEnd; i++)
	{
		nbhdTri(i) = VectorND<int>(1, maxNbhdTri);
	}

	VectorND<int> nbhdNum(ipPoints.iStart, ipPoints.iLength);
	FindNbhdPolygon(Triangles, nbhdNum, nbhdTri);


}

template<class TT>
inline void VoronoiDiagram<TT>::FindNbhdPolygon(const VectorND<Polygon2D>& polygon, VectorND<int> nbhdNum, VectorND<VectorND<int>>& nbhdPoly)
{
	int vIdx;
#pragma omp parallel for private(vIdx)
	for (int i = polygon.iStart; i <= polygon.iEnd; i++)
	{
		for (int j = 1; j <= polygon(i).nGon; j++)
		{
			vIdx = polygon(i).Index(j);
			nbhdNum(vIdx)++;
			nbhdPoly(vIdx, nGon(vIdx)) = i;
		}
	}
}



