#pragma once

#include "CommonDef.h"
#include "Vector2D.h"
#include "VectorND.h"
#include "Field2D.h"
#include "Grid2D.h"

class Polygon2D
{
public:
	int nGon;
	VectorND<Vector2D<double>> Points;
	VectorND<int> Index;

	Polygon2D();
	~Polygon2D();

	Polygon2D(const int & ipNGon);

	inline Vector2D<double> & operator [](const int& i) const;

	inline Vector2D<double> & operator ()(const int& i) const;

	inline void operator = (const Polygon2D & ipPoly);

	inline Vector2D<Vector2D<double>> EdgePoint(const int & edge);
	inline Vector2D<int> EdgeIndex(const int & edge);

	inline void Plot();
	inline double Area();
private:

};

Polygon2D::Polygon2D()
{
}

Polygon2D::~Polygon2D()
{
}

inline Polygon2D::Polygon2D(const int & ipNGon)
{
	nGon = ipNGon;
	Points = VectorND<Vector2D<double>>(1, nGon);
	Index = VectorND<int>(1, nGon);
}

inline Vector2D<double> & Polygon2D::operator[](const int & i) const
{
	return Points(i);
}

inline Vector2D<double> & Polygon2D::operator()(const int & i) const
{
	return Points(i);
}

inline void Polygon2D::operator=(const Polygon2D & ipPoly)
{
	nGon = ipPoly.nGon;
	Points = ipPoly.Points;
	Index = ipPoly.Index;
}

inline Vector2D<Vector2D<double>> Polygon2D::EdgePoint(const int & edge)
{
	assert(edge >= 1 || edge <= nGon);

	Vector2D<Vector2D<double>> returnEdge;
	returnEdge(0) = Points(edge);
	returnEdge(1) = Points(edge % nGon + 1);

	return returnEdge;
}

inline Vector2D<int> Polygon2D::EdgeIndex(const int & edge)
{
	assert(edge >= 1 || edge <= nGon);
	return Vector2D<int>(Index(edge), Index(edge% nGon + 1));
}

inline void Polygon2D::Plot()
{
	VecND2DVariable("polygon", Points);
	MATLAB.Command("plot([polygon(:, 1); polygon(1, 1)], [polygon(:, 2); polygon(1, 2)])");
}

inline double Polygon2D::Area()
{
	assert(3 <= nGon);

	double area = 0;
#pragma omp parallel for reduction(+:area)
	for (int i = 1; i <= nGon; i++)
	{
		area += (Points(i%nGon + 1).x - Points(i).x)*(Points(i%nGon + 1).y + Points(i).y) / 2;
	}

	return abs(area);
}
