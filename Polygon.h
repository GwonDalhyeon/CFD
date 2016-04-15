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
	VectorND<Vector2D<double>> Vertices;
	VectorND<int> VerticesIndex;

	Polygon2D();
	~Polygon2D();

	Polygon2D(const int & ipNGon);

	inline Vector2D<double> & operator [](const int& i) const;

	inline Vector2D<double> & operator ()(const int& i) const;

	inline void operator = (const Polygon2D & ipPoly);

	inline Vector2D<Vector2D<double>> EdgePoint(const int & edge);
	inline Vector2D<int> EdgeIndex(const int & edge);

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
	Vertices = VectorND<Vector2D<double>>(1, nGon);
	VerticesIndex = VectorND<int>(1, nGon);
}

inline Vector2D<double> & Polygon2D::operator[](const int & i) const
{
	return Vertices(i);
}

inline Vector2D<double> & Polygon2D::operator()(const int & i) const
{
	return Vertices(i);
}

inline void Polygon2D::operator=(const Polygon2D & ipPoly)
{
	nGon = ipPoly.nGon;
	Vertices = ipPoly.Vertices;
	VerticesIndex = ipPoly.VerticesIndex;
}

inline Vector2D<Vector2D<double>> Polygon2D::EdgePoint(const int & edge)
{
	assert(edge >= 1 || edge <= nGon);

	Vector2D<Vector2D<double>> returnEdge;
	returnEdge(0) = Vertices(edge);
	returnEdge(1) = Vertices((edge + 1) % nGon);

	return returnEdge;
}

inline Vector2D<int> Polygon2D::EdgeIndex(const int & edge)
{
	assert(edge >= 1 || edge <= nGon);
	return Vector2D<int>(VerticesIndex(edge), VerticesIndex((edge + 1) % nGon));
}

inline double Polygon2D::Area()
{
	assert(3 <= nGon);

	double area = 0;
#pragma omp parallel for reduction(+:area)
	for (int i = 1; i <= nGon; i++)
	{
		area += (Vertices(i%nGon + 1).x - Vertices(i).x)*(Vertices(i%nGon + 1).y + Vertices(i).y) / 2;
	}

	return abs(area);
}
