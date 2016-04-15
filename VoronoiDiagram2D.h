#pragma once
#include "Bitmap.h"

template<class TT>
class Voronoi
{
public:
	Voronoi();
	~Voronoi();


	inline void Make(const Grid2D & grid, const VectorND<Vector2D<TT>> & ipPoints, Field2D<double> & ipSection);

	Field2D<double> R;
	Field2D<double> G;
	Field2D<double> B;

	VectorND<double> colorMapR;
	VectorND<double> colorMapG;
	VectorND<double> colorMapB;
	VectorND<Vector2D<TT>> points;

	inline void CreateColors();

	inline void CreateSites(const Grid2D & grid, Field2D<double> & ipSection);

private:

};

template<class TT>
inline Voronoi<TT>::Voronoi()
{
}

template<class TT>
inline Voronoi<TT>::~Voronoi()
{
}

template<class TT>
inline void Voronoi<TT>::Make(const Grid2D & grid, const VectorND<Vector2D<TT>> & ipPoints, Field2D<double> & ipSection)
{
	R = Field2D<double>(grid);
	G = Field2D<double>(grid);
	B = Field2D<double>(grid);

	colorMapR = VectorND<double>(ipPoints.iStart, ipPoints.iLength);
	colorMapG = VectorND<double>(ipPoints.iStart, ipPoints.iLength);
	colorMapB = VectorND<double>(ipPoints.iStart, ipPoints.iLength);
	points = ipPoints;

	CreateColors();
	CreateSites(grid, ipSection);
	//SetSitesPoints();
}


template<class TT>
inline void Voronoi<TT>::CreateColors()
{
	int tempRand;

	for (int i = points.iStart; i <= points.iEnd; i++)
	{
		tempRand = rand();
		colorMapR(i) = (double)(tempRand) / RAND_MAX;
		tempRand = rand();
		colorMapG(i) = (double)(tempRand) / RAND_MAX;
		tempRand = rand();
		colorMapB(i) = (double)(tempRand) / RAND_MAX;
	}
}

template<class TT>
inline void Voronoi<TT>::CreateSites(const Grid2D & grid, Field2D<double> & ipSection)
{
	int ind = -1;
	TT dist = 1000000;
	TT d;
	Vector2D<TT> p;

#pragma omp parallel for private(ind, dist, p, d)
	for (int i = R.iStart; i <= R.iEnd; i++)
	{
		for (int j = R.jStart; j <= R.jEnd; j++)
		{
			ind = -1;
			dist = INT_MAX;
			for (int it = points.iStart; it <= points.iEnd; it++)
			{
				p = points[it];
				d = (p - grid(i, j)).magnitude();
				//d = DistanceSqrd(p, ww, hh);
				if (d < dist)
				{
					dist = d;
					ind = it;
				}
			}

			if (ind > -1)
			{
				R(i, j) = colorMapR(ind);
				G(i, j) = colorMapG(ind);
				B(i, j) = colorMapB(ind);
				ipSection(i, j) = ind;
			}
		}
	}
	R.Variable("R");
	G.Variable("G");
	B.Variable("B");
	MATLAB.Command("Section(:,:,1)=R;Section(:,:,2)=G;Section(:,:,3)=B;");
	MATLAB.Command("figure, imshow(Section);");
}

