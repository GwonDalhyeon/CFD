#pragma once

#include "Array1D.h"
#include "Array2D.h"
#include "Field1D.h"

class EulerEquation1DSolver
{
public:
	EulerEquation1DSolver();
	~EulerEquation1DSolver();

	Grid1D grid;

	Field1D<double> density;
	Field1D<double> pressure;
	Field1D<double> energy;

	double cflCondition;
	double dt;


	inline void InitialCondition(const int & example);

private:

};

EulerEquation1DSolver::EulerEquation1DSolver()
{
}

EulerEquation1DSolver::~EulerEquation1DSolver()
{
}

inline void EulerEquation1DSolver::InitialCondition(const int & example)
{
	if (example==1)
	{
		grid = Grid1D(-1, 1, 201);
		density = Field1D<double>(grid);
		pressure = Field1D<double>(grid);
		energy = Field1D<double>(grid);

		cflCondition = 0.5;

	}
	if (example==2)
	{

	}
}
