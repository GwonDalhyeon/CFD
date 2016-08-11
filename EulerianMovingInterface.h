#pragma once

#include "AdvectionMethod2D.h"

class MovingInterface
{
public:
	MovingInterface();
	~MovingInterface();

	Grid2D grid;

	Field2D<double> U; // x velocity
	Field2D<double> V; // y velocity
	Field2D<double> Surfactant;
	Field2D<double> SurfactantOld;

	LevelSet2D levelSet;

	Array2D<double>A;
	// CG solver 1
	CSR<double> Acsr;
	// CG solver 2
	VectorND<double> a;
	VectorND<int> row;
	VectorND<int> col;
	int nonzeroNum;

	double dt;
	double totalT;

	double cflCondition;

	int maxIteration;
	int writeOutputIteration;

	inline void InitialCondition(const int& example);
	inline void MovingInterfaceSolver(const int& example);

	inline void OneStepSemiImplicit(const int&example);
	inline void TwoStepSemiImplicit(const int&example);

	inline void SurfactantNormalTerm(const Field2D<double>& ipField, Array2D<double>& term);
	inline double ExactSurfactant(const int& example, const double& x, const double& y, const double& time);
	inline void GenerateLinearSystem1(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystem1(VectorND<double>& vectorB, const double & scaling);
	inline void GenerateLinearSystem2(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystem2(VectorND<double>& vectorB, const double & scaling);

private:

};

MovingInterface::MovingInterface()
{
}

MovingInterface::~MovingInterface()
{
}

inline void MovingInterface::InitialCondition(const int & example)
{
	if (example == 1)
	{
		cout << "*************************" << endl;
		cout << "       An Eulerian Formulation  " << endl;
		cout << " for Solving PDE along Moving Interface" << endl;
		cout << "          --JJ Xu, HK Zhao-- " << endl;
		cout << "               Example 1 " << endl;
		cout << "*************************" << endl;
		int gridSize = 40;
		grid = Grid2D(-2, 2, gridSize + 1, -2, 2, gridSize + 1);

		U = Field2D<double>(grid);
		V = Field2D<double>(grid);

		levelSet = LevelSet2D(grid);
		double radius = 1.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}
		Surfactant = Field2D<double>(grid);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				Surfactant(i, j) = grid(i, j).y / (sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) + DBL_EPSILON) + 2;
			}
		}
		SurfactantOld = Surfactant;

		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		A = Array2D<double>(1, (grid.iRes - 2)*(grid.jRes - 2), 1, (grid.iRes - 2)*(grid.jRes - 2));
		dt = grid.dx / 4.0;
		maxIteration = 100;
	}

	if (example == 2)
	{

	}
}

inline void MovingInterface::MovingInterfaceSolver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	const char* cmd;

	totalT = 0;
	InitialCondition(example);

	grid.Variable();

	Surfactant.Variable("Surfactant");
	MATLAB.Command("figure,surf(X,Y,Surfactant);");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);

	levelSet.ComputeMeanCurvature();
	levelSet.ComputeNormal();

	GenerateLinearSystem1(A, 1);
	A.Variable("A1");
	//// CG solver 1
	Acsr = CSR<double>(A);
	//// CG solver 2
	CGSolver::SparseA(A, a, row, col, nonzeroNum);


	//// Step 1
	cout << "********************************" << endl;
	cout << "       Iteration " << to_string(1) << " : Start" << endl;

	totalT += dt;
	OneStepSemiImplicit(example);

	cout << "       Iteration " << to_string(1) << " : End" << endl;
	cout << "********************************" << endl;

	Surfactant.Variable("SurfactantNew");
	MATLAB.Command("figure,surf(X,Y,SurfactantNew);");
	str = string("title(['iteration : ', num2str(") + to_string(1) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);

	GenerateLinearSystem2(A, 1);
	A.Variable("A2");
	//// CG solver 1
	Acsr = CSR<double>(A);
	//// CG solver 2
	CGSolver::SparseA(A, a, row, col, nonzeroNum);

	for (int i = 2; i <= 5; i++)
	{
		cout << "********************************" << endl;
		cout << "       Iteration " << to_string(i) << " : Start" << endl;

		totalT += dt;
		TwoStepSemiImplicit(example);

		cout << "       Iteration " << to_string(i) << " : End" << endl;
		cout << "********************************" << endl;

		Surfactant.Variable("SurfactantNew");
		MATLAB.Command("figure,surf(X,Y,SurfactantNew);");
		str = string("title(['iteration : ', num2str(") + to_string(i) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
	}

}

inline void MovingInterface::OneStepSemiImplicit(const int&example)
{
	//// Copy Old Data.
	SurfactantOld = Surfactant;

	//// Linear Equation
	VectorND<double> vectorB((grid.iRes - 2)*(grid.jRes - 2));
	GenerateLinearSystem1(vectorB, 1);
	vectorB.Variable("vectorB");

	VectorND<double> tempSur((grid.iRes - 2)*(grid.jRes - 2));

	int solver = 1;
	if (solver == 1)
	{
		tempSur = CGSolver::SolverCSR(Acsr, vectorB, grid.dx*grid.dy);
	}
	else if (solver == 2)
	{
		CGSolver::SolverSparse(A.iRes, a, row, col, vectorB, tempSur);
	}
	tempSur.Variable("sur");

	int index;
#pragma omp parallel for private(index)
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			index = (i - (grid.iStart + 1)) + (j - (grid.jStart + 1))*(grid.iRes - 2);
			Surfactant(i, j) = tempSur(index);
		}
	}

	//// Dirichlet Boundary Condition from the Exact Solution
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		Surfactant(i, grid.jStart) = ExactSurfactant(example, grid(i, grid.jStart).x, grid(i, grid.jStart).y, totalT);
		Surfactant(i, grid.jEnd) = ExactSurfactant(example, grid(i, grid.jEnd).x, grid(i, grid.jEnd).y, totalT);
	}
#pragma omp parallel for
	for (int j = grid.jStart; j <= grid.jEnd; j++)
	{
		Surfactant(grid.iStart, j) = ExactSurfactant(example, grid(grid.iStart, j).x, grid(grid.iStart, j).y, totalT);
		Surfactant(grid.iEnd, j) = ExactSurfactant(example, grid(grid.iEnd, j).x, grid(grid.iEnd, j).y, totalT);
	}


	//	//// Exact Solution
	//#pragma omp parallel for
	//	for (int i = grid.iStart; i <= grid.iEnd; i++)
	//	{
	//		for (int j = grid.jStart; j <= grid.jEnd; j++)
	//		{
	//			Surfactant(i, j) = ExactSurfactant(example, grid(i, j).x, grid(i, j).y, totalT);
	//		}
	//	}
	//	Surfactant.Variable("SurfactantExact");
	//	MATLAB.Command("figure,surf(X,Y,SurfactantExact);");
}

inline void MovingInterface::TwoStepSemiImplicit(const int&example)
{
	//// Linear Equation
	VectorND<double> vectorB((grid.iRes - 2)*(grid.jRes - 2));
	GenerateLinearSystem2(vectorB, 1);
	vectorB.Variable("vectorB");

	VectorND<double> tempSur((grid.iRes - 2)*(grid.jRes - 2));

	int solver = 1;
	if (solver == 1)
	{
		tempSur = CGSolver::SolverCSR(Acsr, vectorB, grid.dx*grid.dy);
	}
	else if (solver == 2)
	{
		CGSolver::SolverSparse(A.iRes, a, row, col, vectorB, tempSur);
	}
	tempSur.Variable("sur");

	int index;
#pragma omp parallel for private(index)
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			index = (i - (grid.iStart + 1)) + (j - (grid.jStart + 1))*(grid.iRes - 2);
			Surfactant(i, j) = tempSur(index);
		}
	}

	//// Dirichlet Boundary Condition from the Exact Solution
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		Surfactant(i, grid.jStart) = ExactSurfactant(example, grid(i, grid.jStart).x, grid(i, grid.jStart).y, totalT);
		Surfactant(i, grid.jEnd) = ExactSurfactant(example, grid(i, grid.jEnd).x, grid(i, grid.jEnd).y, totalT);
	}
#pragma omp parallel for
	for (int j = grid.jStart; j <= grid.jEnd; j++)
	{
		Surfactant(grid.iStart, j) = ExactSurfactant(example, grid(grid.iStart, j).x, grid(grid.iStart, j).y, totalT);
		Surfactant(grid.iEnd, j) = ExactSurfactant(example, grid(grid.iEnd, j).x, grid(grid.iEnd, j).y, totalT);
	}

}

inline void MovingInterface::SurfactantNormalTerm(const Field2D<double>& ipField, Array2D<double>& term)
{
	Field2D<Vector2D<double>> gradient = Field2D<double>::Gradient(ipField);
	Vector2D<double> normal;
	Array2D<double> Hessian(2, 2);
#pragma omp parallel for private(normal, Hessian)
	for (int i = term.iStart; i <= term.iEnd; i++)
	{
		for (int j = term.jStart; j <= term.jEnd; j++)
		{
			normal = levelSet.normal(i, j);
			Hessian = Surfactant.Hessian(i, j);
			term(i, j) = 0;
			term(i, j) += -levelSet.meanCurvature(i, j)*dotProduct(normal, gradient(i, j));
			term(i, j) += -normal(0)*(Hessian(0, 0)*normal(0) + Hessian(0, 1)*normal(1));
			term(i, j) += -normal(1)*(Hessian(1, 0)*normal(0) + Hessian(1, 1)*normal(1));
		}
	}
}

inline double MovingInterface::ExactSurfactant(const int & example, const double & x, const double & y, const double & time)
{
	if (example == 1)
	{
		return exp(-time / (x*x + y*y + DBL_EPSILON))*(y / sqrt(x*x + y*y + DBL_EPSILON)) + 2;
	}
	return 0.0;
}

inline void MovingInterface::GenerateLinearSystem1(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int innerIStart = grid.iStart + 1;
	int innerIEnd = grid.iEnd - 1;
	int innerJStart = grid.jStart + 1;
	int innerJEnd = grid.jEnd - 1;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

	int index, leftIndex, rightIndex, bottomIndex, topIndex;
#pragma omp parallel for private(index, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int j = innerJStart; j <= innerJEnd; j++)
	{
		for (int i = innerIStart; i <= innerIEnd; i++)
		{
			index = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			leftIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart - 1)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			rightIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart + 1)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			bottomIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart - 1)*innerIRes;
			topIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart + 1)*innerIRes;
			// Boundary condition.
			if (j == innerJStart)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
			}
			else if (j > innerJStart && j < innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;

				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;

				}
				else
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
			}
			else if (j == innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + 2 * dt*grid.oneOverdx2 + 2 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
			}
		}
	}
}

inline void MovingInterface::GenerateLinearSystem1(VectorND<double>& vectorB, const double & scaling)
{
	int innerIStart = grid.iStart + 1;
	int innerIEnd = grid.iEnd - 1;
	int innerJStart = grid.jStart + 1;
	int innerJEnd = grid.jEnd - 1;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

	Array2D<double> term(innerIStart, innerIRes, innerJStart, innerJRes);
	SurfactantNormalTerm(Surfactant, term);
	int index;
#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;
			vectorB(index) = Surfactant(i, j);// +dt*term(i, j);

			if (i == innerIStart)
			{
				vectorB(index) += dt*grid.dx2*ExactSurfactant(1, grid(i - 1, j).x, grid(i - 1, j).y, totalT);
			}
			if (i == innerIEnd)
			{
				vectorB(index) += dt*grid.dx2*ExactSurfactant(1, grid(i + 1, j).x, grid(i + 1, j).y, totalT);
			}
			if (j == innerJStart)
			{
				vectorB(index) += dt*grid.dy2*ExactSurfactant(1, grid(i, j - 1).x, grid(i, j - 1).y, totalT);

			}
			if (j == innerJEnd)
			{
				vectorB(index) += dt*grid.dy2*ExactSurfactant(1, grid(i, j + 1).x, grid(i, j + 1).y, totalT);
			}
			vectorB(index) *= scaling;
		}
	}
}

inline void MovingInterface::GenerateLinearSystem2(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int innerIStart = grid.iStart + 1;
	int innerIEnd = grid.iEnd - 1;
	int innerJStart = grid.jStart + 1;
	int innerJEnd = grid.jEnd - 1;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

	int index, leftIndex, rightIndex, bottomIndex, topIndex;
#pragma omp parallel for private(index, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int j = innerJStart; j <= innerJEnd; j++)
	{
		for (int i = innerIStart; i <= innerIEnd; i++)
		{
			index = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			leftIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart - 1)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			rightIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart + 1)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			bottomIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart - 1)*innerIRes;
			topIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart)
				+ (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart + 1)*innerIRes;
			// Boundary condition.
			if (j == innerJStart)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(topIndex) = scaling * -1 * dt*grid.oneOverdy2;
				}
			}
			else if (j > innerJStart && j < innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;

				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;

				}
				else
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
					matrixA(topIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
			}
			else if (j == innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + 1 * dt*grid.oneOverdx2 + 1 * dt*grid.oneOverdy2);
					matrixA(leftIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1.0 / 2.0 * dt*grid.oneOverdy2;
				}
			}
		}
	}
}

inline void MovingInterface::GenerateLinearSystem2(VectorND<double>& vectorB, const double & scaling)
{
	int innerIStart = grid.iStart + 1;
	int innerIEnd = grid.iEnd - 1;
	int innerJStart = grid.jStart + 1;
	int innerJEnd = grid.jEnd - 1;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

	Array2D<double> term(innerIStart, innerIRes, innerJStart, innerJRes);
	Array2D<double> termOld(innerIStart, innerIRes, innerJStart, innerJRes);
	SurfactantNormalTerm(Surfactant, term);
	SurfactantNormalTerm(SurfactantOld, termOld);

	int index;
#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;
			vectorB(index) = Surfactant(i, j) + dt / 2.0*(Surfactant.dxxPhi(i, j) + Surfactant.dyyPhi(i, j));
			//+ 3.0 / 2.0*dt*term(i, j) - 1.0 / 2.0*dt*termOld(i, j);

			if (i == innerIStart)
			{
				vectorB(index) += dt / 2.0*grid.dx2*ExactSurfactant(1, grid(i - 1, j).x, grid(i - 1, j).y, totalT);
			}
			if (i == innerIEnd)
			{
				vectorB(index) += dt / 2.0*grid.dx2*ExactSurfactant(1, grid(i + 1, j).x, grid(i + 1, j).y, totalT);
			}
			if (j == innerJStart)
			{
				vectorB(index) += dt / 2.0*grid.dy2*ExactSurfactant(1, grid(i, j - 1).x, grid(i, j - 1).y, totalT);

			}
			if (j == innerJEnd)
			{
				vectorB(index) += dt / 2.0*grid.dy2*ExactSurfactant(1, grid(i, j + 1).x, grid(i, j + 1).y, totalT);
			}
			vectorB(index) *= scaling;
		}
	}
}


