#pragma once

#include "AdvectionMethod2D.h"

class MovingInterface
{
public:
	MovingInterface();
	~MovingInterface();

	Grid2D grid;

	FD U; // x velocity
	FD V; // y velocity
	FD UOld; // x velocity
	FD VOld; // y velocity
	FD Surfactant;
	FD SurfactantOld;

	LS levelSet;
	LS levelSetOld;

	int ExamNum;


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

	inline void SurfactantNormalTerm(const FD& ipField, LS& ipLevelSet, Array2D<double>& term);
	inline double ExactSurfactant(const double& x, const double& y, const double& time);
	inline void GenerateLinearSystem1(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystem1(VectorND<double>& vectorB, const double & scaling);
	inline void GenerateLinearSystem2(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystem2(VectorND<double>& vectorB, const double & scaling);

	inline void SurfactantExtension(const LS& ipLevelSet, FD& ipField, const int& width);

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
	ExamNum = example;
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

		U = FD(grid);
		V = FD(grid);

		levelSet = LS(grid);
		double radius = 1.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}

		totalT = 1;
		Surfactant = FD(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);

			}
		}
		SurfactantOld = Surfactant;

		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		A = Array2D<double>(1, (grid.iRes - 2)*(grid.jRes - 2), 1, (grid.iRes - 2)*(grid.jRes - 2));
		dt = grid.dx / 4.0;
		maxIteration = 80;
	}

	if (example == 2)
	{
		cout << "*************************" << endl;
		cout << "       An Eulerian Formulation  " << endl;
		cout << " for Solving PDE along Moving Interface" << endl;
		cout << "          --JJ Xu, HK Zhao-- " << endl;
		cout << "               Example 2 " << endl;
		cout << "*************************" << endl;

		int gridSize = 40;
		grid = Grid2D(-3, 5, gridSize + 1, -3, 3, gridSize * 3. / 4. + 1);

		U = FD(grid);
		V = FD(grid);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				U(i, j) = 1;
				V(i, j) = 0;
			}
		}
		UOld = U;
		VOld = V;

		levelSet = LS(grid);
		double radius = 2.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}
		levelSetOld = levelSet;

		Surfactant = FD(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}
		SurfactantOld = Surfactant;

		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		A = Array2D<double>(1, (grid.iRes - 2)*(grid.jRes - 2), 1, (grid.iRes - 2)*(grid.jRes - 2));
		dt = grid.dx / 4.0;
		maxIteration = 80;
	}

	if (example == 3)
	{
		cout << "*************************" << endl;
		cout << "       An Eulerian Formulation  " << endl;
		cout << " for Solving PDE along Moving Interface" << endl;
		cout << "          --JJ Xu, HK Zhao-- " << endl;
		cout << "               Example 4 " << endl;
		cout << "*************************" << endl;

		int gridSize = 50;
		grid = Grid2D(-2, 8, gridSize + 1, -2, 2, gridSize * 4. / 10. + 1);

		U = FD(grid);
		V = FD(grid);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				U(i, j) = 1;
				V(i, j) = 0;
			}
		}
		UOld = U;
		VOld = V;

		levelSet = LS(grid);
		double radius = 1.0;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt(grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y) - radius;
			}
		}
		levelSetOld = levelSet;

		Surfactant = FD(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}
		SurfactantOld = Surfactant;

		AdvectionMethod2D<double>::alpha = 1.5*grid.dx;

		A = Array2D<double>(1, (grid.iRes - 2)*(grid.jRes - 2), 1, (grid.iRes - 2)*(grid.jRes - 2));
		dt = grid.dx / 4.0;
		maxIteration = 100;
	}
}

inline void MovingInterface::MovingInterfaceSolver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	const char* cmd;


	InitialCondition(example);

	grid.Variable();

	//LocalLS LLS(levelSet, 3.0*grid.dx);
	//LLS.levelSet.phi.Variable("phi");
	//LLS.T1.Variable("T1");
	//LLS.T2.Variable("T2");
	//LLS.T3.Variable("T3");

	GenerateLinearSystem1(A, 1);
	//A.Variable("A1");
	//// CG solver 1
	Acsr = CSR<double>(A);
	//// CG solver 2
	CGSolver::SparseA(A, a, row, col, nonzeroNum);

	FD exact(grid);


	//// Initial Surfactant
	Surfactant.Variable("Surfactant");
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1]), subplot(1,2,1),surf(X,Y,Surfactant);");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);

	//// Step 1
	cout << "********************************" << endl;
	cout << "       Iteration " << to_string(1) << " : Start" << endl;

	totalT += dt;
	OneStepSemiImplicit(example);

	cout << "       Iteration " << to_string(1) << " : End" << endl;
	cout << "********************************" << endl;

	Surfactant.Variable("Surfactant");
	MATLAB.Command("surf(X,Y,Surfactant);");
	str = string("title(['iteration : ', num2str(") + to_string(1) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);

	exact.Variable("SurfactantExact");
	str = string("l2(") + to_string(1) + string(")=norm(SurfactantNew-SurfactantExact)");
	cmd = str.c_str();
	MATLAB.Command(cmd);
	MATLAB.Command("subplot(1,2,2),surf(X,Y,SurfactantExact);");
	str = string("title(['iteration : ', num2str(") + to_string(1) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);


	GenerateLinearSystem2(A, 1);
	//A.Variable("A2");
	//// CG solver 1
	Acsr = CSR<double>(A);
	//// CG solver 2
	CGSolver::SparseA(A, a, row, col, nonzeroNum);


	for (int i = 2; i <= maxIteration; i++)
	{
		cout << "********************************" << endl;
		cout << "       Iteration " << to_string(i) << " : Start" << endl;

		totalT += dt;
		//OneStepSemiImplicit(example);
		TwoStepSemiImplicit(example);

		cout << "       Iteration " << to_string(i) << " : End" << endl;
		cout << "********************************" << endl;

		Surfactant.Variable("Surfactant");
		MATLAB.Command("subplot(1,2,1),surf(X,Y,Surfactant);");
		str = string("title(['iteration : ', num2str(") + to_string(i) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				exact(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			}
		}

		exact.Variable("SurfactantExact");
		str = string("l2(") + to_string(i) + string(")=norm(Surfactant-SurfactantExact)");
		cmd = str.c_str();
		MATLAB.Command(cmd);
		MATLAB.Command("subplot(1,2,2),surf(X,Y,SurfactantExact);");
		str = string("title(['iteration : ', num2str(") + to_string(i) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
	}
	//levelSet.phi.Variable("levelSet");
	//SurfactantExtension(levelSet, Surfactant, 10);

}

inline void MovingInterface::OneStepSemiImplicit(const int&example)
{
	//// Copy Old Data.
	//SurfactantOld = Surfactant;

	//// Linear Equation
	VectorND<double> vectorB((grid.iRes - 2)*(grid.jRes - 2));
	GenerateLinearSystem1(vectorB, 1);
	//vectorB.Variable("vectorB");

	VectorND<double> tempSur((grid.iRes - 2)*(grid.jRes - 2));

	int solver = 2;
	if (solver == 1)
	{
		tempSur = CGSolver::SolverCSR(Acsr, vectorB, grid.dx*grid.dy);
	}
	else if (solver == 2)
	{
		CGSolver::SolverSparse(A.iRes, a, row, col, vectorB, tempSur);
	}
	//tempSur.Variable("sur");

	int index;
#pragma omp parallel for private(index)
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			index = (i - (grid.iStart + 1)) + (j - (grid.jStart + 1))*(grid.iRes - 2);
			//if (grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y < 0.8*0.8)
			//{
			//	Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			//}
			//else
			{
				Surfactant(i, j) = tempSur(index);
			}
		}
	}

	//// Dirichlet Boundary Condition from the Exact Solution

#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		Surfactant(i, grid.jStart) = ExactSurfactant(grid(i, grid.jStart).x, grid(i, grid.jStart).y, totalT);
		Surfactant(i, grid.jEnd) = ExactSurfactant(grid(i, grid.jEnd).x, grid(i, grid.jEnd).y, totalT);
	}
#pragma omp parallel for
	for (int j = grid.jStart; j <= grid.jEnd; j++)
	{
		Surfactant(grid.iStart, j) = ExactSurfactant(grid(grid.iStart, j).x, grid(grid.iStart, j).y, totalT);
		Surfactant(grid.iEnd, j) = ExactSurfactant(grid(grid.iEnd, j).x, grid(grid.iEnd, j).y, totalT);
	}



	levelSetOld = levelSet;
	AdvectionMethod2D<double>::LSPropagatingTVDRK3(levelSet, U, V, dt);
	levelSet.phi.Variable("phi");
}

inline void MovingInterface::TwoStepSemiImplicit(const int&example)
{
	//// Copy Old Data.
	SurfactantOld = Surfactant;

	//// Linear Equation
	VectorND<double> vectorB((grid.iRes - 2)*(grid.jRes - 2));
	GenerateLinearSystem2(vectorB, 1);
	//vectorB.Variable("vectorB");

	VectorND<double> tempSur((grid.iRes - 2)*(grid.jRes - 2));

	int solver = 2;
	if (solver == 1)
	{
		tempSur = CGSolver::SolverCSR(Acsr, vectorB, grid.dx*grid.dy);
	}
	else if (solver == 2)
	{
		CGSolver::SolverSparse(A.iRes, a, row, col, vectorB, tempSur);
	}
	//tempSur.Variable("sur");

	int index;
#pragma omp parallel for private(index)
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			index = (i - (grid.iStart + 1)) + (j - (grid.jStart + 1))*(grid.iRes - 2);
			//if (grid(i, j).x*grid(i, j).x + grid(i, j).y*grid(i, j).y < 0.8*0.8)
			//{
			//	Surfactant(i, j) = ExactSurfactant(grid(i, j).x, grid(i, j).y, totalT);
			//}
			//else
			{
				Surfactant(i, j) = tempSur(index);
			}
		}
	}

	//// Dirichlet Boundary Condition from the Exact Solution
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		Surfactant(i, grid.jStart) = ExactSurfactant(grid(i, grid.jStart).x, grid(i, grid.jStart).y, totalT);
		Surfactant(i, grid.jEnd) = ExactSurfactant(grid(i, grid.jEnd).x, grid(i, grid.jEnd).y, totalT);
	}
#pragma omp parallel for
	for (int j = grid.jStart; j <= grid.jEnd; j++)
	{
		Surfactant(grid.iStart, j) = ExactSurfactant(grid(grid.iStart, j).x, grid(grid.iStart, j).y, totalT);
		Surfactant(grid.iEnd, j) = ExactSurfactant(grid(grid.iEnd, j).x, grid(grid.iEnd, j).y, totalT);
	}

	levelSetOld = levelSet;
	AdvectionMethod2D<double>::LSPropagatingTVDRK3(levelSet, U, V, dt);
	//levelSet.phi.Variable("phi");
}

inline void MovingInterface::SurfactantNormalTerm(const FD& ipField, LS& ipLevelSet, Array2D<double>& term)
{
	ipLevelSet.ComputeMeanCurvature();
	ipLevelSet.ComputeNormal();

	FV gradientF(FD::Gradient(ipField));
	//FV Normal(grid);
	VT normal;
	Array2D<double> Hessian(2, 2);

	FD wenoDxMinus(grid);
	FD wenoDxPlus(grid);
	FD wenoDyMinus(grid);
	FD wenoDyPlus(grid);
	AdvectionMethod2D<double>::WENO5thDerivation(ipField, wenoDxMinus, wenoDxPlus, wenoDyMinus, wenoDyPlus);

	FV gradientU(FD::Gradient(U));
	FV gradientV(FD::Gradient(V));

	double curvatureThreshold = 3.0;
	double curvature;
#pragma omp parallel for private(normal, Hessian, curvature)
	for (int i = term.iStart; i <= term.iEnd; i++)
	{
		for (int j = term.jStart; j <= term.jEnd; j++)
		{
			normal = ipLevelSet.normal(i, j);
			//Normal(i, j) = normal;
			Hessian = Surfactant.Hessian(i, j);
			curvature = -ipLevelSet.meanCurvature(i, j); //// LEVELSET.MEANCURVATURE has a nagative sign. so mutiple -1.
			term(i, j) = 0;
			if (abs(curvature) < curvatureThreshold)
			{
				term(i, j) += -curvature*dotProduct(normal, gradientF(i, j));
			}
			else
			{
				term(i, j) += -AdvectionMethod2D<double>::sign(curvature)*curvatureThreshold*dotProduct(normal, gradientF(i, j));
			}
			term(i, j) += -normal(0)*(Hessian(0, 0)*normal(0) + Hessian(0, 1)*normal(1));
			term(i, j) += -normal(1)*(Hessian(1, 0)*normal(0) + Hessian(1, 1)*normal(1));

			//// Upwind WENO
			term(i, j) += -(AdvectionMethod2D<double>::Plus(U(i, j))*wenoDxMinus(i, j) + AdvectionMethod2D<double>::Minus(U(i, j))*wenoDxPlus(i, j));
			term(i, j) += -(AdvectionMethod2D<double>::Plus(V(i, j))*wenoDyMinus(i, j) + AdvectionMethod2D<double>::Minus(V(i, j))*wenoDyPlus(i, j));

			////
			term(i, j) += normal(0)*(gradientU(i, j).x*normal(0) + gradientU(i, j).y*normal(1))*ipField(i, j);
			term(i, j) += normal(1)*(gradientV(i, j).x*normal(0) + gradientV(i, j).y*normal(1))*ipField(i, j);
		}
	}
	//term.Variable("addedTerm");
	//ArrayVec2DVariable("normal", Normal.dataArray);
	//ArrayVec2DVariable("gradient", gradient.dataArray);
}

inline double MovingInterface::ExactSurfactant(const double & x, const double & y, const double & time)
{
	if (ExamNum == 1)
	{
		return exp(-time / (x*x + y*y + DBL_EPSILON))*(y / sqrt(x*x + y*y + DBL_EPSILON)) + 2;
	}
	if (ExamNum == 2 || ExamNum == 3)
	{
		return exp(-time / ((x - totalT)*(x - totalT) + y*y + DBL_EPSILON))*(y / sqrt((x - totalT)*(x - totalT) + y*y + DBL_EPSILON)) + 2;
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
	SurfactantNormalTerm(Surfactant, levelSet, term);
	//term.Variable("addedTerm");
	int index;
#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;
			vectorB(index) = Surfactant(i, j) + dt*term(i, j);

			if (i == innerIStart)
			{
				vectorB(index) += dt*grid.oneOverdx2*ExactSurfactant(grid(i - 1, j).x, grid(i - 1, j).y, totalT);
			}
			if (i == innerIEnd)
			{
				vectorB(index) += dt*grid.oneOverdx2*ExactSurfactant(grid(i + 1, j).x, grid(i + 1, j).y, totalT);
			}
			if (j == innerJStart)
			{
				vectorB(index) += dt*grid.oneOverdy2*ExactSurfactant(grid(i, j - 1).x, grid(i, j - 1).y, totalT);

			}
			if (j == innerJEnd)
			{
				vectorB(index) += dt*grid.oneOverdy2*ExactSurfactant(grid(i, j + 1).x, grid(i, j + 1).y, totalT);
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
			//// Boundary condition.
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
					matrixA(topIndex) = scaling * -1 / 2.0 * dt*grid.oneOverdy2;
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
	SurfactantNormalTerm(Surfactant, levelSet, term);
	SurfactantNormalTerm(SurfactantOld, levelSetOld, termOld);

	int index;
#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;
			vectorB(index) = Surfactant(i, j) + dt / 2.0*(Surfactant.dxxPhi(i, j) + Surfactant.dyyPhi(i, j))
				+ 3.0 / 2.0*dt*term(i, j) - 1.0 / 2.0*dt*termOld(i, j);

			if (i == innerIStart)
			{
				vectorB(index) += dt / 2.0*grid.oneOverdx2*ExactSurfactant(grid(i - 1, j).x, grid(i - 1, j).y, totalT);
			}
			if (i == innerIEnd)
			{
				vectorB(index) += dt / 2.0*grid.oneOverdx2*ExactSurfactant(grid(i + 1, j).x, grid(i + 1, j).y, totalT);
			}
			if (j == innerJStart)
			{
				vectorB(index) += dt / 2.0*grid.oneOverdy2*ExactSurfactant(grid(i, j - 1).x, grid(i, j - 1).y, totalT);

			}
			if (j == innerJEnd)
			{
				vectorB(index) += dt / 2.0*grid.oneOverdy2*ExactSurfactant(grid(i, j + 1).x, grid(i, j + 1).y, totalT);
			}
			vectorB(index) *= scaling;
		}
	}
}

inline void MovingInterface::SurfactantExtension(const LS & ipLevelSet, FD& ipField, const int& width)
{
	VT normal;
	FD k1(ipField.grid);
	FD k2(ipField.grid);
	FD k3(ipField.grid);

	FD wenoXMinus(ipField.grid);
	FD wenoXPlus(ipField.grid);
	FD wenoYMinus(ipField.grid);
	FD wenoYPlus(ipField.grid);

	Array2D<VT> Normal(ipField.grid);


	ipField.Variable("extension");
	MATLAB.Command("subplot(1,2,2),contour(X,Y,levelSet,[0 0],'r');");
	MATLAB.Command("subplot(1,2,2), contour(X,Y,extension),hold on, contour(X,Y,levelSet,[0 0]), hold off;");
	string str = string("title(['iteration : ', num2str(") + to_string(0) + string(")]);");
	MATLAB.Command(str.c_str());
	for (int i = 1; i <= width; i++)
	{
		FD originField = ipField;
		AdvectionMethod2D<double>::WENO5thDerivation(ipField, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(normal)
		for (int i = ipField.grid.iStart; i <= ipField.grid.iEnd; i++)
		{
			for (int j = ipField.grid.jStart; j <= ipField.grid.jEnd; j++)
			{
				Normal(i, j) = ipLevelSet.normal(i, j);
				normal = ipLevelSet.normal(i, j);
				k1(i, j) = 0;
				k1(i, j) += (AdvectionMethod2D<double>::Plus(normal(0))*wenoXMinus(i, j) + AdvectionMethod2D<double>::Minus(normal(0))*wenoXPlus(i, j));
				k1(i, j) += (AdvectionMethod2D<double>::Plus(normal(1))*wenoYMinus(i, j) + AdvectionMethod2D<double>::Minus(normal(1))*wenoYPlus(i, j));
				k1(i, j) *= -dt*ipLevelSet(i, j) / sqrt(ipLevelSet(i, j)*ipLevelSet(i, j) + grid.dx2);
				ipField(i, j) = originField(i, j) + k1(i, j);
			}
		}

		AdvectionMethod2D<double>::WENO5thDerivation(ipField, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(normal)
		for (int i = ipField.grid.iStart; i <= ipField.grid.iEnd; i++)
		{
			for (int j = ipField.grid.jStart; j <= ipField.grid.jEnd; j++)
			{
				normal = ipLevelSet.normal(i, j);
				k2(i, j) = 0;
				k2(i, j) += (AdvectionMethod2D<double>::Plus(normal(0))*wenoXMinus(i, j) + AdvectionMethod2D<double>::Minus(normal(0))*wenoXPlus(i, j));
				k2(i, j) += (AdvectionMethod2D<double>::Plus(normal(1))*wenoYMinus(i, j) + AdvectionMethod2D<double>::Minus(normal(1))*wenoYPlus(i, j));
				k2(i, j) *= -dt*ipLevelSet(i, j) / sqrt(ipLevelSet(i, j)*ipLevelSet(i, j) + grid.dx2);
				ipField(i, j) = 3.0 / 4.0*originField(i, j) + 1.0 / 4.0*(ipField(i, j) + k2(i, j));
			}
		}

		AdvectionMethod2D<double>::WENO5thDerivation(ipField, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for private(normal)
		for (int i = ipField.grid.iStart; i <= ipField.grid.iEnd; i++)
		{
			for (int j = ipField.grid.jStart; j <= ipField.grid.jEnd; j++)
			{
				normal = ipLevelSet.normal(i, j);
				k3(i, j) = 0;
				k3(i, j) += (AdvectionMethod2D<double>::Plus(normal(0))*wenoXMinus(i, j) + AdvectionMethod2D<double>::Minus(normal(0))*wenoXPlus(i, j));
				k3(i, j) += (AdvectionMethod2D<double>::Plus(normal(1))*wenoYMinus(i, j) + AdvectionMethod2D<double>::Minus(normal(1))*wenoYPlus(i, j));
				k3(i, j) *= -dt*ipLevelSet(i, j) / sqrt(ipLevelSet(i, j)*ipLevelSet(i, j) + grid.dx2);
				ipField(i, j) = 1.0 / 3.0*originField(i, j) + 2.0 / 3.0*(ipField(i, j) + k3(i, j));
			}
		}
		ipField.Variable("extension");
		MATLAB.Command("subplot(1,2,1), surf(X,Y,extension);");
		MATLAB.Command("subplot(1,2,2), contour(X,Y,extension),hold on, contour(X,Y,levelSet,[0 0]), hold off;");
		string str = string("title(['iteration : ', num2str(") + to_string(i) + string(")]);");
		MATLAB.Command(str.c_str());
	}


}


