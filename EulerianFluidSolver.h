#pragma once


#include "AdvectionMethod2D.h"


class EulerianFluidSolver2D
{
public:
	Grid2D gridP;
	Grid2D gridPinner;
	Grid2D gridU;
	Grid2D gridUinner;
	Grid2D gridV;
	Grid2D gridVinner;
	
	FD originU;
	FD originV;
	FD U; // x velocity
	FD Ustar;
	FD V; // y velocity
	FD Vstar;
	FD P; // Pressure
	FD Rho;
	FD Mu;
	FD Nu; // Viscosity

	FD advectionU;
	FD advectionV;

	FD diffusionU;
	FD diffusionV;

	LS levelSet;

	int CGsolverNum;
	Array2D<double>poissonMatrix;
	// CG solver 1
	CSR<double> poissonCSR;
	// CG solver 2
	VectorND<double> a;
	VectorND<int> row;
	VectorND<int> col;
	int nonzeroNum;
	VectorND<double> vectorB;
	VectorND<double> tempP;

	int iteration;
	double reynoldNum;
	double dt;
	
	double cflCondition;

	int accuracyOrder;

	int ghostWidth;

	int maxIteration;
	int writeOutputIteration;

	EulerianFluidSolver2D();
	~EulerianFluidSolver2D();

	inline void InitialCondition(const int& example);

	// Chorin's Projection Method : 1st order accuracy.
	inline void FluidSolver(const int& example);

	inline void GenerateLinearSystem(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystem(VectorND<double>& vectorB, const double & scaling);
	inline void TVDRK3TimeAdvection();
	inline void EulerMethod();
	inline void AdvectionTerm(FD& U, FD& V, FD& TermU, FD& TermV);
	inline void DiffusionTerm(const FD& U, const FD& V, FD& TermU, FD& TermV);

	inline double AdaptiveTimeStep(const FD& velocity1, const FD& velocity2);



	// Second Order Projection Method
	FD oldU;
	FD oldV;
	FD gradientPx;
	FD gradientPy;
	FD Phi;
	FD Phixxyy;
	FD advectionU1;
	FD advectionV1;
	FD advectionU2;
	FD advectionV2;

	Array2D<double> UCNMatrix;
	// CG solver 1
	CSR<double> UCN_CSR;
	// CG solver 2
	VectorND<double> Ua;
	VectorND<int> Urow;
	VectorND<int> Ucol;
	int UnonzeroNum;
	VectorND<double> Ub;
	VectorND<double> tempU;

	Array2D<double> VCNMatrix;
	// CG solver 1
	CSR<double> VCN_CSR;
	// CG solver 2
	VectorND<double> Va;
	VectorND<int> Vrow;
	VectorND<int> Vcol;
	int VnonzeroNum;
	VectorND<double> Vb;
	VectorND<double> tempV;

	Array2D<double> PhiCNMatrix;
	// CG solver 1
	CSR<double> PhiCN_CSR;
	// CG solver 2
	VectorND<double> Phia;
	VectorND<int> Phirow;
	VectorND<int> Phicol;
	int PhinonzeroNum;
	VectorND<double> Phib;
	VectorND<double> tempPhi;

	inline void NSSolver2ndOrder(const int& example);
	inline void EulerMethod2ndOrder();
	inline void EulerMethod2ndOrder1();
	inline void EulerMethod2ndOrder2();
	inline void EulerMethod2ndOrder3();
	inline void EulerMethod2ndOrder1stIteration1();
	inline void GenerateLinearSystemUV(Array2D<double>& matrixA, const Grid2D& ipGrid, const double & scaling);
	inline void GenerateLinearSystemUV(VectorND<double>& vectorB, const FD& vel, const FD& gradP, const FD& advec, const Grid2D& ipGrid, const double & scaling);
	inline void GenerateLinearSystemPhi(Array2D<double>& matrixA, const double & scaling);
	inline void GenerateLinearSystemPhi(VectorND<double>& vectorB, const double & scaling);

private:

};




EulerianFluidSolver2D::EulerianFluidSolver2D()
{
}

EulerianFluidSolver2D::~EulerianFluidSolver2D()
{
}

inline void EulerianFluidSolver2D::InitialCondition(const int & example)
{
	if (example==1)
	{
		cout << "*************************" << endl;
		cout << "    Navier-Stokes equation" << endl;
		cout << "    Chorin's Projection Method" << endl;
		cout << "    Cavity Flow" << endl;
		cout << "*************************" << endl;
		cout << endl;

		int numP = 51;
		double ddx = 0.01;
		gridP = Grid2D(0, ddx*double(numP-1), numP, 0, ddx*double(numP - 1), numP);
		gridPinner = Grid2D(gridP.xMin + gridP.dx, gridP.xMax - gridP.dx, 1, gridP.iRes - 2,
			gridP.yMin + gridP.dy, gridP.yMax - gridP.dy, 1, gridP.jRes - 2);
		gridU = Grid2D(gridP.xMin - gridP.dx/2, gridP.xMax + gridP.dx/2, gridP.iRes + 1,
			gridP.yMin, gridP.yMax, gridP.jRes);
		gridUinner = Grid2D(gridU.xMin + gridU.dx, gridU.xMax - gridU.dx, 1, gridU.iRes - 2,
			gridU.yMin + gridU.dy, gridU.yMax - gridU.dy, 1, gridU.jRes - 2);
		gridV = Grid2D(gridP.xMin, gridP.xMax, gridP.iRes, 
			gridP.yMin - gridP.dy / 2, gridP.yMax + gridP.dy / 2, gridP.jRes + 1);
		gridVinner = Grid2D(gridV.xMin + gridV.dx, gridV.xMax - gridV.dx, 1, gridV.iRes - 2,
			gridV.yMin + gridV.dy, gridV.yMax - gridV.dy, 1, gridV.jRes - 2);
		P = FD(gridP);
		U = FD(gridU);
		V = FD(gridV);

		Rho = FD(gridP);

		// initial condition
#pragma omp parallel for
		for (int i = gridU.iStart; i <= gridU.iEnd; i++)
		{
			U(i, gridU.jEnd) = 1;
		}
		Ustar = U;
		Vstar = V;

#pragma omp parallel for
		for (int i = Rho.iStart; i <= Rho.iEnd; i++)
		{
			for (int j = Rho.jStart; j <= Rho.jEnd; j++)
			{
				Rho(i, j) = 1;
			}
		}
		advectionU = FD(gridUinner);
		advectionV = FD(gridVinner);

		diffusionU = FD(gridUinner);
		diffusionV = FD(gridVinner);

		vectorB = VectorND<double>(gridPinner.iRes*gridPinner.jRes);
		tempP = VectorND<double>(gridPinner.iRes*gridPinner.jRes);
		poissonMatrix = Array2D<double>(1, gridPinner.iRes*gridPinner.jRes, 1, gridPinner.iRes*gridPinner.jRes);

		accuracyOrder = 1;
		reynoldNum = 1000;
		cflCondition = 0.5;

		dt = cflCondition*gridP.dx;
		maxIteration = 2000;
		writeOutputIteration = 10;


	}

	if (example==2)
	{
		cout << "*************************" << endl;
		cout << "    Navier-Stokes equation" << endl;
		cout << "    2nd Order Projection Method" << endl;
		cout << "    with a Lagged Pressure Term" << endl;
		cout << "    -Cavity Flow-" << endl;
		cout << "*************************" << endl;
		cout << endl;

		int numP = 51;
		double ddx = 0.01;
		gridP = Grid2D(0, ddx*double(numP - 1), numP, 0, ddx*double(numP - 1), numP);
		gridPinner = Grid2D(gridP.xMin + gridP.dx, gridP.xMax - gridP.dx, 1, gridP.iRes - 2,
			gridP.yMin + gridP.dy, gridP.yMax - gridP.dy, 1, gridP.jRes - 2);
		gridU = Grid2D(gridP.xMin - gridP.dx / 2, gridP.xMax + gridP.dx / 2, gridP.iRes + 1,
			gridP.yMin, gridP.yMax, gridP.jRes);
		gridUinner = Grid2D(gridU.xMin + gridU.dx, gridU.xMax - gridU.dx, 1, gridU.iRes - 2,
			gridU.yMin + gridU.dy, gridU.yMax - gridU.dy, 1, gridU.jRes - 2);
		gridV = Grid2D(gridP.xMin, gridP.xMax, gridP.iRes,
			gridP.yMin - gridP.dy / 2, gridP.yMax + gridP.dy / 2, gridP.jRes + 1);
		gridVinner = Grid2D(gridV.xMin + gridV.dx, gridV.xMax - gridV.dx, 1, gridV.iRes - 2,
			gridV.yMin + gridV.dy, gridV.yMax - gridV.dy, 1, gridV.jRes - 2);
		//P = FD(gridP);
		Phi = FD(gridP);
		Phixxyy = FD(gridP);
		U = FD(gridU);
		V = FD(gridV);

		gradientPx = FD(gridUinner);
		gradientPy = FD(gridVinner);

		// initial condition
#pragma omp parallel for
		for (int i = gridU.iStart; i <= gridU.iEnd; i++)
		{
			U(i, gridU.jEnd) = 1;
		}
		Ustar = U;
		Vstar = V;

		oldU = U;
		oldV = V;

		advectionU = FD(gridUinner);
		advectionV = FD(gridVinner);
		advectionU1 = FD(gridUinner);
		advectionV1 = FD(gridVinner);
		advectionU2 = FD(gridUinner);
		advectionV2 = FD(gridVinner);
		diffusionU = FD(gridUinner);
		diffusionV = FD(gridVinner);
		
		Ub = VectorND<double>(gridUinner.iRes*gridUinner.jRes);
		tempU = VectorND<double>(gridUinner.iRes*gridUinner.jRes);
		UCNMatrix = Array2D<double>(1, gridUinner.iRes*gridUinner.jRes, 1, gridUinner.iRes*gridUinner.jRes);
		
		Vb = VectorND<double>(gridVinner.iRes*gridVinner.jRes);
		tempV = VectorND<double>(gridVinner.iRes*gridVinner.jRes);
		VCNMatrix = Array2D<double>(1, gridVinner.iRes*gridVinner.jRes, 1, gridVinner.iRes*gridVinner.jRes);

		Phib = VectorND<double>(gridPinner.iRes*gridPinner.jRes);
		tempPhi = VectorND<double>(gridPinner.iRes*gridPinner.jRes);
		PhiCNMatrix = Array2D<double>(1, gridPinner.iRes*gridPinner.jRes, 1, gridPinner.iRes*gridPinner.jRes);


		accuracyOrder = 2;
		reynoldNum = 100;
		cflCondition = 0.5;

		maxIteration = 2000;
		writeOutputIteration = 10;

		dt = cflCondition*gridP.dx;
		CGsolverNum = 2;
	}
}

inline void EulerianFluidSolver2D::FluidSolver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	const char* cmd;

	InitialCondition(example);
	gridP.Variable("Xp", "Yp");
	gridU.Variable("Xu", "Yu");
	gridV.Variable("Xv", "Yv");
	U.Variable("U");
	V.Variable("V");
	P.Variable("P");

	GenerateLinearSystem(poissonMatrix, -gridP.dx*gridP.dx);
	poissonMatrix.Variable("poisson");

	//int solver = 2;
	//// CG solver 1
	poissonCSR = CSR<double>(poissonMatrix);

	//// CG solver 2
	CGSolver::SparseA(poissonMatrix, a, row, col, nonzeroNum);

	P.Variable("P");
	U.Variable("U");
	V.Variable("V");
	str = string("quiver(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),axis([Xp(1)-(Xp(end)-Xp(1))/10 Xp(end)+(Xp(end)-Xp(1))/10 Yp(1)-(Yp(end)-Yp(1))/10 Yp(end)+(Yp(end)-Yp(1))/10]);");
	str = str + string("hold on,streamline(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,-100:0.1:100,-100:0.1:100),hold off;");
	MATLAB.Command(str.c_str());
	str = string("title(['iteration : ', num2str(") + to_string(0) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);
	double totalT = 0.0;
	for (int i = 1; i <= maxIteration; i++)
	{
		cout << endl;
		cout << "********************************" << endl;
		cout << "       Iteration " << to_string(i) << " : Start" << endl;
		TVDRK3TimeAdvection();
		totalT += dt;
		cout << "       Iteration " << to_string(i) << " : End" << endl;
		cout << "********************************" << endl;
		P.Variable("P");
		U.Variable("U");
		V.Variable("V");
		MATLAB.Command("axis([Xp(1)-(Xp(end)-Xp(1))/10 Xp(end)+(Xp(end)-Xp(1))/10 Yp(1)-(Yp(end)-Yp(1))/10 Yp(end)+(Yp(end)-Yp(1))/10])");

		str = string("quiver(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),axis([Xp(1)-(Xp(end)-Xp(1))/10 Xp(end)+(Xp(end)-Xp(1))/10 Yp(1)-(Yp(end)-Yp(1))/10 Yp(end)+(Yp(end)-Yp(1))/10]);");
		str = str + string("hold on,streamline(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,-100:0.1:100,-100:0.1:100),hold off;");
		MATLAB.Command(str.c_str());
		str = string("title(['iteration : ', num2str(") + to_string(i) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
		MATLAB.Command("divU =U(:,2:end)-U(:,1:end-1),divV =V(2:end,:)-V(1:end-1,:);div=divU+divV;");

		if (writeFile && i%writeOutputIteration == 0)
		{
			fileName = "pressure" + to_string(i);
			P.WriteFile(fileName);
			fileName = "xVelocity" + to_string(i);
			U.WriteFile(fileName);
			fileName = "yVelocity" + to_string(i);
			V.WriteFile(fileName);
		}
	}
}

inline void EulerianFluidSolver2D::GenerateLinearSystem(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int innerIStart = gridPinner.iStart;
	int innerIEnd = gridPinner.iEnd;
	int innerJStart = gridPinner.jStart;
	int innerJEnd = gridPinner.jEnd;
	int innerIRes = gridPinner.iRes;
	int innerJRes = gridPinner.jRes;

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
					matrixA(index) = scaling*(-1 * gridP.oneOverdx2 - 1 * gridP.oneOverdy2);
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(-1 * gridP.oneOverdx2 - 1 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 1 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				
			}
			else if (j>innerJStart && j<innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(-1 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;

				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(-1 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;

				}
				else
				{
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
					matrixA(topIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				
			}
			else if (j == innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 1 * gridP.oneOverdy2);
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(-1 * gridP.oneOverdx2 - 1 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(-2 * gridP.oneOverdx2 - 1 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * 1 * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * 1 * gridP.oneOverdy2;
				}
				
			}
		}
	}
}

inline void EulerianFluidSolver2D::GenerateLinearSystem(VectorND<double>& vectorB, const double & scaling)
{
	int innerIStart = gridPinner.iStart;
	int innerIEnd = gridPinner.iEnd;
	int innerJStart = gridPinner.jStart;
	int innerJEnd = gridPinner.jEnd;
	int innerIRes = gridPinner.iRes;
	int innerJRes = gridPinner.jRes;

	int index;
	P(gridP.iStart, gridP.jEnd - 1) = 1;
#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;

			vectorB(index) = Rho(i, j) / dt*((U(i + 1, j) - U(i, j))*gridU.oneOverdx 
				+ (V(i, j + 1) - V(i, j))*gridV.oneOverdy);
			if (i == innerIStart && j==innerJEnd)
			{
				vectorB(index) += -P(i - 1, j)*P.oneOverdx2;
			}
			//cout << endl;
			//cout << "(i,j) = (" << i << "," << j << ")" << endl;
			//cout << "index = " << index << endl;
			//cout << "U " << U(i + 1, j) << " " << U(i, j) << endl;
			//cout << "V " << V(i, j + 1) << " " << V(i, j) << endl;
			//cout << "B" << vectorB(index) << endl;

			//if (i == innerIStart)
			//{
			//	vectorB(index) += -P(i - 1, j)*P.oneOverdx2;
			//}
			//if (i == innerIEnd)
			//{
			//	vectorB(index) += -P(i + 1, j)*P.oneOverdx2;
			//}
			//if (j == innerJStart)
			//{
			//	vectorB(index) += -P(i, j - 1)*P.oneOverdy2;

			//}
			//if (j == innerJEnd)
			//{
			//	vectorB(index) += -P(i, j + 1)*P.oneOverdy2;
			//}
			//vectorB.Variable("vecB");
			vectorB(index) *= scaling;
			//if (abs(vectorB(index))>0)
			//{
			//	cout << i << " " << j << " " << index << " " << vectorB(index) << endl;
			//}
		}
	}

}

inline void EulerianFluidSolver2D::TVDRK3TimeAdvection()
{
	originU = U;
	originV = V;

	//dt = AdaptiveTimeStep(U, V);
	
	/////////////////
	//// Step 1  ////
	/////////////////
	if (accuracyOrder == 1)
	{
		EulerMethod();
	}
	else if (accuracyOrder == 2)
	{
		EulerMethod2ndOrder();
	}
	//U.Variable("U1");
	//V.Variable("V1");
	//MATLAB.Command("quiver(Xp,Yp,U1(:,1:end-1),V1(1:end-1,:))");

	/////////////////
	//// Step 2  ////
	/////////////////
	if (accuracyOrder == 1)
	{
		EulerMethod();
	}
	else if (accuracyOrder == 2)
	{
		EulerMethod2ndOrder();
	}
	//U.Variable("U21");
	//V.Variable("V21");
	//MATLAB.Command("quiver(Xp,Yp,U21(:,1:end-1),V21(1:end-1,:))");
#pragma omp parallel for
	for (int i = gridU.iStart; i <= gridU.iEnd; i++)
	{
		for (int j = gridU.jStart; j <= gridU.jEnd; j++)
		{
			U(i, j) = 3. / 4. * originU(i, j) + 1. / 4. * U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = gridV.iStart; i <= gridV.iEnd; i++)
	{
		for (int j = gridV.jStart; j <= gridV.jEnd; j++)
		{
			V(i, j) = 3. / 4. * originV(i, j) + 1. / 4. * V(i, j);
		}
	}
	//U.Variable("U2");
	//V.Variable("V2");
	//MATLAB.Command("quiver(Xp,Yp,U22(:,1:end-1),V22(1:end-1,:))");
	
	/////////////////
	//// Step 3  ////
	/////////////////
	if (accuracyOrder == 1)
	{
		EulerMethod();
	}
	else if (accuracyOrder == 2)
	{
		EulerMethod2ndOrder();
	}
	//U.Variable("U31");
	//V.Variable("V31");
	//MATLAB.Command("quiver(Xp,Yp,U31(:,1:end-1),V31(1:end-1,:))");
#pragma omp parallel for
	for (int i = gridU.iStart; i <= gridU.iEnd; i++)
	{
		for (int j = gridU.jStart; j <= gridU.jEnd; j++)
		{
			U(i, j) = 1. / 3. * originU(i, j) + 2. / 3. * U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = gridV.iStart; i <= gridV.iEnd; i++)
	{
		for (int j = gridV.jStart; j <= gridV.jEnd; j++)
		{
			V(i, j) = 1. / 3. * originV(i, j) + 2. / 3. * V(i, j);
		}
	}
	//U.Variable("U3");
	//V.Variable("V3");
	//MATLAB.Command("quiver(Xp,Yp,U32(:,1:end-1),V32(1:end-1,:))");

}

inline void EulerianFluidSolver2D::EulerMethod()
{
	//// Boundary : Linear extension.
#pragma omp parallel for
	for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
	{
		U(gridU.iStart, j) = -U(gridU.iStart + 1, j);
		U(gridU.iEnd, j) = -U(gridU.iEnd - 1, j);
	}
#pragma omp parallel for
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		V(i, gridV.jStart) = -V(i, gridV.jStart + 1);
		V(i, gridV.jEnd) = -V(i, gridV.jEnd - 1);
	}

	////////////////////////////////////////////////
	////     Projection Method 1 : advection    ////
	////////////////////////////////////////////////
	
	Array2D<double>& K1U = U.K1;
	Array2D<double>& K1V = V.K1;

	AdvectionTerm(U, V, advectionU, advectionV);
	DiffusionTerm(U, V, diffusionU, diffusionV);
	//advectionU.Variable("advectionU");
	//diffusionU.Variable("diffusionU");
	//advectionV.Variable("advectionV");
	//diffusionV.Variable("diffusionV");

#pragma omp parallel for
	for (int i = K1U.iStart; i <= K1U.iEnd; i++)
	{
		for (int j = K1U.jStart; j <= K1U.jEnd; j++)
		{
			K1U(i, j) = dt*(-advectionU(i, j) + 1. / reynoldNum*diffusionU(i, j));
			U(i, j) = U(i, j) + K1U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = K1V.iStart; i <= K1V.iEnd; i++)
	{
		for (int j = K1V.jStart; j <= K1V.jEnd; j++)
		{
			K1V(i, j) = dt*(-advectionV(i, j) + 1. / reynoldNum*diffusionV(i, j));
			V(i, j) = V(i, j) + K1V(i, j);
		}
	}
	//K1U.Variable("k1u");
	//K1V.Variable("k1v");
	//// Boundary : Linear extension.
#pragma omp parallel for
	for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
	{
		U(gridU.iStart, j) = -U(gridU.iStart + 1, j);
		U(gridU.iEnd, j) = -U(gridU.iEnd - 1, j);
	}
#pragma omp parallel for
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		V(i, gridV.jStart) = -V(i, gridV.jStart + 1);
		V(i, gridV.jEnd) = -V(i, gridV.jEnd - 1);
	}


	U.Variable("Ustar");
	V.Variable("Vstar");
	MATLAB.Command("divUstar =Ustar(:,2:end)-Ustar(:,1:end-1),divVstar =Vstar(2:end,:)-Vstar(1:end-1,:);divstar=divUstar+divVstar;");
	//MATLAB.Command("quiver(Xp,Yp,Ustar(:,1:end-1),Vstar(1:end-1,:)");


	////////////////////////////////////////////////
	////     Projection Method 2 : Poisson Eq   ////
	////////////////////////////////////////////////
	GenerateLinearSystem(vectorB, -gridP.dx*gridP.dx);
	vectorB.Variable("vectorB");
	int solver = 2;
	if (solver == 1)
	{
		tempP = CGSolver::SolverCSR(poissonCSR, vectorB, gridP.dx*gridP.dy);
	}
	else if (solver == 2)
	{
		CGSolver::SolverSparse(poissonMatrix.iRes, a, row, col, vectorB, tempP);
	}
	//tempP.Variable("tempP");

	int index;
#pragma omp parallel for private(index)
	for (int i = gridPinner.iStart; i <= gridPinner.iEnd; i++)
	{
		for (int j = gridPinner.jStart; j <= gridPinner.jEnd; j++)
		{
			index = (i - gridPinner.iStart) + (j - gridPinner.jStart)*gridPinner.iRes;
			P(i, j) = tempP(index);
		}
	}
#pragma omp parallel for
	for (int i = gridP.iStart; i <= gridP.iEnd; i++)
	{
		P(i, P.jStart) = P(i, P.jStart + 1);
		P(i, P.jEnd) = P(i, P.jEnd - 1);
	}
#pragma omp parallel for
	for (int j = gridP.jStart; j <= gridP.jEnd; j++)
	{
		P(P.iStart, j) = P(P.iStart + 1, j);
		P(P.iEnd, j) = P(P.iEnd - 1, j);
	}
	P(gridP.iStart, gridP.jEnd - 1) = 1;
	P.Variable("P");


	//////////////////////////////////////////////
	////     Projection Method 3 : New U,V    ////
	//////////////////////////////////////////////
#pragma omp parallel for
	for (int i = gridUinner.iStart; i <= gridUinner.iEnd; i++)
	{
		for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
		{
			U(i, j) = U(i, j) - dt * 1. / Rho(i, j)*(P(i, j) - P(i - 1, j))*P.oneOverdx;
		}
	}
#pragma omp parallel for
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		for (int j = gridVinner.jStart; j <= gridVinner.jEnd; j++)
		{
			V(i, j) = V(i, j) - dt * 1. / Rho(i, j)*(P(i, j) - P(i, j - 1))*P.oneOverdy;
		}
	}
	//// Boundary : Linear extension.
#pragma omp parallel for
	for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
	{
		U(gridU.iStart, j) = -U(gridU.iStart + 1, j);
		U(gridU.iEnd, j) = -U(gridU.iEnd - 1, j);
	}
#pragma omp parallel for
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		V(i, gridV.jStart) = -V(i, gridV.jStart + 1);
		V(i, gridV.jEnd) = -V(i, gridV.jEnd - 1);
	}

	U.Variable("Unew");
	V.Variable("Vnew");
	MATLAB.Command("divUnew =Unew(:,2:end)-Unew(:,1:end-1),divVnew =Vnew(2:end,:)-Vnew(1:end-1,:);divnew=divUnew+divVnew;");
	//MATLAB.Command("quiver(Xp,Yp,U(1:end-1,:),V(:,1:end-1)");
}


inline void EulerianFluidSolver2D::AdvectionTerm(FD& U, FD& V, FD& TermU, FD& TermV)
{

	Array2D<double>& dUdxM = U.dfdxM;
	Array2D<double>& dUdxP = U.dfdxP;
	Array2D<double>& dUdyM = U.dfdyM;
	Array2D<double>& dUdyP = U.dfdyP;
	AdvectionMethod2D<double>::ENO3rdDerivation(U, dUdxM, dUdxP, dUdyM, dUdyP);
	Array2D<double>& dVdxM = V.dfdxM;
	Array2D<double>& dVdxP = V.dfdxP;
	Array2D<double>& dVdyM = V.dfdyM;
	Array2D<double>& dVdyP = V.dfdyP;
	AdvectionMethod2D<double>::ENO3rdDerivation(V, dVdxM, dVdxP, dVdyM, dVdyP);


	double Ux, Uy;
	double aveV;
#pragma omp parallel for private(aveV, Ux, Uy)
	for (int i = TermU.iStart; i <= TermU.iEnd; i++)
	{
		for (int j = TermU.jStart; j <= TermU.jEnd; j++)
		{
			aveV = (V(i - 1, j) + V(i, j) + V(i - 1, j + 1) + V(i, j + 1)) / 4.;
			if (U(i,j)>0)
			{
				Ux = dUdxM(i, j);
			}
			else
			{
				Ux = dUdxP(i, j);
			}

			if (aveV>0)
			{
				Uy = dUdyM(i, j);
			}
			else
			{
				Uy = dUdyP(i, j);
			}
			TermU(i, j) = U(i, j)*Ux + aveV*Uy;
		}
	}

	double aveU;
	double Vx, Vy;
#pragma omp parallel for private(aveU, Vx, Vy)
	for (int i = TermV.iStart; i <= TermV.iEnd; i++)
	{
		for (int j = TermV.jStart; j <= TermV.jEnd; j++)
		{
			aveU = (U(i, j - 1) + U(i + 1, j - 1) + U(i, j) + U(i + 1, j)) / 4.;
			if (aveU>0 && i > TermV.iStart)
			{
				Vx = dVdxM(i, j);
			}
			else
			{
				Vx = dVdxP(i, j);
			}

			if (V(i,j)>0)
			{
				Vy = dVdyM(i, j);
			}
			else
			{
				Vy = dVdyP(i, j);
			}
			TermV(i, j) = aveU*Vx + V(i, j)*Vy;
			//if ((i==TermV.iStart || i == TermV.iStart+1) && j==TermV.jEnd)
			//{
			//	cout << endl;
			//	cout << aveU << " " << Vx << endl;
			//	cout << V(i, j) << " " << Vy << endl;
			//	cout << TermV(i, j) << endl;
			//	cout << endl;
			//}
		}
	}
	//dVdxM.Variable("dvdxM");
	//dVdxP.Variable("dvdxP");
}

inline void EulerianFluidSolver2D::DiffusionTerm(const FD& U, const FD& V, FD& TermU, FD& TermV)
{
#pragma omp parallel for
	for (int i = TermU.iStart; i <= TermU.iEnd; i++)
	{
		for (int j = TermU.jStart; j <= TermU.jEnd; j++)
		{
			TermU(i, j) = U.dxxPhi(i, j) + U.dyyPhi(i, j);
		}
	}

#pragma omp parallel for 
	for (int i = TermV.iStart; i <= TermV.iEnd; i++)
	{
		for (int j = TermV.jStart; j <= TermV.jEnd; j++)
		{
			TermV(i, j) = V.dxxPhi(i, j) + V.dyyPhi(i, j);
		}
	}
}

inline double EulerianFluidSolver2D::AdaptiveTimeStep(const FD& velocity1, const FD& velocity2)
{
	double maxVel = 0;
	for (int i = velocity1.iStart; i <= velocity1.iEnd; i++)
	{
		for (int j = velocity1.jStart; j <= velocity1.jEnd; j++)
		{
			if (abs(velocity1(i,j))>maxVel)
			{
				maxVel = abs(velocity1(i, j));
			}
		}
	}
	for (int i = velocity2.iStart; i <= velocity2.iEnd; i++)
	{
		for (int j = velocity2.jStart; j <= velocity2.jEnd; j++)
		{
			if (abs(velocity2(i, j))>maxVel)
			{
				maxVel = abs(velocity2(i, j));
			}
		}
	}
	return cflCondition*min(gridU.dx, gridV.dy) / maxVel;
}

//////////////////////////////////////////////////////////////////////////
//
// "Accurate Projection Methods for the Incompressible Navier-Stokes Equations "
//                 Brown, Cortez, Minion (2001)
//     Sec 5.2 Projection Methods with a Lagged Pressure Term
//
/////////////////////////////////////////////////////////////////////////
inline void EulerianFluidSolver2D::NSSolver2ndOrder(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	const char* cmd;

	InitialCondition(example);
	gridP.Variable("Xp", "Yp");
	gridU.Variable("Xu", "Yu");
	gridV.Variable("Xv", "Yv");
	U.Variable("U");
	V.Variable("V");

	GenerateLinearSystemUV(UCNMatrix, gridUinner, 1);
	GenerateLinearSystemUV(VCNMatrix, gridVinner, 1);
	GenerateLinearSystemPhi(PhiCNMatrix, -gridP.dx2 / dt);
	//// CG solver 1
	UCN_CSR = CSR<double>(UCNMatrix);
	VCN_CSR = CSR<double>(VCNMatrix);
	PhiCN_CSR = CSR<double>(PhiCNMatrix);
	//// CG solver 2
	CGSolver::SparseA(UCNMatrix, Ua, Urow, Ucol, UnonzeroNum);
	CGSolver::SparseA(VCNMatrix, Va, Vrow, Vcol, VnonzeroNum);
	CGSolver::SparseA(PhiCNMatrix, Phia, Phirow, Phicol, PhinonzeroNum);

	str = string("quiver(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),axis([Xp(1)-(Xp(end)-Xp(1))/10 Xp(end)+(Xp(end)-Xp(1))/10 Yp(1)-(Yp(end)-Yp(1))/10 Yp(end)+(Yp(end)-Yp(1))/10]);");
	str = str + string("hold on,streamline(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,-100:0.1:100,-100:0.1:100),hold off;");
	MATLAB.Command(str.c_str());
	str = string("title(['iteration : ', num2str(") + to_string(0) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);
	
	oldU = U;
	oldV = V;
	double totalT = 0.0;
	for (iteration = 1; iteration <= maxIteration; iteration++)
	{
		cout << endl;
		cout << "********************************" << endl;
		cout << "       Iteration " << to_string(iteration) << " : Start" << endl;
		TVDRK3TimeAdvection();
		oldU.dataArray = originU.dataArray;
		oldV.dataArray = originV.dataArray;
		totalT += dt;
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;
		cout << "********************************" << endl;
		//P.Variable("P");
		U.Variable("U");
		V.Variable("V");
		MATLAB.Command("figure(1),axis([Xp(1)-(Xp(end)-Xp(1))/10 Xp(end)+(Xp(end)-Xp(1))/10 Yp(1)-(Yp(end)-Yp(1))/10 Yp(end)+(Yp(end)-Yp(1))/10])");

		str = string("quiver(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),axis([Xp(1)-(Xp(end)-Xp(1))/10 Xp(end)+(Xp(end)-Xp(1))/10 Yp(1)-(Yp(end)-Yp(1))/10 Yp(end)+(Yp(end)-Yp(1))/10]);");
		str = str + string("hold on,streamline(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,-100:0.1:100,-100:0.1:100),hold off;");
		MATLAB.Command(str.c_str());
		str = string("title(['iteration : ', num2str(") + to_string(iteration) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
		MATLAB.Command("divU =U(:,2:end)-U(:,1:end-1),divV =V(2:end,:)-V(1:end-1,:);div=divU+divV;");
		//MATLAB.Command("figure(2),plot(V(51,2:end-1),'-o'),grid on,");
		if (writeFile && iteration%writeOutputIteration == 0)
		{
			fileName = "pressure" + to_string(iteration);
			P.WriteFile(fileName);
			fileName = "xVelocity" + to_string(iteration);
			U.WriteFile(fileName);
			fileName = "yVelocity" + to_string(iteration);
			V.WriteFile(fileName);
		}
	}
}

inline void EulerianFluidSolver2D::EulerMethod2ndOrder()
{
	U.SaveOld();
	V.SaveOld();

	////////////////////////////////////////////////
	////     Projection Method 1 : advection    ////
	////////////////////////////////////////////////
	//if (iteration>=2)
	//{
	//	EulerMethod2ndOrder1();
	//}
	//else
	{
		EulerMethod2ndOrder1stIteration1();
	}

	///////////////////////////////////////////////////////////////
	////     Projection Method 2 : Recover U from Projection   ////
	///////////////////////////////////////////////////////////////
	EulerMethod2ndOrder2();

	////////////////////////////////////////////////////////////
	////     Projection Method 3 : New Gradient Pressure    ////
	////////////////////////////////////////////////////////////
	EulerMethod2ndOrder3();

	oldU.dataArray = U.dataArrayOld;
	oldV.dataArray = V.dataArrayOld;

	//MATLAB.Command("VVB = reshape(Vb, 49, 50)';");
	//MATLAB.Command("figure(2),subplot(1,2,1),plot(Vnew(51,2:end-1),'-o'),grid on, subplot(1,2,2),plot(advectionV(end,:),'-o'),grid on");
	//MATLAB.Command("figure(2),subplot(2,2,1),surf(Xu(2:end-1,2:end-1),Yu(2:end-1,2:end-1),Unew(2:end-1,2:end-1))");
	//MATLAB.Command("subplot(2,2,2),surf(Xv(2:end-1,2:end-1),Yv(2:end-1,2:end-1),Vnew(2:end-1,2:end-1))");
	//MATLAB.Command("subplot(2,2,3),surf(Xp(2:end-1,2:end-1),Yp(2:end-1,2:end-1),Phi(2:end-1,2:end-1))");
	//MATLAB.Command("subplot(2,2,4),surf(VVB)");
}

inline void EulerianFluidSolver2D::EulerMethod2ndOrder1()
{
	AdvectionTerm(U, V, advectionU1, advectionV1);
	AdvectionTerm(oldU, oldV, advectionU2, advectionV2);
	DiffusionTerm(U, V, diffusionU, diffusionV);

	double viscosity = 1;
	// 2nd-order Adams-Bashforth formula
#pragma omp parallel for
	for (int i = advectionU.iStart; i <= advectionU.iEnd; i++)
	{
		for (int j = advectionU.jStart; j <= advectionU.jEnd; j++)
		{
			advectionU(i, j) = 1. / 2.*(3 * advectionU1(i, j) - advectionU2(i, j));
		}
	}
#pragma omp parallel for
	for (int i = advectionV.iStart; i <= advectionV.iEnd; i++)
	{
		for (int j = advectionV.jStart; j <= advectionV.jEnd; j++)
		{
			advectionV(i, j) = 1. / 2.*(3 * advectionV1(i, j) - advectionV2(i, j));
		}
	}
	//advectionU.Variable("advectionU");
	//diffusionU.Variable("diffusionU");
	//advectionU1.Variable("advectionU1");
	//advectionU2.Variable("advectionU2");

	//advectionV.Variable("advectionV");
	//diffusionV.Variable("diffusionV");
	//advectionV1.Variable("advectionV1");
	//advectionV2.Variable("advectionV2");

	// Crank-Nicolson
	GenerateLinearSystemUV(Ub, U, gradientPx, advectionU, gridUinner, 1);
	GenerateLinearSystemUV(Vb, V, gradientPy, advectionV, gridVinner, 1);

	//Ub.Variable("Ub");
	//Vb.Variable("Vb");
	
	if (CGsolverNum == 1)
	{
		tempU = CGSolver::SolverCSR(UCN_CSR, Ub, gridU.dx*gridU.dy);
		tempV = CGSolver::SolverCSR(VCN_CSR, Vb, gridV.dx*gridV.dy);
	}
	else if (CGsolverNum == 2)
	{
		CGSolver::SolverSparse(UCNMatrix.iRes, Ua, Urow, Ucol, Ub, tempU);
		CGSolver::SolverSparse(VCNMatrix.iRes, Va, Vrow, Vcol, Vb, tempV);
	}
	//tempU.Variable("tempU");
	//tempV.Variable("tempV");

	int index;
#pragma omp parallel for private(index)
	for (int i = gridUinner.iStart; i <= gridUinner.iEnd; i++)
	{
		for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
		{
			index = (i - gridUinner.iStart) + (j - gridUinner.jStart)*gridUinner.iRes;
			U(i, j) = tempU(index);
		}
	}
#pragma omp parallel for private(index)
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		for (int j = gridVinner.jStart; j <= gridVinner.jEnd; j++)
		{
			index = (i - gridVinner.iStart) + (j - gridVinner.jStart)*gridVinner.iRes;
			V(i, j) = tempV(index);
		}
	}

	// Boundary : Linear extension. (But, Dirichlet로 줄 방법은 없나???)
#pragma omp parallel for
	for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
	{
		U(gridU.iStart, j) = -U(gridU.iStart + 1, j);
		U(gridU.iEnd, j) = -U(gridU.iEnd - 1, j);
	}
#pragma omp parallel for
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		V(i, gridV.jStart) = -V(i, gridV.jStart + 1);
		V(i, gridV.jEnd) = -V(i, gridV.jEnd - 1);
	}

	//U.Variable("Ustar");
	//V.Variable("Vstar");
	//MATLAB.Command("divUstar =Ustar(:,2:end)-Ustar(:,1:end-1),divVstar =Vstar(2:end,:)-Vstar(1:end-1,:);divstar=divUstar+divVstar;");
	//MATLAB.Command("quiver(Xp,Yp,Ustar(:,1:end-1),Vstar(1:end-1,:))");
}

inline void EulerianFluidSolver2D::EulerMethod2ndOrder2()
{
	GenerateLinearSystemPhi(Phib, -gridP.dx2 / dt);
	//Phib.Variable("phiB");
	
	if (CGsolverNum == 1)
	{
		tempPhi = CGSolver::SolverCSR(PhiCN_CSR, Phib, gridP.dx*gridP.dy);
	}
	else if (CGsolverNum == 2)
	{
		CGSolver::SolverSparse(PhiCNMatrix.iRes, Phia, Phirow, Phicol, Phib, tempPhi);
	}
	//tempPhi.Variable("tempPhi");

	int index;
#pragma omp parallel for private(index)
	for (int i = gridPinner.iStart; i <= gridPinner.iEnd; i++)
	{
		for (int j = gridPinner.jStart; j <= gridPinner.jEnd; j++)
		{
			index = (i - gridPinner.iStart) + (j - gridPinner.jStart)*gridPinner.iRes;
			Phi(i, j) = tempPhi(index);
		}
	}
#pragma omp parallel for
	for (int i = gridP.iStart; i <= gridP.iEnd; i++)
	{
		Phi(i, Phi.jStart) = Phi(i, Phi.jStart + 1);
		Phi(i, Phi.jEnd) = Phi(i, Phi.jEnd - 1);
	}
#pragma omp parallel for
	for (int j = gridP.jStart; j <= gridP.jEnd; j++)
	{
		Phi(Phi.iStart, j) = Phi(Phi.iStart + 1, j);
		Phi(Phi.iEnd, j) = Phi(Phi.iEnd - 1, j);
	}
	Phi(gridP.iStart, gridP.jEnd - 1) = 1;
	//Phi.Variable("Phi");

#pragma omp parallel for private(index)
	for (int i = gridUinner.iStart; i <= gridUinner.iEnd; i++)
	{
		for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
		{
			U(i, j) = U(i, j) - dt*(Phi(i, j) - Phi(i - 1, j))*Phi.oneOverdx;
		}
	}
#pragma omp parallel for private(index)
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		for (int j = gridVinner.jStart; j <= gridVinner.jEnd; j++)
		{
			V(i, j) = V(i, j) - dt*(Phi(i, j) - Phi(i, j - 1))*Phi.oneOverdy;
		}
	}
#pragma omp parallel for
	for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
	{
		U(gridU.iStart, j) = -U(gridU.iStart + 1, j);
		U(gridU.iEnd, j) = -U(gridU.iEnd - 1, j);
	}
#pragma omp parallel for
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		V(i, gridV.jStart) = -V(i, gridV.jStart + 1);
		V(i, gridV.jEnd) = -V(i, gridV.jEnd - 1);
	}

	//U.Variable("Unew");
	//V.Variable("Vnew");
	//MATLAB.Command("divUnew =Unew(:,2:end)-Unew(:,1:end-1),divVnew =Vnew(2:end,:)-Vnew(1:end-1,:);divnew=divUnew+divVnew;");
	//MATLAB.Command("quiver(Xp,Yp,Unew(:,1:end-1),Unew(1:end-1,:))");


}

inline void EulerianFluidSolver2D::EulerMethod2ndOrder3()
{
	double viscosity = 1;
#pragma omp parallel for
	for (int i = Phixxyy.grid.iStart; i <= Phixxyy.grid.iEnd; i++)
	{
		for (int j = Phixxyy.grid.jStart; j <= Phixxyy.grid.jEnd; j++)
		{
			Phixxyy(i, j) = Phi.dxxPhi(i, j) + Phi.dyyPhi(i, j);
		}
	}
	//Phixxyy.Variable("Phixxyy");

#pragma omp parallel for
	for (int i = gradientPx.grid.iStart; i <= gradientPx.grid.iEnd; i++)
	{
		for (int j = gradientPx.grid.jStart; j <= gradientPx.grid.jEnd; j++)
		{
			gradientPx(i, j) = gradientPx(i, j) + (Phi(i, j) - Phi(i - 1, j))*Phi.oneOverdx
				- viscosity / reynoldNum*dt / 2.0 * (Phixxyy(i, j) - Phixxyy(i - 1, j))*Phixxyy.oneOverdx;
		}
	}
#pragma omp parallel for
	for (int i = gradientPy.grid.iStart; i <= gradientPy.grid.iEnd; i++)
	{
		for (int j = gradientPy.grid.jStart; j <= gradientPy.grid.jEnd; j++)
		{
			gradientPy(i, j) = gradientPy(i, j) + (Phi(i, j) - Phi(i, j - 1))*Phi.oneOverdy
				- viscosity / reynoldNum*dt / 2.0 * (Phixxyy(i, j) - Phixxyy(i, j - 1))*Phixxyy.oneOverdy;
		}
	}

	//gradientPx.Variable("gradientPx");
	//gradientPy.Variable("gradientPy");
}

inline void EulerianFluidSolver2D::EulerMethod2ndOrder1stIteration1()
{
	AdvectionTerm(U, V, advectionU, advectionV);
	DiffusionTerm(U, V, diffusionU, diffusionV);

	double viscosity = 1;
	// 2nd-order Adams-Bashforth formula
	//advectionU.Variable("advectionU");
	//diffusionU.Variable("diffusionU");

	//advectionV.Variable("advectionV");
	//diffusionV.Variable("diffusionV");

	// Crank-Nicolson
	GenerateLinearSystemUV(Ub, U, gradientPx, advectionU, gridUinner, 1);
	GenerateLinearSystemUV(Vb, V, gradientPy, advectionV, gridVinner, 1);

	//Ub.Variable("Ub");
	//Vb.Variable("Vb");

	if (CGsolverNum == 1)
	{
		tempU = CGSolver::SolverCSR(UCN_CSR, Ub, gridU.dx*gridU.dy);
		tempV = CGSolver::SolverCSR(VCN_CSR, Vb, gridV.dx*gridV.dy);
	}
	else if (CGsolverNum == 2)
	{
		CGSolver::SolverSparse(UCNMatrix.iRes, Ua, Urow, Ucol, Ub, tempU);
		CGSolver::SolverSparse(VCNMatrix.iRes, Va, Vrow, Vcol, Vb, tempV);
	}
	//tempU.Variable("tempU");
	//tempV.Variable("tempV");

	int index;
#pragma omp parallel for private(index)
	for (int i = gridUinner.iStart; i <= gridUinner.iEnd; i++)
	{
		for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
		{
			index = (i - gridUinner.iStart) + (j - gridUinner.jStart)*gridUinner.iRes;
			U(i, j) = tempU(index);
		}
	}
#pragma omp parallel for private(index)
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		for (int j = gridVinner.jStart; j <= gridVinner.jEnd; j++)
		{
			index = (i - gridVinner.iStart) + (j - gridVinner.jStart)*gridVinner.iRes;
			V(i, j) = tempV(index);
		}
	}

	// Boundary : Linear extension. (But, Dirichlet로 줄 방법은 없나???)
#pragma omp parallel for
	for (int j = gridUinner.jStart; j <= gridUinner.jEnd; j++)
	{
		U(gridU.iStart, j) = -U(gridU.iStart + 1, j);
		U(gridU.iEnd, j) = -U(gridU.iEnd - 1, j);
	}
#pragma omp parallel for
	for (int i = gridVinner.iStart; i <= gridVinner.iEnd; i++)
	{
		V(i, gridV.jStart) = -V(i, gridV.jStart + 1);
		V(i, gridV.jEnd) = -V(i, gridV.jEnd - 1);
	}

	//U.Variable("Ustar");
	//V.Variable("Vstar");
	//MATLAB.Command("divUstar =Ustar(:,2:end)-Ustar(:,1:end-1),divVstar =Vstar(2:end,:)-Vstar(1:end-1,:);divstar=divUstar+divVstar;");
	//MATLAB.Command("quiver(Xp,Yp,Ustar(:,1:end-1),Vstar(1:end-1,:))");
}



inline void EulerianFluidSolver2D::GenerateLinearSystemUV(Array2D<double>& matrixA, const Grid2D& ipGrid, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int innerIStart = ipGrid.iStart;
	int innerIEnd = ipGrid.iEnd;
	int innerJStart = ipGrid.jStart;
	int innerJEnd = ipGrid.jEnd;
	int innerIRes = ipGrid.iRes;
	int innerJRes = ipGrid.jRes;

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
					matrixA(index) = scaling*(1 + dt / reynoldNum * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(rightIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdx2;
					matrixA(topIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + dt / reynoldNum * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -1. / 2. * dt / reynoldNum *ipGrid.oneOverdx2;
					matrixA(topIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + dt / reynoldNum * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -1. / 2. * dt / reynoldNum *ipGrid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1. / 2. * dt / reynoldNum* ipGrid.oneOverdx2;
					matrixA(topIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
				}

			}
			else if (j>innerJStart && j<innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + dt / reynoldNum * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(rightIndex) = scaling * -1. / 2. * dt / reynoldNum* ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
					matrixA(topIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;

				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + dt / reynoldNum * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -1. / 2. * dt / reynoldNum *ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
					matrixA(topIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;

				}
				else
				{
					matrixA(index) = scaling*(1 + dt / reynoldNum * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -1. / 2. * dt / reynoldNum *ipGrid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1. / 2. * dt / reynoldNum* ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
					matrixA(topIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
				}

			}
			else if (j == innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + dt / reynoldNum * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(rightIndex) = scaling * -1. / 2. * dt / reynoldNum* ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + dt / reynoldNum * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -1. / 2. * dt / reynoldNum *ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + dt / reynoldNum * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -1. / 2. * dt / reynoldNum *ipGrid.oneOverdx2;
					matrixA(rightIndex) = scaling * -1. / 2. * dt / reynoldNum* ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -1. / 2. * dt / reynoldNum * ipGrid.oneOverdy2;
				}

			}
		}
	}
}

inline void EulerianFluidSolver2D::GenerateLinearSystemUV(VectorND<double>& vectorB, const FD & vel, const FD & gradP, const FD & advec, const Grid2D& ipGrid, const double & scaling)
{
	int innerIStart = ipGrid.iStart;
	int innerIEnd = ipGrid.iEnd;
	int innerJStart = ipGrid.jStart;
	int innerJEnd = ipGrid.jEnd;
	int innerIRes = ipGrid.iRes;
	int innerJRes = ipGrid.jRes;

	int index;
#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;

			vectorB(index) = dt*(-gradP(i, j) - advec(i, j) + 1. / 2./reynoldNum*(vel.dxxPhi(i, j) + vel.dyyPhi(i, j)));

			if (i == innerIStart)
			{
				vectorB(index) += (2 * vel(i - 1, j) - vel.dataArrayOld(i - 1, i))*ipGrid.oneOverdx2* dt / reynoldNum / 2.;
			}
			if (i == innerIEnd)
			{
				vectorB(index) += (2 * vel(i + 1, j) - vel.dataArrayOld(i + 1, j))*ipGrid.oneOverdx2* dt / reynoldNum / 2.;
			}
			if (j == innerJStart)
			{
				vectorB(index) += (2 * vel(i, j - 1) - vel.dataArrayOld(i, j - 1))*ipGrid.oneOverdy2* dt / reynoldNum / 2.;
			}
			if (j == innerJEnd)
			{
				vectorB(index) += (2 * vel(i, j + 1) - vel.dataArrayOld(i, j + 1))*ipGrid.oneOverdy2* dt / reynoldNum / 2.;
			}

			//cout << endl;
			//cout << "(i,j) = (" << i << "," << j << ")" << endl;
			//cout << "index = " << index << endl;
			//cout << "U " << U(i + 1, j) << " " << U(i, j) << endl;
			//cout << "V " << V(i, j + 1) << " " << V(i, j) << endl;
			//cout << "B" << vectorB(index) << endl;
			vectorB(index) *= scaling;
		}
	}
}

inline void EulerianFluidSolver2D::GenerateLinearSystemPhi(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int innerIStart = gridPinner.iStart;
	int innerIEnd = gridPinner.iEnd;
	int innerJStart = gridPinner.jStart;
	int innerJEnd = gridPinner.jEnd;
	int innerIRes = gridPinner.iRes;
	int innerJRes = gridPinner.jRes;

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
					matrixA(index) = scaling * dt * (-1 * gridP.oneOverdx2 + -1 * gridP.oneOverdy2);
					matrixA(rightIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(topIndex) = scaling * dt * gridP.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling * dt * (-1 * gridP.oneOverdx2 + -1 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(topIndex) = scaling * dt * gridP.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling * dt * (-2 * gridP.oneOverdx2 + -1 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(topIndex) = scaling * dt * gridP.oneOverdy2;
				}
			}
			else if (j>innerJStart && j<innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling * dt * (-1 * gridP.oneOverdx2 + -2 * gridP.oneOverdy2);
					matrixA(rightIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * dt * gridP.oneOverdy2;
					matrixA(topIndex) = scaling * dt * gridP.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling * dt * (-1 * gridP.oneOverdx2 + -2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * dt * gridP.oneOverdy2;
					matrixA(topIndex) = scaling * dt * gridP.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling * dt * (-2 * gridP.oneOverdx2 + -2 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * dt * gridP.oneOverdy2;
					matrixA(topIndex) = scaling * dt * gridP.oneOverdy2;
				}
			}
			else if (j == innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling * dt * (-2 * gridP.oneOverdx2 + -1 * gridP.oneOverdy2);
					matrixA(rightIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * dt * gridP.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling * dt * (-1 * gridP.oneOverdx2 + -1 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * dt * gridP.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling * dt * (-2 * gridP.oneOverdx2 + -1 * gridP.oneOverdy2);
					matrixA(leftIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(rightIndex) = scaling * dt * gridP.oneOverdx2;
					matrixA(bottomIndex) = scaling * dt * gridP.oneOverdy2;
				}
			}
		}
	}
}

inline void EulerianFluidSolver2D::GenerateLinearSystemPhi(VectorND<double>& vectorB, const double & scaling)
{
	int innerIStart = gridPinner.iStart;
	int innerIEnd = gridPinner.iEnd;
	int innerJStart = gridPinner.jStart;
	int innerJEnd = gridPinner.jEnd;
	int innerIRes = gridPinner.iRes;
	int innerJRes = gridPinner.jRes;

	int index;
	Phi(gridP.iStart, gridP.jEnd - 1) = 1;
#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;

			vectorB(index) = (U(i + 1, j) - U(i, j))*U.grid.oneOverdx + (V(i, j + 1) - V(i, j))*V.grid.oneOverdy;
			
			if (i == innerIStart && j == innerJEnd)
			{
				vectorB(index) += - dt * Phi(gridP.iStart, gridP.jEnd - 1) * Phi.oneOverdx2;
			}
			//cout << endl;
			//cout << "(i,j) = (" << i << "," << j << ")" << endl;
			//cout << "index = " << index << endl;
			//cout << "U " << U(i + 1, j) << " " << U(i, j) << endl;
			//cout << "V " << V(i, j + 1) << " " << V(i, j) << endl;
			//cout << "B" << vectorB(index) << endl;

			vectorB(index) *= scaling;
		}
	}
}

