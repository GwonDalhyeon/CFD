#pragma once
#include "FluidSolver2D.h"
#include "EulerianMovingInterface.h"

class CoalescingDrop
{
public:
	CoalescingDrop(FluidSolver2D & ipFluid, MovingInterface& ipInterfaceSurfactant);
	~CoalescingDrop();

	int ExamNum;
	int iteration;


	int Nf; // The Num of Interfaces per Front.
	double rho1, rho2, rhoF;
	double mu1, mu2;
	double densityRatio;// 1 : Bubble, 0.1 : Liquid drop
	double viscosityRatio; // 1 : Bubble, 0.1 : Liquid drop
	double lengthscale;
	double gamma0;
	double timescale;


	///////////////////////////////
	//// Incompressible Fluid  ////
	///////////////////////////////
	FluidSolver2D& Fluid;

	Grid2D& grid = Fluid.grid;
	Grid2D& gridU = Fluid.gridU;
	Grid2D& gridV = Fluid.gridV;
	FD& U = Fluid.U;
	FD& V = Fluid.V;

	bool singularForce = false;
	FD SurfaceForceX;
	FD SurfaceForceY;
	FV SurfGradSurfTension;

	int& ProjectionOrder = Fluid.ProjectionOrder;
	double& reynoldNum = Fluid.reynoldNum;
	double& cflCondition = Fluid.cflCondition;
	double& Pe = Fluid.Pe;
	double& Oh = Fluid.Oh;
	double& We = Fluid.We;
	double& dt = Fluid.dt;

	double& Ca = Fluid.Ca;
	double& Xi = Fluid.Xi;
	double& El = Fluid.El;

	int& maxIteration = Fluid.maxIteration;
	int& writeOutputIteration = Fluid.writeOutputIteration;

	///////////////////////////////
	//// Insoluble Surfactant  ////
	///////////////////////////////
	MovingInterface& InterfaceSurfactant;
	FD& Surfactant = InterfaceSurfactant.Surfactant;
	FD& SurfaceTension = InterfaceSurfactant.SurfaceTension;

	LS& levelSet = InterfaceSurfactant.levelSet;
	double& totalT = InterfaceSurfactant.totalT;


	inline void InitialCondition(const int& example);

	/////////////////////////////////////
	//// Coalescing Solver : Bubble  ////
	/////////////////////////////////////
	void CoalescingBubbleSolver(const int& example);

	// Navier-Stokes equation solver
	inline void NSSolver();
	inline void EulerMethod();
	inline void EulerMethod1();
	inline void EulerMethod2ndOrder();
	inline void EulerMethod2ndOrder1();
	inline void EulerMethod2ndOrder1stIteration1();

	inline void ComputeSurfaceForce();
	inline void GenerateLinearSystemUV(Array2D<double>& matrixA, const Grid2D& ipGrid, const double & scaling, CSR<double>& csrForm);
	inline void GenerateLinearSystemUV(VectorND<double>& vectorB, const FD & vel, const FD & gradP, const FD & advec, const FD& force, const Grid2D& ipGrid, const double & scaling);


	inline void PlotSurfactant();
	inline void PlotVelocity();
private:

};

CoalescingDrop::CoalescingDrop(FluidSolver2D & ipFluid, MovingInterface& ipInterfaceSurfactant)
	:Fluid(ipFluid), InterfaceSurfactant(ipInterfaceSurfactant)
{
}

CoalescingDrop::~CoalescingDrop()
{
}

inline void CoalescingDrop::InitialCondition(const int & example)
{
	ExamNum = example;
	if (example == 1)
	{
		cout << "*************************" << endl;
		cout << "        Simulations of surfactant effects  " << endl;
		cout << " on the dynamics of coalescing drops and bubbles" << endl;
		cout << "      --DW Martin and F Blanchette-- " << endl;
		cout << "               Example 1 " << endl;
		cout << "*************************" << endl;

		int gridSize = 100;
		grid = Grid2D(-2.5, 2.5, gridSize + 1, 0, 5, gridSize + 1);

		Nf = 2;
		rho2 = 1.25; rho1 = rho2;
		mu2 = 1.81*pow(10, -5); mu1 = mu2;
		densityRatio = rho2 / rho1; // 1 : Bubble, 0.1 : Liquid drop
		viscosityRatio = mu2 / mu1; // 1 : Bubble, 0.1 : Liquid drop
		lengthscale = 0.01;
		gamma0 = 3.0*pow(10, -3); //??????????
		timescale = sqrt(rho1*pow(lengthscale, 3) / (Nf*gamma0));

		//// Initialize Surfactant Fields
		levelSet = LS(grid, 3 * grid.dx);
		LS levelSet1(grid, 3 * grid.dx);
		LS levelSet2(grid, 3 * grid.dx);
		double radius = 0.5;
		VT center(0, 2.5 - grid.dy / 10);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet1(i, j) = sqrt((grid(i, j).x - center.x)*(grid(i, j).x - center.x) + (grid(i, j).y - center.y)*(grid(i, j).y - center.y)) - radius;
			}
		}
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet2(i, j) = grid(i, j).y - 2 - grid.dy / 10;
			}
		}

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (levelSet1(i, j) * levelSet2(i, j) <= 0)
				{
					levelSet(i, j) = min(levelSet1(i, j), levelSet2(i, j));
				}
				else if (abs(levelSet1(i, j)) < abs(levelSet2(i, j)))
				{
					levelSet(i, j) = levelSet1(i, j);
				}
				else
				{
					levelSet(i, j) = levelSet2(i, j);
				}
			}
		}

		InterfaceSurfactant.InitialCondition(8);


		//// Initialize Velocity Fields
		Fluid.InitialCondition(5);
		
		We = rho2 * lengthscale * lengthscale / gamma0;
		Oh = mu1 / sqrt(lengthscale*rho1*Nf*gamma0);

#pragma omp parallel for
		for (int i = gridV.iStart; i <= gridV.iEnd; i++)
		{
			for (int j = gridV.jStart; j <= gridV.jEnd - 1; j++)
			{
				if (levelSet1(i, j) <= 0)
				{
					V(i, j) = -We;
					Fluid.originV(i, j) = V(i, j);
				}
			}
		}

		cflCondition = 1.0 / 8.0;
		//dt = cflCondition*min(grid.dx, grid.dy);
		dt = 0.001;

		singularForce = false;
		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		double finalTime = 10.0;
		maxIteration = ceil(finalTime / dt); // Up to 2 sec.
		totalT = 0;
	}
}

inline void CoalescingDrop::CoalescingBubbleSolver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	const char* cmd;


	InitialCondition(example);

	grid.Variable();
	levelSet.phi.Variable("phi0");
	levelSet.tube.Variable("Tube");
	Surfactant.Variable("Surfactant0");

	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	//MATLAB.Command("subplot(2,1,1),surf(X,Y,Surfactant),subplot(2,1,2),contour(X,Y,phi0,[0 0],'b'),grid on,axis equal");
	MATLAB.Command("IntSur0 = sum(sum(Surfactant0.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");

	PlotVelocity();
	MATLAB.WriteImage("surfactant", 0, "fig");
	MATLAB.WriteImage("surfactant", 0, "png");

	int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy)) * 2;
	int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.2*min(grid.dx, grid.dy)));
	for (iteration = 1; iteration <= maxIteration; iteration++)
	{
		cout << "*******************************************************************" << endl;
		cout << "       Iteration " << to_string(iteration) << " : Start" << endl;
		totalT += dt;
		//// Step 1-1 : Surfactant Diffusion
		cout << "Diffusion Start" << endl;
		InterfaceSurfactant.LSurfactantDiffusion(iteration);
		cout << "Diffusion End" << endl;

		//// Step 1-2 : New Surface Tension
		//InterfaceSurfactant.DimlessNonlinearLangmu1rEOS(2);
		//InterfaceSurfactant.SurfaceTension.Variable("SurfaceTension");

		//// Step 2 : Navier-Stokes equation
		NSSolver();



		AdvectionMethod2D<double>::LLSPropagatingTVDRK3(levelSet, U, V, dt);
		AdvectionMethod2D<double>::LLSReinitializationTVDRK3(levelSet, dt, reinitialIter);

		//InterfaceSurfactant.ConserveSurfactantFactorBeta();

		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);
		levelSet.UpdateInterface();
		levelSet.UpdateLLS();


		//MATLAB.Command("subplot(2,1,1)");
		//PlotSurfactant();
		//MATLAB.Command("subplot(2,1,2)");
		PlotVelocity();
		if (iteration % 10000000 == 0)
		{
			MATLAB.WriteImage("surfactant", iteration, "fig");
			MATLAB.WriteImage("surfactant", iteration, "png");
		}
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;
		cout << "*******************************************************************" << endl;
	}
}

inline void CoalescingDrop::NSSolver()
{
	Fluid.originU.dataArray = U.dataArray;
	Fluid.originV.dataArray = V.dataArray;

	/////////////////
	//// Step 1  ////
	/////////////////
	if (ProjectionOrder == 1)
	{
		EulerMethod();
	}
	else if (ProjectionOrder == 2)
	{
		EulerMethod2ndOrder();
	}
	//U.Variable("U1");
	//V.Variable("V1");
	//MATLAB.Command("quiver(Xp,Yp,U1(:,1:end-1),V1(1:end-1,:))");

	/////////////////
	//// Step 2  ////
	/////////////////
	if (ProjectionOrder == 1)
	{
		EulerMethod();
	}
	else if (ProjectionOrder == 2)
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
			U(i, j) = 3. / 4. * Fluid.originU(i, j) + 1. / 4. * U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = gridV.iStart; i <= gridV.iEnd; i++)
	{
		for (int j = gridV.jStart; j <= gridV.jEnd; j++)
		{
			V(i, j) = 3. / 4. * Fluid.originV(i, j) + 1. / 4. * V(i, j);
		}
	}
	//U.Variable("U2");
	//V.Variable("V2");
	//MATLAB.Command("quiver(Xp,Yp,U22(:,1:end-1),V22(1:end-1,:))");

	/////////////////
	//// Step 3  ////
	/////////////////
	if (ProjectionOrder == 1)
	{
		EulerMethod();
	}
	else if (ProjectionOrder == 2)
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
			U(i, j) = 1. / 3. * Fluid.originU(i, j) + 2. / 3. * U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = gridV.iStart; i <= gridV.iEnd; i++)
	{
		for (int j = gridV.jStart; j <= gridV.jEnd; j++)
		{
			V(i, j) = 1. / 3. * Fluid.originV(i, j) + 2. / 3. * V(i, j);
		}
	}
	//U.Variable("U3");
	//V.Variable("V3");
	//MATLAB.Command("quiver(Xp,Yp,U32(:,1:end-1),V32(1:end-1,:))");
}

inline void CoalescingDrop::EulerMethod()
{
	////////////////////////////////////////////////
	////     Projection Method 1 : advection    ////
	////////////////////////////////////////////////

	EulerMethod1();

	////////////////////////////////////////////////
	////     Projection Method 2 : Poisson Eq   ////
	////////////////////////////////////////////////
	Fluid.EulerMethodStep2();

	//////////////////////////////////////////////
	////     Projection Method 3 : New U,V    ////
	//////////////////////////////////////////////
	Fluid.EulerMethodStep3();
}

inline void CoalescingDrop::EulerMethod1()
{
	//// Compute Surface Force
	if (singularForce)
	{
		ComputeSurfaceForce();
	}

	Array2D<double>& K1U = U.K1;
	Array2D<double>& K1V = V.K1;
	Fluid.AdvectionTerm(U, V, Fluid.AdvectionU, Fluid.AdvectionV);
	Fluid.DiffusionTerm(U, V, Fluid.DiffusionU, Fluid.DiffusionV);
	//AdvectionU.Variable("AdvectionU");
	//DiffusionU.Variable("DiffusionU");
	//AdvectionV.Variable("AdvectionV");
	//DiffusionV.Variable("DiffusionV");

#pragma omp parallel for
	for (int i = K1U.iStart; i <= K1U.iEnd; i++)
	{
		for (int j = K1U.jStart; j <= K1U.jEnd; j++)
		{
			K1U(i, j) = dt*(-Fluid.AdvectionU(i, j) + 1. / reynoldNum*Fluid.DiffusionU(i, j) + SurfaceForceX(i, j));
			U(i, j) = U(i, j) + K1U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = K1V.iStart; i <= K1V.iEnd; i++)
	{
		for (int j = K1V.jStart; j <= K1V.jEnd; j++)
		{
			K1V(i, j) = dt*(-Fluid.AdvectionV(i, j) + 1. / reynoldNum*Fluid.DiffusionV(i, j) + SurfaceForceY(i, j));
			V(i, j) = V(i, j) + K1V(i, j);
		}
	}
	//K1U.Variable("k1u");
	//K1V.Variable("k1v");
	Fluid.TreatVelocityBC(U, V);


	//U.Variable("Ustar");
	//V.Variable("Vstar");
	//MATLAB.Command("divUstar =Ustar(:,2:end)-Ustar(:,1:end-1),divVstar =Vstar(2:end,:)-Vstar(1:end-1,:);divstar=divUstar+divVstar;");
	//MATLAB.Command("quiver(Xp,Yp,Ustar(:,1:end-1),Vstar(1:end-1,:)");
}

inline void CoalescingDrop::EulerMethod2ndOrder()
{
	U.SaveOld();
	V.SaveOld();

	////////////////////////////////////////////////
	////     Projection Method 1 : advection    ////
	////////////////////////////////////////////////

	//if (iteration >= 2)
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
	Fluid.EulerMethod2ndOrder2();

	////////////////////////////////////////////////////////////
	////     Projection Method 3 : New Gradient Pressure    ////
	////////////////////////////////////////////////////////////
	Fluid.EulerMethod2ndOrder3();

	Fluid.oldU.dataArray = U.dataArrayOld;
	Fluid.oldV.dataArray = V.dataArrayOld;

	//MATLAB.Command("VVB = reshape(Vb, 49, 50)';");
	//AdvectionV.Variable("AdvectionV");
	//MATLAB.Command("figure(2),subplot(1,2,1),plot(Vnew(51,2:end-1),'-o'),grid on, subplot(1,2,2),plot(AdvectionV(end,:),'-o'),grid on");
	//MATLAB.Command("figure(2),subplot(2,2,1),surf(Xu(2:end-1,2:end-1),Yu(2:end-1,2:end-1),Unew(2:end-1,2:end-1))");
	//MATLAB.Command("subplot(2,2,2),surf(Xv(2:end-1,2:end-1),Yv(2:end-1,2:end-1),Vnew(2:end-1,2:end-1))");
	//MATLAB.Command("subplot(2,2,3),surf(Xp(2:end-1,2:end-1),Yp(2:end-1,2:end-1),Phi(2:end-1,2:end-1))");
	//MATLAB.Command("subplot(2,2,4),surf(VVB)");
}

inline void CoalescingDrop::EulerMethod2ndOrder1()
{
}

inline void CoalescingDrop::EulerMethod2ndOrder1stIteration1()
{
	if (singularForce)
	{
		ComputeSurfaceForce();
		//SurfaceForceX.Variable("SurfaceForceX");
		//SurfaceForceY.Variable("SurfaceForceY");
	}


	Fluid.AdvectionTerm(U, V, Fluid.AdvectionU, Fluid.AdvectionV);
	Fluid.DiffusionTerm(U, V, Fluid.DiffusionU, Fluid.DiffusionV);

	double viscosity = 1;
	//// 2nd-order Adams-Bashforth formula
	//Fluid.AdvectionU.Variable("AdvectionU");
	//Fluid.DiffusionU.Variable("DiffusionU");

	//Fluid.AdvectionV.Variable("AdvectionV");
	//Fluid.DiffusionV.Variable("DiffusionV");

	//// Crank-Nicolson
	GenerateLinearSystemUV(Fluid.Ub, U, Fluid.gradientPx, Fluid.AdvectionU, SurfaceForceX, Fluid.U.innerGrid, 1);
	GenerateLinearSystemUV(Fluid.Vb, V, Fluid.gradientPy, Fluid.AdvectionV, SurfaceForceY, Fluid.V.innerGrid, 1);

	//Fluid.Ub.Variable("Ub");
	//Fluid.Vb.Variable("Vb");


	CGSolver::Solver(Fluid.UCN_CSR, Fluid.Ub, Fluid.tempU);
	CGSolver::Solver(Fluid.VCN_CSR, Fluid.Vb, Fluid.tempV);


	//tempU.Variable("tempU");
	//tempV.Variable("tempV");

	int index;
#pragma omp parallel for private(index)
	for (int i = Fluid.U.innerGrid.iStart; i <= Fluid.U.innerGrid.iEnd; i++)
	{
		for (int j = Fluid.U.innerGrid.jStart; j <= Fluid.U.innerGrid.jEnd; j++)
		{
			index = (i - Fluid.U.innerGrid.iStart) + (j - Fluid.U.innerGrid.jStart)*Fluid.U.innerGrid.iRes;
			U(i, j) = Fluid.tempU(index);
		}
	}
#pragma omp parallel for private(index)
	for (int i = Fluid.V.innerGrid.iStart; i <= Fluid.V.innerGrid.iEnd; i++)
	{
		for (int j = Fluid.V.innerGrid.jStart; j <= Fluid.V.innerGrid.jEnd; j++)
		{
			index = (i - Fluid.V.innerGrid.iStart) + (j - Fluid.V.innerGrid.jStart)*Fluid.V.innerGrid.iRes;
			V(i, j) = Fluid.tempV(index);
		}
	}

	Fluid.TreatVelocityBC(U, V);

	//U.Variable("Ustar");
	//V.Variable("Vstar");
	//MATLAB.Command("divUstar =Ustar(:,2:end)-Ustar(:,1:end-1),divVstar =Vstar(2:end,:)-Vstar(1:end-1,:);divstar=divUstar+divVstar;");
	//MATLAB.Command("quiver(Xp,Yp,Ustar(:,1:end-1),Vstar(1:end-1,:))");
}

inline void CoalescingDrop::ComputeSurfaceForce()
{
	FD& meanCurvature = levelSet.meanCurvature;
	FV& unitNormal = levelSet.unitNormal;
	//Array2D<VT>& gradient = Surfactant.gradient;

	int ComputedTubeRange = 2;
	levelSet.LComputeMeanCurvature(ComputedTubeRange);
	levelSet.LComputeUnitNormal(ComputedTubeRange);

	int i, j;
	//// Surface Tension Surface Gradient
#pragma omp parallel for private(i, j)
	for (int k = 1; k <= levelSet.numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (levelSet.tube(i, j) <= ComputedTubeRange)
		{
			//gradient(i, j) = VT(Surfactant.dxPhi(i, j), Surfactant.dyPhi(i, j));
			//SurfGradSurfTension(i, j) = Surfactant.gradient(i, j)
			//	- DotProduct(unitNormal(i, j), gradient(i, j))*unitNormal(i, j);
			SurfaceTension(i, j) = gamma0;
			SurfaceForceX(i, j) = 2 * meanCurvature(i, j)*SurfaceTension(i, j)*unitNormal(i, j).x + SurfGradSurfTension(i, j).x;
			SurfaceForceX(i, j) *= AdvectionMethod2D<double>::DeltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude();
			SurfaceForceY(i, j) = 2 * meanCurvature(i, j)*SurfaceTension(i, j)*unitNormal(i, j).y + SurfGradSurfTension(i, j).y;
			SurfaceForceY(i, j) *= -AdvectionMethod2D<double>::DeltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude();

		}
		else
		{
			SurfGradSurfTension(i, j) = 0;

			SurfaceForceX(i, j) = 0;
			SurfaceForceY(i, j) = 0;
		}
	}
}

inline void CoalescingDrop::GenerateLinearSystemUV(Array2D<double>& matrixA, const Grid2D & ipGrid, const double & scaling, CSR<double>& csrForm)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	matrixA.initialize(1, ipGrid.iRes*ipGrid.jRes, 1, ipGrid.iRes*ipGrid.jRes);

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
					matrixA(index) = scaling*(1 + dt * Oh * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(rightIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdx2;
					matrixA(topIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + dt * Oh * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -0.5 * dt * Oh *ipGrid.oneOverdx2;
					matrixA(topIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + dt * Oh * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -0.5 * dt * Oh *ipGrid.oneOverdx2;
					matrixA(rightIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdx2;
					matrixA(topIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
				}

			}
			else if (j > innerJStart && j < innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + dt * Oh * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(rightIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
					matrixA(topIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;

				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + dt * Oh * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -0.5 * dt * Oh *ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
					matrixA(topIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;

				}
				else
				{
					matrixA(index) = scaling*(1 + dt * Oh * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -0.5 * dt * Oh *ipGrid.oneOverdx2;
					matrixA(rightIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
					matrixA(topIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
				}

			}
			else if (j == innerJEnd)
			{
				if (i == innerIStart)
				{
					matrixA(index) = scaling*(1 + dt * Oh * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(rightIndex) = scaling * -0.5 * dt * Oh* ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
				}
				else if (i == innerIEnd)
				{
					matrixA(index) = scaling*(1 + dt * Oh * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -0.5 * dt * Oh *ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
				}
				else
				{
					matrixA(index) = scaling*(1 + dt * Oh * (ipGrid.oneOverdx2 + ipGrid.oneOverdy2));
					matrixA(leftIndex) = scaling * -0.5 * dt * Oh *ipGrid.oneOverdx2;
					matrixA(rightIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdx2;
					matrixA(bottomIndex) = scaling * -0.5 * dt * Oh * ipGrid.oneOverdy2;
				}

			}
		}
	}
	csrForm = CSR<double>(matrixA);
	matrixA.Delete();
}

inline void CoalescingDrop::GenerateLinearSystemUV(VectorND<double>& vectorB, const FD & vel, const FD & gradP, const FD & advec, const FD & force, const Grid2D & ipGrid, const double & scaling)
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
			vectorB(index) = vel(i, j) + dt*(-gradP(i, j) - advec(i, j) + 0.5 * Oh * (vel.dxxPhi(i, j) + vel.dyyPhi(i, j)) + force(i, j));

			if (i == innerIStart)
			{
				vectorB(index) += -0.5 * dt * Oh * (2 * vel(i - 1, j) - vel.dataArrayOld(i - 1, j))*ipGrid.oneOverdx2;
			}
			if (i == innerIEnd)
			{
				vectorB(index) += -0.5 * dt * Oh * (2 * vel(i + 1, j) - vel.dataArrayOld(i + 1, j))*ipGrid.oneOverdx2;
			}
			if (j == innerJStart)
			{
				vectorB(index) += -0.5 * dt * Oh * (2 * vel(i, j - 1) - vel.dataArrayOld(i, j - 1))*ipGrid.oneOverdy2;
			}
			if (j == innerJEnd)
			{
				vectorB(index) += -0.5 * dt * Oh * (2 * vel(i, j + 1) - vel.dataArrayOld(i, j + 1))*ipGrid.oneOverdy2;
			}

			vectorB(index) *= scaling;
		}
	}
}

inline void CoalescingDrop::PlotSurfactant()
{
	string str;

	Surfactant.Variable("Surfactant");
	levelSet.tube.Variable("Tube");
	MATLAB.Command("SurTube1 = Surfactant.*(Tube<=1);");
	//MATLAB.Command("contour(X,Y,Tube,'r'),axis equal,axis([X(1) X(end) Y(1) Y(end)]), hold on,surf(X,Y,SurTube1), h=colorbar,h.Limits=[0 max(max(SurTube1))],hold off;set(gca,'fontsize',20)");
	MATLAB.Command("surf(X,Y,SurTube1), h=colorbar,h.Limits=[0 max(max(SurTube1))],axis equal,axis([X(1) X(end) Y(1) Y(end)]), hold on,hold off;set(gca,'fontsize',20)");

	//// Measure Surfactant Loss ////
	MATLAB.Command("IntSur = sum(sum(SurTube1.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");
	MATLAB.Command("loss = (IntSur-IntSur0)/IntSur0*100;");
	MATLAB.Variable("i", iteration);
	MATLAB.Variable("totalT", totalT);
	str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', error(%)  :',num2str(loss),'%']);");
	MATLAB.Command(str.c_str());
}

inline void CoalescingDrop::PlotVelocity()
{
	string str;

	U.Variable("U");
	V.Variable("V");
	levelSet.phi.Variable("phi");

	str = string("quiver(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),axis([X(1)-(X(end)-X(1))/10 X(end)+(X(end)-X(1))/10 Y(1)-(Y(end)-Y(1))/10 Y(end)+(Y(end)-Y(1))/10]);");
	//str = str + string("hold on,streamline(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,-100:0.1:100,-100:0.1:100),hold off;");
	str = str + string("hold on, contour(X,Y,phi,[0 0],'r'), hold off;axis([X(1)-(X(end)-X(1))/10 X(end)+(X(end)-X(1))/10 Y(1)-(Y(end)-Y(1))/10 Y(end)+(Y(end)-Y(1))/10]),axis equal");
	MATLAB.Command(str.c_str());
	str = string("title(['iteration : ', num2str(") + to_string(iteration) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	MATLAB.Command(str.c_str());
	MATLAB.Command("divU =U(:,2:end)-U(:,1:end-1),divV =V(2:end,:)-V(1:end-1,:);div=divU+divV;");
}





