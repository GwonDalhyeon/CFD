#pragma once

#include "AdvectionMethod2D.h"
#include "FluidSolver2D.h"
#include "EulerianMovingInterface.h"

class InsolubleSurfactant
{
public:
	InsolubleSurfactant(FluidSolver2D & ipFluid, MovingInterface& ipInterfaceSurfactant);
	~InsolubleSurfactant();

	int ExamNum;
	int iteration = 0;
	///////////////////////////////
	//// Incompressible Fluid  ////
	///////////////////////////////
	FluidSolver2D& Fluid;

	Grid2D& grid = Fluid.grid;
	Grid2D& gridU = Fluid.gridU;
	Grid2D& gridV = Fluid.gridV;
	FD& Pressure = Fluid.Pressure;
	FD& U = Fluid.U;
	FD& V = Fluid.V;

	bool singularForce = true;
	FD SurfaceForceX;
	FD SurfaceForceY;
	FV SurfGradSurfTension;


	int& ProjectionOrder = Fluid.ProjectionOrder;
	double& Re = Fluid.Re;
	double& cflCondition = Fluid.cflCondition;
	double& Ca = Fluid.Ca;
	double& Xi = Fluid.Xi;
	double& El = Fluid.El;
	double& Pe = Fluid.Pe;

	int& VspatialOrder = Fluid.spatialOrder; // Velocity spatial order
	int& temporalOrder = Fluid.temporalOrder;

	double& dt = Fluid.dt;

	int& maxIteration = Fluid.maxIteration;
	int& writeOutputIteration = Fluid.writeOutputIteration;


	///////////////////////////////
	//// Insoluble Surfactant  ////
	///////////////////////////////
	MovingInterface& InterfaceSurfactant;
	FD& Surfactant = InterfaceSurfactant.Surfactant;
	FD& SurfaceTension = InterfaceSurfactant.SurfaceTension;
	LS& levelSet = InterfaceSurfactant.levelSet;

	int& LspatialOrder = InterfaceSurfactant.LspatialOrder; // Level Set spatial order

	double totalT = 0;



	inline void InitialCondition(const int& example);

	/////////////////////////////////////////////////////
	//// Surfactant + Navier-Stokes equation solver  ////
	/////////////////////////////////////////////////////
	inline void ContinuumMethodWithSurfactantSolver(const int& example);

	inline void NSSolver();
	inline void EulerMethod();
	inline void EulerMethod1();
	inline void EulerMethod2ndOrder();
	inline void EulerMethod2ndOrder1();
	inline void EulerMethod2ndOrder1stIteration1();

	inline void ComputeSurfaceForce();
	inline void GenerateLinearSystemUV(VectorND<double>& vectorB, const FD & vel, const FD & gradP, const FD & advec, const FD& force, const Grid2D& ipGrid, const double & scaling);
	inline void PlotSurfactant();
	inline void PlotVelocity();


	/////////////////////////////////////////////////////
	////              CoalescingDrop                 ////
	/////////////////////////////////////////////////////
	int Nf = 1; // The Num of Interfaces per Front. a Liquid Drop : 1, a Soap Bubble : 2
	double densityI = Fluid.densityI, densityE = Fluid.densityE; // Air : 1.25 Kg/m^3
	double densityF = 1; // Film Density. 10^3 Kg/m^3, while the thickness has a typical value of h0 = 10E-6 m.
	double viscosityI = Fluid.viscosityI, viscosityE = Fluid.viscosityE; // Air : 1.81*E-5 Pa s
	double densityRatio = 1;// = densityE/densityI;  a Liquid Drop : 0.1, a Soap Bubble : 1
	double viscosityRatio = 1; // = muE/muI;  a Liquid Drop : 0.1, a Soap Bubble : 1
	
	double lengthscale = 1; // Standard Length Scale
	double timescale = 1; // 

	double& Oh = Fluid.Oh;
	double& We = Fluid.We;
	double& Bo = Fluid.Bo;
	const double& gravity = Fluid.gravity;

	double BoF = 0; // a Film Bond Number. a Liquid Drop : 0, a Soap Bubble : nonzero.

	double& gamma0 = Fluid.gamma0;

private:

};

inline InsolubleSurfactant::InsolubleSurfactant(FluidSolver2D & ipFluid, MovingInterface & ipInterfaceSurfactant)
	:Fluid(ipFluid), InterfaceSurfactant(ipInterfaceSurfactant)
{
}

InsolubleSurfactant::~InsolubleSurfactant()
{
}

inline void InsolubleSurfactant::InitialCondition(const int & example)
{
	ExamNum = example;
	if (example == 1)
	{
		cout << "*******************************************************************" << endl;
		cout << "			       A level-set continuum method  " << endl;
		cout << "		 for two-phase flows with insoluble surfactant" << endl;
		cout << "				--JJ Xu, Y Yang, J Lowengrub-- " << endl;
		cout << "						 Example 1 " << endl;
		cout << "*******************************************************************" << endl;

		int gridSize = 200;
		grid = Grid2D(-5, 5, gridSize + 1, -2, 2, gridSize * 2. / 5. + 1);

		// Initialize Velocity Fields
		Fluid.InitialCondition(3);
		ProjectionOrder = 1;
		Re = 100;
		Ca = 0.3;
		Xi = 0.3;
		El = 0.2;
		Pe = 10;
		cflCondition = 0.4;
		dt = cflCondition*min(grid.dx, grid.dy);

		//singularForce = false;
		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		// Initialize Surfactant Fields
		InterfaceSurfactant.InitialCondition(7);
		InterfaceSurfactant.cflCondition = cflCondition;
		InterfaceSurfactant.dt = dt;

		double finalT = 5;// Up to 5 sec.
		maxIteration = ceil(finalT / dt);
		writeOutputIteration = 30;
	}

	if (example == 2)
	{
		cout << "*******************************************************************" << endl;
		cout << "			       A level-set continuum method  " << endl;
		cout << "		 for two-phase flows with insoluble surfactant" << endl;
		cout << "				--JJ Xu, Y Yang, J Lowengrub-- " << endl;
		cout << "						 Example 4.3 " << endl;
		cout << "*******************************************************************" << endl;

		int gridSize = 100;
		double xLength = 2, yLength = 1.2;
		grid = Grid2D(-xLength, xLength, xLength * gridSize + 1, -yLength, yLength, yLength*gridSize + 1);

		// Initialize Velocity Fields
		//Fluid.isPCG = true;
		Fluid.InitialCondition(3);
		ProjectionOrder = 1;
		Re = 100;
		Ca = 0.1;
		Xi = 0.3;
		El = 0.2;
		Pe = 10;
		cflCondition = 0.4;
		dt = cflCondition*min(grid.dx, grid.dy);
		
		//singularForce = false;
		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		// Initialize Surfactant Fields
		InterfaceSurfactant.InitialCondition(8);
		InterfaceSurfactant.cflCondition = cflCondition;
		InterfaceSurfactant.dt = dt;

		double finalT = 10;// Up to 5 sec.
		maxIteration = ceil(finalT / dt);
		writeOutputIteration = 30;
	}

	if (example == 3)
	{
		cout << "*************************" << endl;
		cout << "        Simulations of surfactant effects  " << endl;
		cout << " on the dynamics of coalescing drops and bubbles" << endl;
		cout << "      --DW Martin and F Blanchette-- " << endl;
		cout << "               Example 1 " << endl;
		cout << "*************************" << endl;

		int gridSize = 100;
		double xLength = 2.5;
		grid = Grid2D(-xLength, xLength, gridSize + 1, 0, 2*xLength, gridSize + 1);

		Nf = 2; // a Liquid Drop : 1, a Soap Bubble : 2
		densityE = 1.25; //
		densityI = densityE;
		viscosityE = 1.81*pow(10, -5);
		viscosityI = viscosityE;
		densityRatio = densityE / densityI; // 1 : Bubble, 0.1 : Liquid drop
		viscosityRatio = viscosityE / viscosityI; // 1 : Bubble, 0.1 : Liquid drop
		lengthscale = 0.01;
		gamma0 = 3.0*pow(10, -3); //??????????
		timescale = sqrt(densityI*pow(lengthscale, 3) / (Nf*gamma0));

		//// Initialize Surfactant Fields
		levelSet = LS(grid, 3 * grid.dx);
		LS levelSet1(grid, 3 * grid.dx);
		LS levelSet2(grid, 3 * grid.dx);
		double radius = 0.5;
		VT center(0, 2.5);
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
				levelSet2(i, j) = grid(i, j).y - 2 + grid.dx;
			}
		}

		levelSet.CombineLevelSet(levelSet1, levelSet2);
		levelSet.InitialTube();

		InterfaceSurfactant.InitialCondition(9);

		//// Initialize Velocity Fields
		Fluid.InitialCondition(5);

		We = densityE * lengthscale * lengthscale / gamma0;
		Oh = viscosityI / sqrt(lengthscale*densityI*Nf*gamma0);

#pragma omp parallel for
		for (int i = gridV.iStart; i <= gridV.iEnd; i++)
		{
			for (int j = gridV.jStart; j <= gridV.jEnd - 1; j++)
			{
				if (j > gridV.jStart)
				{
					if (levelSet1(i, j) + levelSet1(i, j - 1) <= 0)
					{
						V(i, j) = -We;
						Fluid.originV(i, j) = -We;
					}
				}
			}
		}

		cflCondition = 1.0 / 8.0;
		dt = cflCondition*min(grid.dx, grid.dy);

		singularForce = true;
		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		double finalTime = 10.0;
		maxIteration = ceil(finalTime / dt); // Up to 2 sec.
	}

	ofstream conditionFile;
	conditionFile.open("D:\\Data/Condition.txt", ios::binary);
	conditionFile << "Re = " << Re << endl;
	conditionFile << "Ca = " << Ca << endl;
	conditionFile << "Xi = " << Xi << endl;
	conditionFile << "El = " << El << endl;
	conditionFile << "Pe = " << Pe << endl;
	conditionFile << "cflCondition = " << cflCondition << endl;
	conditionFile << "dx = " << grid.dx << endl;
	conditionFile << "dt = " << dt << endl;
	conditionFile.close();
	

}

inline void InsolubleSurfactant::ContinuumMethodWithSurfactantSolver(const int & example)
{
	bool isPlot = true;
	bool writeFile = false;
	string fileName;
	string str;
	clock_t startTime = clock();
	double  before = 0, after = 0, timeCheck = 0;

	InitialCondition(example);

	grid.Variable();
	levelSet.phi.Variable("phi0");
	levelSet.tube.Variable("Tube");
	Surfactant.Variable("Surfactant0");

	if (isPlot)
	{
		MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
		Surfactant.Variable("SurTube1");
		//MATLAB.Command("subplot(2,1,1)");
		//PlotSurfactant();
		//MATLAB.Command("subplot(2,1,2)");
		PlotVelocity();
		MATLAB.Command("IntSur0 = sum(sum(SurTube1.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");

		//MATLAB.WriteImage("surfactant", iteration, "fig");
		MATLAB.WriteImage("surfactant", iteration, "png");
	}
	


	int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy)) * 2;
	int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.2*min(grid.dx, grid.dy)));
	for (iteration = 1; iteration <= maxIteration; iteration++)
	{
		before = clock();
		cout << "*******************************************************************" << endl;
		cout << "       Iteration " << to_string(iteration) << " : Start" << endl;
		totalT += dt;
		//// Step 1-1 : Surfactant Diffusion
		//cout << "Diffusion Start" << endl;
		InterfaceSurfactant.LSurfactantDiffusion(iteration);
		//cout << "Diffusion End" << endl;
		cout << endl;

		//// Step 1-2 : New Surface Tension
		InterfaceSurfactant.DimlessNonlinearLangmuirEOS(1);
		//InterfaceSurfactant.SurfaceTension.Variable("SurfaceTension");

		//// Step 2 : Navier-Stokes equation
		NSSolver();
		//Pressure.Variable("Pressure");

		//// Step 3 : Level Set Propagation
		AdvectionMethod2D<double>::LLSPropagatingTVDRK3MACGrid(levelSet, U, V, dt, LspatialOrder);
		AdvectionMethod2D<double>::LLSReinitializationTVDRK3(levelSet, dt, reinitialIter, LspatialOrder);

		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, temporalOrder, LspatialOrder, extensionIter);
		levelSet.UpdateInterface();
		levelSet.UpdateLLS();

		InterfaceSurfactant.ConserveSurfactantFactorBeta();

		if (iteration % 10 == 0 && isPlot)
		{
			//MATLAB.Command("subplot(2,1,1)");
			//PlotSurfactant();
			//MATLAB.Command("subplot(2,1,2)");
			PlotVelocity();
			//MATLAB.WriteImage("surfactant", iteration, "fig");
			MATLAB.WriteImage("surfactant", iteration, "png");
			
		}

		timeCheck += ((after = clock()) - before) / CLOCKS_PER_SEC;
		cout << "Consuming Time : " + to_string((after - before) / CLOCKS_PER_SEC) + " / " + to_string(timeCheck) << endl;
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;
		cout << "*******************************************************************" << endl;
	}
}

inline void InsolubleSurfactant::NSSolver()
{
	Fluid.originU.dataArray = U.dataArray;
	Fluid.originV.dataArray = V.dataArray;

	double* uVal(U.dataArray.values);
	double* vVal(V.dataArray.values);
	double* uOriginVal(Fluid.originU.dataArray.values);
	double* vOriginVal(Fluid.originV.dataArray.values);
	int uRes = U.dataArray.ijRes;
	int vRes = V.dataArray.ijRes;

	//// Compute Surface Force
	if (singularForce)
	{
		ComputeSurfaceForce();
	}

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
#pragma omp parallel for
	for (int i = 0; i < uRes; i++)
	{
		uVal[i] = 1. / 4. * (3 * uOriginVal[i] + uVal[i]);
	}
#pragma omp parallel for
	for (int i = 0; i < vRes; i++)
	{
		vVal[i] = 1. / 4. * (3 * vOriginVal[i] + vVal[i]);
	}

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
#pragma omp parallel for
	for (int i = 0; i < uRes; i++)
	{
		uVal[i] = 1. / 3. * (uOriginVal[i] + 2. * uVal[i]);
	}
#pragma omp parallel for
	for (int i = 0; i < vRes; i++)
	{
		vVal[i] = 1. / 3. * (vOriginVal[i] + 2. * vVal[i]);
	}
}

inline void InsolubleSurfactant::EulerMethod()
{
	////////////////////////////////////////////////
	////     Projection Method 1 : advection    ////
	////////////////////////////////////////////////

	EulerMethod1();
	Fluid.TreatBCAlongXaxis(U);
	Fluid.TreatBCAlongYaxis(U);
	Fluid.TreatBCAlongYaxis(V);
	Fluid.TreatBCAlongXaxis(V);
	////////////////////////////////////////////////
	////     Projection Method 2 : Poisson Eq   ////
	////////////////////////////////////////////////
	Fluid.EulerMethodStep2();

	//////////////////////////////////////////////
	////     Projection Method 3 : New U,V    ////
	//////////////////////////////////////////////
	Fluid.EulerMethodStep3();
	Fluid.TreatBCAlongXaxis(U);
	Fluid.TreatBCAlongYaxis(U);
	Fluid.TreatBCAlongYaxis(V);
	Fluid.TreatBCAlongXaxis(V);
}

inline void InsolubleSurfactant::EulerMethod1()
{
	//Array2D<double>& K1U = U.K1;
	//Array2D<double>& K1V = V.K1;
	
	Array2D<int>& UBC = U.BC;
	Array2D<int>& VBC = V.BC;
	Array2D<double>& AdvectionU = Fluid.AdvectionU.dataArray;
	Array2D<double>& AdvectionV = Fluid.AdvectionV.dataArray;
	Array2D<double>& DiffusionU = Fluid.DiffusionU.dataArray;
	Array2D<double>& DiffusionV = Fluid.DiffusionV.dataArray;

	Fluid.AdvectionTerm(U, V, Fluid.AdvectionU, Fluid.AdvectionV);
	Fluid.DiffusionTerm(U, V, Fluid.DiffusionU, Fluid.DiffusionV);
	//AdvectionU.Variable("AdvectionU");
	//DiffusionU.Variable("DiffusionU");
	//AdvectionV.Variable("AdvectionV");
	//DiffusionV.Variable("DiffusionV");
	double oneOverRe = 1 / Re;
#pragma omp parallel for
	for (int i = U.iStart; i <= U.iEnd; i++)
	{
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			if (UBC(i, j)<0) continue;

			//K1U(i, j) = dt*(-AdvectionU(i, j) + oneOverRe*DiffusionU(i, j) + SurfaceForceX(i, j));
			U(i, j) += dt*(- AdvectionU(i, j) + oneOverRe*DiffusionU(i, j) + SurfaceForceX(i, j));
		}
	}
#pragma omp parallel for
	for (int i = V.iStart; i <= V.iEnd; i++)
	{
		for (int j = V.jStart; j <= V.jEnd; j++)
		{
			if (VBC(i, j)<0) continue;
			
			//K1V(i, j) = dt*(-AdvectionV(i, j) + oneOverRe*DiffusionV(i, j) + SurfaceForceY(i, j));
			V(i, j) += dt*(- AdvectionV(i, j) + oneOverRe*DiffusionV(i, j) + SurfaceForceY(i, j));
		}
	}
	//K1U.Variable("k1u");
	//K1V.Variable("k1v");

	//U.Variable("Ustar");
	//V.Variable("Vstar");
	//MATLAB.Command("divUstar =Ustar(:,2:end)-Ustar(:,1:end-1),divVstar =Vstar(2:end,:)-Vstar(1:end-1,:);divstar=divUstar+divVstar;");
	//MATLAB.Command("quiver(Xp,Yp,Ustar(:,1:end-1),Vstar(1:end-1,:)");
}

inline void InsolubleSurfactant::EulerMethod2ndOrder()
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

inline void InsolubleSurfactant::EulerMethod2ndOrder1()
{
}

inline void InsolubleSurfactant::EulerMethod2ndOrder1stIteration1()
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
	GenerateLinearSystemUV(Fluid.Ub, U, Fluid.gradientPx, Fluid.AdvectionU, SurfaceForceX, Fluid.U.grid, 1);
	GenerateLinearSystemUV(Fluid.Vb, V, Fluid.gradientPy, Fluid.AdvectionV, SurfaceForceY, Fluid.V.grid, 1);

	//Fluid.Ub.Variable("Ub");
	//Fluid.Vb.Variable("Vb");

	CGSolver::Solver(Fluid.UCN_CSR, Fluid.Ub, Fluid.tempU);
	CGSolver::Solver(Fluid.VCN_CSR, Fluid.Vb, Fluid.tempV);

	//tempU.Variable("tempU");
	//tempV.Variable("tempV");

	Fluid.VectorToGrid(Fluid.tempU, U);
	Fluid.VectorToGrid(Fluid.tempV, V);

	Fluid.TreatVelocityBC(U, V);

	//U.Variable("Ustar");
	//V.Variable("Vstar");
	//MATLAB.Command("divUstar =Ustar(:,2:end)-Ustar(:,1:end-1),divVstar =Vstar(2:end,:)-Vstar(1:end-1,:);divstar=divUstar+divVstar;");
	//MATLAB.Command("quiver(Xp,Yp,Ustar(:,1:end-1),Vstar(1:end-1,:))");
}

inline void InsolubleSurfactant::ComputeSurfaceForce()
{
	//////////////////////////
	// ST : Surface Tension
	// L : Level Set
	// S : Surfactant
	// G : Gradient
	// U : Unit
	// N : Normal
	//////////////////////////
	Array2D<double>& meanCurvature = levelSet.meanCurvature.dataArray;
	Array2D<VT>& LunitNormal = levelSet.unitNormal.dataArray;
	Array2D<VT>& STgrad = SurfaceTension.gradient;
	Array2D<int>& tube = levelSet.tube;

	SurfaceTension.Gradient();

	int ComputedTubeRange = 1;

	levelSet.LComputeMeanCurvature(ComputedTubeRange);
	levelSet.LComputeUnitNormal(ComputedTubeRange);
	
	double oneOverReCa = 1 / (Re*Ca);
	int numTube = levelSet.numTube;
	double oneOverdx = grid.oneOverdx, oneOverdy = grid.oneOverdy;
	int i, j;
	VT STG, LUN, LG, STSurfaceG;
	double LGMag, deltaL, curvature, ST;
	double STval, STvalLeft, STvalBottom;
	double Lval, LvalLeft, LvalBottom;

#pragma omp parallel for private(i, j, STG, LUN, LG, STSurfaceG, LGMag, deltaL, curvature, ST, STval, STvalLeft, STvalBottom, Lval, LvalLeft, LvalBottom)
	for (int k = 1; k <= numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (tube(i, j) <= ComputedTubeRange)
		{
			STval = SurfaceTension(i, j), STvalLeft = SurfaceTension(i - 1, j), STvalBottom = SurfaceTension(i, j - 1);
			Lval = levelSet(i, j), LvalLeft = levelSet(i - 1, j), LvalBottom = levelSet(i, j - 1);
			/////////////////////////////////
			// SurfaceForceX on MAC grid.  //
			/////////////////////////////////
			LG.x = (Lval - LvalLeft) * oneOverdx;
			LG.y = 0.5 * (levelSet.dyPhi(i - 1, j) + levelSet.dyPhi(i, j));
			LGMag = LG.magnitude();
			LUN = LG / LGMag;
			STG.x = (STval - STvalLeft) * oneOverdx;
			STG.y = 0.5 * (SurfaceTension.dyPhi(i - 1, j) + SurfaceTension.dyPhi(i, j));
			STSurfaceG = STG - DotProduct(LUN, STG)*LUN;
			deltaL = AdvectionMethod2D<double>::DeltaFt(0.5 * (Lval + LvalLeft));
			curvature = 0.5 * (- 2 * meanCurvature(i, j) - 2 * meanCurvature(i - 1, j));
			ST = 0.5 * (STval + STvalLeft);

			SurfaceForceX(i, j) = curvature * ST * LUN.x - STSurfaceG.x;
			SurfaceForceX(i, j) *= - deltaL * LGMag * oneOverReCa;
			
			/////////////////////////////////
			// SurfaceForceY on MAC grid.  //
			/////////////////////////////////
			LG.x = 0.5 * (levelSet.dxPhi(i, j) + levelSet.dxPhi(i, j - 1));
			LG.y = (Lval - LvalBottom) * oneOverdy;
			LGMag = LG.magnitude();
			LUN = LG / LGMag;
			STG.x = 0.5 * (SurfaceTension.dxPhi(i, j) + SurfaceTension.dxPhi(i, j - 1));
			STG.y = (STval - STvalBottom) * oneOverdy;
			STSurfaceG = STG - DotProduct(LUN, STG)*LUN;
			deltaL = AdvectionMethod2D<double>::DeltaFt(0.5 * (Lval + LvalBottom));
			curvature = 0.5 * (- 2 * meanCurvature(i, j) - 2 * meanCurvature(i, j - 1));
			ST = 0.5 * (STval + STvalBottom);

			SurfaceForceY(i, j) = curvature * ST * LUN.y - STSurfaceG.y;
			SurfaceForceY(i, j) *= - deltaL * LGMag * oneOverReCa;
		}
		else
		{
			SurfaceForceX(i, j) = 0;
			SurfaceForceY(i, j) = 0;
		}
	}
	//SurfaceForceX.Variable("ForceX");
	//MATLAB.Command("sumX= sum(sum(ForceX));");
	//SurfaceForceY.Variable("ForceY");
	//MATLAB.Command("sumY= sum(sum(ForceY));");
}

inline void InsolubleSurfactant::GenerateLinearSystemUV(VectorND<double>& vectorB, const FD & vel, const FD & gradP, const FD & advec, const FD& force, const Grid2D & ipGrid, const double & scaling)
{
	int iStart = vel.iStart, iEnd = vel.iEnd, jStart = vel.jStart, jEnd = vel.jEnd;
	double oneOverdx2 = ipGrid.oneOverdx2, oneOverdy2 = ipGrid.oneOverdy2;
	double oneOver2Re = 1. / (2. * Re);
	double dtOver2Re = dt / (2 * Re);
	double* bVal(vectorB.values);

	int index;
#pragma omp parallel for private(index)
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			index = vel.BC(i, j);
			if (index < 0) continue;

			bVal[index] = vel(i, j) + dt*(-gradP(i, j) - advec(i, j) + oneOver2Re*(vel.dxxPhi(i, j) + vel.dyyPhi(i, j)) + force(i, j));
			
			if (i == iStart+1)
			{
				bVal[index] += dtOver2Re * oneOverdx2 * (2 * vel(i - 1, j) - vel.dataArrayOld(i - 1, j));
			}
			if (i == iEnd - 1)
			{
				bVal[index] += dtOver2Re * oneOverdx2 * (2 * vel(i + 1, j) - vel.dataArrayOld(i + 1, j));
			}
			if (j == jStart + 1)
			{
				bVal[index] += dtOver2Re * oneOverdy2 * (2 * vel(i, j - 1) - vel.dataArrayOld(i, j - 1));
			}
			if (j == jEnd - 1)
			{
				bVal[index] += dtOver2Re * oneOverdy2 * (2 * vel(i, j + 1) - vel.dataArrayOld(i, j + 1));
			}
		}
	}
//	int index;
//#pragma omp parallel for private(index)
//	for (int i = innerIStart; i <= innerIEnd; i++)
//	{
//		for (int j = innerJStart; j <= innerJEnd; j++)
//		{
//			index = (i - innerIStart) + (j - innerJStart)*innerIRes;
//			vectorB(index) = vel(i, j) + dt*(-gradP(i, j) - advec(i, j) + 1. / (2. * Re)*(vel.dxxPhi(i, j) + vel.dyyPhi(i, j)) + force(i, j));
//
//			if (i == innerIStart)
//			{
//				vectorB(index) += (2 * vel(i - 1, j) - vel.dataArrayOld(i - 1, j))*ipGrid.oneOverdx2* dt / (2 * Re);
//			}
//			if (i == innerIEnd)
//			{
//				vectorB(index) += (2 * vel(i + 1, j) - vel.dataArrayOld(i + 1, j))*ipGrid.oneOverdx2* dt / (2 * Re);
//			}
//			if (j == innerJStart)
//			{
//				vectorB(index) += (2 * vel(i, j - 1) - vel.dataArrayOld(i, j - 1))*ipGrid.oneOverdy2* dt / (2 * Re);
//			}
//			if (j == innerJEnd)
//			{
//				vectorB(index) += (2 * vel(i, j + 1) - vel.dataArrayOld(i, j + 1))*ipGrid.oneOverdy2* dt / (2 * Re);
//			}
//
//			vectorB(index) *= scaling;
//		}
//	}
}

inline void InsolubleSurfactant::PlotSurfactant()
{
	bool lookDown = true;
	string str;

	Surfactant.Variable("Surfactant");
	levelSet.tube.Variable("Tube");
	MATLAB.Command("SurTube1 = Surfactant.*(Tube<=1);");
	//MATLAB.Command("surf(X,Y,SurTube1), h=colorbar,h.Limits=[0 max(max(SurTube1))],axis equal,axis([X(1) X(end) Y(1) Y(end)]), set(gca,'fontsize',20)");
	if (lookDown)
	{
		MATLAB.Command("surf(X,Y,SurTube1), h=colorbar,h.Limits=[0 max(max(SurTube1))],axis([X(1) X(end) Y(1) Y(end)]),axis equal tight,set(gca,'fontsize',20);");
	}
	else
	{
		MATLAB.Command("surf(X,Y,SurTube1), h=colorbar,h.Limits=[0 max(max(SurTube1))],axis equal tight,set(gca,'fontsize',20);");
	}
	//MATLAB.Command("surf(X,Y,Tube),axis equal,set(gca,'fontsize',20)");

	//// Measure Surfactant Loss ////
	MATLAB.Command("IntSur = sum(sum(SurTube1.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");
	//MATLAB.Command("loss = (IntSur-IntSur0)/IntSur0*100;");
	MATLAB.Variable("i", iteration);
	MATLAB.Variable("totalT", totalT);
	//str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', error(%)  :',num2str(loss),'%']);");
	str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT)]);");
	MATLAB.Command(str.c_str());
}

inline void InsolubleSurfactant::PlotVelocity()
{
	string str;

	U.Variable("U");
	V.Variable("V");
	levelSet.phi.Variable("phi");

	//str = string("quiver(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),set(gca,'fontsize',20);");
	//str = str + string("hold on,streamslice(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,'g'),hold off;set(gca,'fontsize',20);axis equal tight;");
	
	str = str + string("plot(0,0),streamslice(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,'g'),hold off;set(gca,'fontsize',20);axis equal tight;");

	//MATLAB.Command(" [sX,sY]=meshgrid(-1.5:.04:1.5, -0.5:.1:0.5)");
	//str = str + string("plot(0,0),streamline(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,sX,sY),hold off;set(gca,'fontsize',20);axis equal tight;");
	
	str = str + string("hold on, contour(X,Y,phi,[0 0],'r'), grid on,hold off;axis equal tight;");
	MATLAB.Command(str.c_str());
	str = string("title(['iteration : ', num2str(") + to_string(iteration) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	MATLAB.Command(str.c_str());
	//MATLAB.Command("divU =U(:,2:end)-U(:,1:end-1),divV =V(2:end,:)-V(1:end-1,:);div=divU+divV;");
}
