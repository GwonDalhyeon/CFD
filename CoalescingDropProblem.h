#pragma once

#include "AdvectionMethod2D.h"
#include "FluidSolver2D.h"
#include "EulerianMovingInterface.h"

class Coalescence
{
public:
	Coalescence(FluidSolver2D & ipFluid, MovingInterface& ipInterfaceSurfactant);
	~Coalescence();

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

	bool& isCSFmodel = Fluid.isCSFmodel;
	FD& SurfaceForce = Fluid.SurfaceForce;
	FD& SurfaceForceX = Fluid.SurfaceForceX;
	FD& SurfaceForceY = Fluid.SurfaceForceY;
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
	double& finalT = Fluid.finalT;

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
	inline void DropCoalescenceSolver(const int& example);

	inline void PlotSurfactant();
	inline void PlotVelocity();


	/////////////////////////////////////////////////////
	////              CoalescingDrop                 ////
	/////////////////////////////////////////////////////
	int Nf = 1; // The Num of Interfaces per Front. a Liquid Drop : 1, a Soap Bubble : 2
	double& densityI = Fluid.densityI, densityE = Fluid.densityE;
	double densityF = pow(10, 3); // Film Density. 10^3 Kg/m^3.
	double& viscosityI = Fluid.viscosityI, viscosityE = Fluid.viscosityE;
	double densityRatio = 1;// = densityE/densityI;  a Liquid Drop : 0.1, a Soap Bubble : 1
	double viscosityRatio = 1; // = muE/muI;  a Liquid Drop : 0.1, a Soap Bubble : 1

	double lengthscale = 1; // Standard Length Scale
	double timescale = 1; // 

	double& Oh = Fluid.Oh;
	double& We = Fluid.We;
	double& Bo = Fluid.Bo;
	const double& gravity = Fluid.gravity;

	double& BoF = Fluid.BoF; // a Film Bond Number. a Liquid Drop : 0, a Soap Bubble : nonzero.

	double& gamma0 = Fluid.gamma0;

	double& thickness0 = Fluid.thickness0; // the thickness has a typical value of h0 = 10E-6 m.


private:

};

inline Coalescence::Coalescence(FluidSolver2D & ipFluid, MovingInterface & ipInterfaceSurfactant)
	:Fluid(ipFluid), InterfaceSurfactant(ipInterfaceSurfactant)
{
}

Coalescence::~Coalescence()
{
}

inline void Coalescence::InitialCondition(const int & example)
{
	ExamNum = example;

	if (example == 1)
	{
		cout << "*************************" << endl;
		cout << " Dynamics of Drop Coalescence at Fluid Interfaces " << endl;
		cout << "      -- F Blanchette & T Bigioni-- " << endl;
		cout << "*************************" << endl;

		int gridSize = 250;
		double xLength = 2.5;
		grid = Grid2D(-xLength, xLength, gridSize + 1, -xLength, xLength, gridSize + 1);

		Re = 1;

		densityI = 1;
		densityRatio = 10;
		densityE = 1 / densityRatio;

		viscosityI = 1;
		viscosityRatio = 10;
		viscosityE = 1 / viscosityRatio;


		gamma0 = 1; Fluid.gammaWater;

		//// Initialize Surfactant Fields
		levelSet = LS(grid, 3 * grid.dx);
		LS levelSet1(grid, 3 * grid.dx);
		LS levelSet2(grid, 3 * grid.dx);
		double radius = 0.5;
		VT center(0, radius);
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
				levelSet2(i, j) = grid(i, j).y - 0.01*radius;
			}
		}

		levelSet.CombineLevelSet(levelSet1, levelSet2);
		levelSet.InitialTube();
		levelSet.LComputeMeanCurvature(1);

		InterfaceSurfactant.InitialCondition(9);

		//// Initialize Velocity Fields
		Fluid.InitialCondition(7);
		Fluid.isGravity = false;
		Fluid.isGravity = true;
		isCSFmodel = false;
		//We = densityE * lengthscale * lengthscale / gamma0;
		Oh = 0.0025;
		Bo = 0.67;
		BoF = 0;
		
//#pragma omp parallel for
//		for (int i = gridV.iStart; i <= gridV.iEnd; i++)
//		{
//			for (int j = gridV.jStart; j <= gridV.jEnd - 1; j++)
//			{
//				if (j > gridV.jStart)
//				{
//					if (levelSet1(i, j) + levelSet1(i, j - 1) <= 0)
//					{
//						V(i, j) = -We;
//						Fluid.originV(i, j) = -We;
//					}
//				}
//			}
//		}

		SurfaceForce = FD(grid);
		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		finalT = 10;



	}

	if (example == 3 || example == 4)
	{
		cout << "*************************" << endl;
		cout << "        Simulations of surfactant effects  " << endl;
		cout << " on the dynamics of coalescing drops and bubbles" << endl;
		cout << "      --DW Martin and F Blanchette-- " << endl;
		cout << "               Example 1 " << endl;
		cout << "*************************" << endl;

		int gridSize = 250;
		double xLength = 0.0025;
		grid = Grid2D(-xLength, xLength, gridSize + 1, -xLength, xLength, gridSize + 1);

		Re = 1;

		// a Liquid Drop : 1, a Soap Bubble : 2
		if (ExamNum == 3) Nf = 2;
		if (ExamNum == 4) Nf = 1;

		densityE = Fluid.densityAir;
		if (ExamNum == 3) densityI = Fluid.densityAir;
		if (ExamNum == 4) densityI = Fluid.densityWater;
		densityRatio = densityE / densityI; // 1 : Bubble, 0.1 : Liquid drop


		viscosityE = Fluid.viscosityAir;
		if (ExamNum == 3) viscosityI = Fluid.viscosityAir;
		if (ExamNum == 4) viscosityI = Fluid.viscosityWater;
		viscosityRatio = viscosityE / viscosityI; // 1 : Bubble, 0.1 : Liquid drop


		gamma0 = Fluid.gammaWater;

		//timescale = sqrt(densityI*pow(lengthscale, 3) / (Nf*gamma0));

		//// Initialize Surfactant Fields
		levelSet = LS(grid, 3 * grid.dx);
		LS levelSet1(grid, 3 * grid.dx);
		LS levelSet2(grid, 3 * grid.dx);
		double radius = 0.0005;
		VT center(0, radius);
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
				levelSet2(i, j) = grid(i, j).y - 0.01*radius;
			}
		}

		levelSet.CombineLevelSet(levelSet1, levelSet2);
		levelSet.InitialTube();
		levelSet.LComputeMeanCurvature(1);

		InterfaceSurfactant.InitialCondition(9);

		//// Initialize Velocity Fields
		Fluid.InitialCondition(7);

		if (ExamNum == 3) thickness0 = pow(10, -6);

		//lengthscale = 0.001;
		//We = densityE * lengthscale * lengthscale / gamma0;
		//Oh = viscosityI / sqrt(lengthscale*densityI*Nf*gamma0);
		//Bo = radius*radius * gravity * (densityI - densityE) / (Nf*gamma0);
		//BoF = thickness0 * radius * gravity * (densityF - densityE) / (Nf*gamma0);
		//densityI = 1;
		//densityE = densityRatio;
		//viscosityI = 1;
		//viscosityE = viscosityRatio;
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

		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		finalT = 5;
	}

	ofstream conditionFile;
	conditionFile.open("D:\\Data/Condition.txt", ios::binary);
	conditionFile << "Re = " << Re << endl;
	conditionFile << "Ca = " << Ca << endl;
	conditionFile << "Xi = " << Xi << endl;
	conditionFile << "El = " << El << endl;
	conditionFile << "Pe = " << Pe << endl;
	conditionFile << "We = " << We << endl;
	conditionFile << "Oh = " << Oh << endl;
	conditionFile << "Bo = " << Bo << endl;
	conditionFile << "BoF = " << BoF << endl;
	conditionFile << "cflCondition = " << cflCondition << endl;
	conditionFile << "dx = " << grid.dx << endl;
	conditionFile << "finalT = " << finalT << endl;
	conditionFile.close();


}

inline void Coalescence::DropCoalescenceSolver(const int & example)
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
		//Surfactant.Variable("SurTube1");
		//MATLAB.Command("subplot(1,2,1)");
		//PlotSurfactant();
		//MATLAB.Command("subplot(1,2,2)");
		//Fluid.PlotVelocity();
		PlotVelocity();
		MATLAB.Command("IntSur0 = sum(sum(SurTube1.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");

		//MATLAB.WriteImage("Coalescence", iteration, "fig");
		MATLAB.WriteImage("Coalescence", iteration, "png");
	}



	int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy)) * 3;
	int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.2*min(grid.dx, grid.dy)));
	while (totalT < finalT)
	{
		iteration++;
		before = clock();
		cout << "*******************************************************************" << endl;
		cout << "       Iteration " << to_string(iteration) << " : Start" << endl;

		dt = Fluid.AdaptiveTimeStep();
		if (totalT + dt>finalT) dt = finalT - totalT;
		totalT += dt;

		cout << "dt : " + to_string(dt) << endl;
		if (dt<DBL_EPSILON)
		{
			cout << "Blow Up !!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			system("pause");
			break;
		}

		//// Step 1-1 : Surfactant Diffusion
		//cout << "Diffusion Start" << endl;
		//InterfaceSurfactant.LSurfactantDiffusion(iteration);
		//Surfactant.Variable("surfactant");
		//cout << "Diffusion End" << endl;
		//cout << endl;

		//// Step 1-2 : New Surface Tension
		InterfaceSurfactant.DimlessNonlinearLangmuirEOS(1);

		//// Step 2 : Navier-Stokes equation
		Fluid.TVDRK3TimeAdvection();

		//// Step 3 : Level Set Propagation
		AdvectionMethod2D<double>::LLSPropagatingTVDRK3MACGrid(levelSet, U, V, dt, LspatialOrder);
		AdvectionMethod2D<double>::LLSReinitializationTVDRK3(levelSet, 0.5*grid.dx, reinitialIter, LspatialOrder);
		//AdvectionMethod2D<double>::LLSReinitializationTVDRK3SubcellFixSecondOrder(levelSet, 0.5*grid.dx, reinitialIter);
		//AdvectionMethod2D<double>::LLSReinitializationTVDRK3usingSubcellFix(levelSet, 0.5*grid.dx, reinitialIter, LspatialOrder);

		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, temporalOrder, LspatialOrder, extensionIter);
		levelSet.UpdateInterface();
		levelSet.UpdateLLS();

		//InterfaceSurfactant.ConserveSurfactantFactorBeta();

		if (iteration % 25 == 0 && isPlot)
		{
			SurfaceTension.Variable("SurfaceTension");
			Pressure.Variable("Pressure");
			//MATLAB.Command("subplot(1,2,1)");
			//PlotSurfactant();
			//MATLAB.Command("subplot(1,2,2)");
			//Fluid.PlotVelocity();
			PlotVelocity();
			//MATLAB.WriteImage("Coalescence", iteration, "fig");
			MATLAB.WriteImage("Coalescence", iteration, "png");

		}

		timeCheck += ((after = clock()) - before) / CLOCKS_PER_SEC;
		cout << "Consuming Time : " + to_string((after - before) / CLOCKS_PER_SEC) + " / " + to_string(timeCheck) << endl;
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;
		cout << "*******************************************************************" << endl;
	}
}

inline void Coalescence::PlotSurfactant()
{
	bool lookDown = true;
	string str;
	int iter = iteration;
	Surfactant.Variable("Surfactant");
	levelSet.tube.Variable("Tube");
	MATLAB.Command("SurTube1 = Surfactant.*(Tube<=1);");
	//MATLAB.Command("surf(X,Y,SurTube1), h=colorbar,h.Limits=[0 max(max(SurTube1))],axis equal,axis([X(1) X(end) Y(1) Y(end)]), set(gca,'fontsize',20)");
	if (lookDown)
	{
		MATLAB.Command("surf(X,Y,SurTube1), h=colorbar,h.Limits=[min(min(SurTube1)) max(max(SurTube1))],axis([X(1) X(end) Y(1) Y(end)]),axis equal tight,set(gca,'fontsize',20);");
	}
	else
	{
		MATLAB.Command("surf(X,Y,SurTube1), h=colorbar,h.Limits=[min(min(SurTube1)) max(max(SurTube1))],axis equal tight,set(gca,'fontsize',20);");
	}
	//MATLAB.Command("surf(X,Y,Tube),axis equal,set(gca,'fontsize',20)");

	//// Measure Surfactant Loss ////
	MATLAB.Command("IntSur = sum(sum(SurTube1.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");
	//MATLAB.Command("loss = (IntSur-IntSur0)/IntSur0*100;");
	MATLAB.Variable("i", iter);
	MATLAB.Variable("totalT", totalT);
	//str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', error(%)  :',num2str(loss),'%']);");
	str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT)]);");
	MATLAB.Command(str.c_str());
}

inline void Coalescence::PlotVelocity()
{
	string str;
	int iter = iteration;

	U.Variable("U");
	V.Variable("V");
	levelSet.phi.Variable("phi");

	str = string("quiver(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),set(gca,'fontsize',20);");
	str = str + string("hold on,streamslice(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,'g'),plot([X(1) X(end)],[1 1],'g'),hold off;set(gca,'fontsize',20);");

	//str = str + string("plot(0,0),streamslice(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,'g'),hold off;set(gca,'fontsize',20);axis equal tight;");

	//MATLAB.Command(" [sX,sY]=meshgrid(-1.5:.04:1.5, -0.5:.1:0.5)");
	//str = str + string("plot(0,0),streamline(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,sX,sY),hold off;set(gca,'fontsize',20);axis equal tight;");

	str = str + string("hold on, contour(X,Y,phi,[0 0],'r'), grid on,hold off;axis equal;axis([-1.5 1.5 -0.5 1.5]);");
	MATLAB.Command(str.c_str());
	str = string("title(['iteration : ', num2str(") + to_string(iteration) + string("),', time : ', num2str(") + to_string(totalT) + string("),', dt : ', num2str(") + to_string(dt) + string(")]);");
	MATLAB.Command(str.c_str());
	//MATLAB.Command("divU =U(:,2:end)-U(:,1:end-1),divV =V(2:end,:)-V(1:end-1,:);div=divU+divV;");
}
