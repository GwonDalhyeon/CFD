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

	bool isSurfaceForce = true;
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
	inline void ContinuumMethodWithSurfactantSolver(const int& example);

	inline void PlotSurfactant();
	inline void PlotVelocity();


	/////////////////////////////////////////////////////
	////              CoalescingDrop                 ////
	/////////////////////////////////////////////////////
	int Nf = 1; // The Num of Interfaces per Front. a Liquid Drop : 1, a Soap Bubble : 2
	double& densityI = Fluid.densityI;
	double& densityE = Fluid.densityE;
	double densityF = pow(10, 3); // Film Density. 10^3 Kg/m^3.
	double& viscosityI = Fluid.viscosityI;
	double& viscosityE = Fluid.viscosityE;
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
		Re = 10;
		Ca = 0.5;
		Xi = 0.3;
		El = 0.5;
		Pe = 10;
		cflCondition = 0.4;
		dt = cflCondition*min(grid.dx, grid.dy);

		//isSurfaceForce = false;
		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		// Initialize Surfactant Fields
		InterfaceSurfactant.InitialCondition(7);
		InterfaceSurfactant.cflCondition = cflCondition;
		InterfaceSurfactant.dt = dt;

		Fluid.isCSFmodel = false;
		Fluid.isDeltaFunction = false;
		Fluid.dimensionlessForm = false;
		Fluid.isMultiPhase = true;
		Fluid.isGravity = false;

		finalT = 5;// Up to 5 sec.
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
		double xLength = 2, yLength = 1;
		grid = Grid2D(-xLength, xLength, xLength * gridSize + 1, -yLength, yLength, yLength*gridSize + 1);

		// Initialize Velocity Fields
		//Fluid.isPCG = true;
		Fluid.InitialCondition(3);
		ProjectionOrder = 1;
		Re = 10;
		Ca = 0.1;
		Xi = 0.1;
		El = 0.2;
		Pe = 10;
		cflCondition = 0.2;
		dt = cflCondition*min(grid.dx, grid.dy);
		
		//isSurfaceForce = false;
		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		// Initialize Surfactant Fields
		InterfaceSurfactant.InitialCondition(8);
		InterfaceSurfactant.cflCondition = cflCondition;
		InterfaceSurfactant.dt = dt;

		Fluid.isCSFmodel = true;
		Fluid.isDeltaFunction = false;
		Fluid.dimensionlessForm = false;
		Fluid.isMultiPhase = true;
		Fluid.isGravity = false;
		

		finalT = 10;// Up to 5 sec.
		writeOutputIteration = 30;
	}

	ofstream conditionFile;
	conditionFile.open("D:\\Data/Condition.txt", ios::binary);
	conditionFile << "csf = " << Fluid.isCSFmodel << endl;
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
		//Surfactant.Variable("SurTube1");
		//MATLAB.Command("subplot(1,2,1)");
		//PlotSurfactant();
		//MATLAB.Command("subplot(1,2,2)");
		//Fluid.PlotVelocity();
		PlotVelocity();
		MATLAB.Command("IntSur0 = sum(sum(SurTube1.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");

		//MATLAB.WriteImage("surfactant", iteration, "fig");
		MATLAB.WriteImage("surfactant", iteration, "png");
	}

	int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy)) * 3;
	int extensionIter = (int)ceil((levelSet.gamma2 - levelSet.gamma1) / (0.5*min(grid.dx, grid.dy)));
	//U.dataArray = 0;
	//V.dataArray = 0;
	while (totalT < finalT)
	{
		iteration++;
		Fluid.iteration = iteration;
		before = clock();
		cout << "*******************************************************************" << endl;
		cout << "       Iteration " << to_string(iteration) << " : Start" << endl;
		
		dt = Fluid.AdaptiveTimeStep();
		if (totalT + dt > finalT) dt = finalT - totalT + DBL_EPSILON;
		totalT += dt;

		cout << "dt : " + to_string(dt) << endl;
		if (dt<DBL_EPSILON)
		{
			cout << "Blow Up !!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;
			system("pause");
			break;
		}
		
		//// Step 1-1 : Surfactant Diffusion
		InterfaceSurfactant.LSurfactantDiffusion(iteration);

		//// Step 1-2 : New Surface Tension
		InterfaceSurfactant.DimlessNonlinearLangmuirEOS(1);

		Fluid.DetermineViscosity();

		//// Step 2 : Level Set Propagation
		AdvectionMethod2D<double>::LLSPropagatingTVDRK3MACGrid(levelSet, U, V, dt, LspatialOrder);
		AdvectionMethod2D<double>::LLSReinitializationTVDRK3(levelSet, 0.5*grid.dx, reinitialIter, LspatialOrder);
		//AdvectionMethod2D<double>::LLSReinitializationTVDRK3SubcellFixSecondOrder(levelSet, 0.5*grid.dx, reinitialIter);
		//AdvectionMethod2D<double>::LLSReinitializationTVDRK3usingSubcellFix(levelSet, 0.5*grid.dx, reinitialIter, LspatialOrder);

		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, temporalOrder, LspatialOrder, extensionIter);
		levelSet.UpdateInterface();
		levelSet.UpdateLLS();

		levelSet.LComputeMeanCurvature(1);
		Fluid.DetermineDensity();

		//// Step 3 : Navier-Stokes equation
		Fluid.TVDRK3TimeAdvection();

		//InterfaceSurfactant.ConserveSurfactantFactorBeta();

		if (iteration % 10 == 0 && isPlot)
		{
			//SurfaceTension.Variable("SurfaceTension");
			//Pressure.Variable("Pressure");
			MATLAB.Command("subplot(2,1,1)");
			PlotSurfactant();
			MATLAB.Command("subplot(2,1,2)");
			PlotVelocity();
			//MATLAB.Command("subplot(2,2,[3 4])");
			//MATLAB.Command("contour(X, Y, phi, [0 0], 'r'), grid on, hold off; axis equal, axis([X(1) X(end) Y(1) Y(end)]);");
			//MATLAB.Command("hold on,quiver(X, Y, surfaceGradientX.*(abs(phi) <= Y(2)-Y(1)), surfaceGradientY.*(abs(phi) <= Y(2)-Y(1))), hold off");
			//MATLAB.WriteImage("surfactant", iteration, "fig");
			MATLAB.WriteImage("surfactant", iteration, "png");
		}

		timeCheck += ((after = clock()) - before) / CLOCKS_PER_SEC;
		cout << "Consuming Time : " + to_string((after - before) / CLOCKS_PER_SEC) + " / " + to_string(timeCheck) << endl;
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;
		cout << "*******************************************************************" << endl;
	}
	//SurfaceTension.Variable("SurfaceTension");
	//Pressure.Variable("Pressure");
	MATLAB.Command("subplot(2,1,1)");
	PlotSurfactant();
	MATLAB.Command("subplot(2,1,2)");
	PlotVelocity();
	//MATLAB.Command("subplot(2,2,[3 4])");
	//MATLAB.Command("contour(X, Y, phi, [0 0], 'r'), grid on, hold off; axis equal, axis([X(1) X(end) Y(1) Y(end)]);");
	//MATLAB.Command("hold on,quiver(X, Y, surfaceGradientX.*(abs(phi) <= Y(2)-Y(1)), surfaceGradientY.*(abs(phi) <= Y(2)-Y(1))), hold off");
	//MATLAB.WriteImage("surfactant", iteration, "fig");
	MATLAB.WriteImage("surfactant", iteration, "png");
	levelSet.phi.WriteFile("phi_my");
}

inline void InsolubleSurfactant::PlotSurfactant()
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
		MATLAB.Command("surf(X,Y,SurTube1, 'EdgeColor','none'),h=colorbar,h.Limits=[min(min(SurTube1(Tube==1))) max(max(SurTube1(Tube==1)))],axis([X(1) X(end) Y(1) Y(end)]),axis equal tight,set(gca,'fontsize',20);");
	}
	else
	{
		MATLAB.Command("surf(X,Y,SurTube1),h=colorbar,h.Limits=[min(min(SurTube1(Tube==1))) max(max(SurTube1(Tube==1)))],axis equal tight,set(gca,'fontsize',20);");
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

inline void InsolubleSurfactant::PlotVelocity()
{
	string str; 
	int iter = iteration;

	U.Variable("U");
	V.Variable("V");
	levelSet.phi.Variable("phi");

	str = string("quiver(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),set(gca,'fontsize',20);");
	str = str + string("hold on,streamslice(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,'g'),hold off;set(gca,'fontsize',20);axis equal tight;");
	
	//str = str + string("plot(0,0),streamslice(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,'g'),hold off;set(gca,'fontsize',20);axis equal tight;");

	//MATLAB.Command(" [sX,sY]=meshgrid(-1.5:.04:1.5, -0.5:.1:0.5)");
	//str = str + string("plot(0,0),streamline(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,sX,sY),hold off;set(gca,'fontsize',20);axis equal tight;");
	
	str = str + string("hold on, contour(X,Y,phi,[0 0],'r'), grid on,hold off;axis equal tight;");
	MATLAB.Command(str.c_str());
	str = string("title(['iteration : ', num2str(") + to_string(iteration) + string("),', time : ', num2str(") + to_string(totalT) + string("),', dt : ', num2str(") + to_string(dt) + string(")]);");
	MATLAB.Command(str.c_str());
	//MATLAB.Command("divU =U(:,2:end)-U(:,1:end-1),divV =V(2:end,:)-V(1:end-1,:);div=divU+divV;");
}
