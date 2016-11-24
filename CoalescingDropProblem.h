#pragma once
#include "EulerianFluidSolver.h"
#include "EulerianMovingInterface.h"

class CoalescingDrop
{
public:
	CoalescingDrop(EulerianFluidSolver2D & ipFluid, MovingInterface& ipInterfaceSurfactant);
	~CoalescingDrop();

	int ExamNum;
	int iteration;


	int Nf; // The Num of Interfaces per Front.
	double rhoI, rhoE;
	double muI, muE;
	double densityRatio = rhoE / rhoI; // 1 : Bubble, 0.1 : Liquid drop
	double viscosityRatio = muE / muI; // 1 : Bubble, 0.1 : Liquid drop


	///////////////////////////////
	//// Incompressible Fluid  ////
	///////////////////////////////
	EulerianFluidSolver2D& Fluid;

	Grid2D& grid = Fluid.gridP;
	Grid2D& gridU = Fluid.gridU;
	Grid2D& gridV = Fluid.gridV;
	FD& U = Fluid.U;
	FD& V = Fluid.V;

	bool singularForce = false;
	FD SurfaceForceX;
	FD SurfaceForceY;
	FV SurfGradSurfTension;

	int& accuracyOrder = Fluid.accuracyOrder;
	double& reynoldNum = Fluid.reynoldNum;
	double& cflCondition = Fluid.cflCondition;
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
	double totalT;


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

	inline void PlotSurfactant();
	inline void PlotVelocity();
private:

};

CoalescingDrop::CoalescingDrop(EulerianFluidSolver2D & ipFluid, MovingInterface& ipInterfaceSurfactant)
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
		
		int gridSize = 150;
		grid = Grid2D(-3, 3, gridSize + 1, 0, 5, gridSize * 5. / 3. + 1);

		// Initialize Velocity Fields
		Fluid.InitialCondition(5);
		Fluid.CGsolverNum = 2;
		accuracyOrder = 2;
		reynoldNum = 10;
		rhoE = 1; rhoI = 1;
		muE = 1; muI = 1;

		cflCondition = 1.0 / 8.0;
		dt = cflCondition*min(grid.dx, grid.dy);

		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		// Initialize Surfactant Fields
		InterfaceSurfactant.InitialCondition(8);
		InterfaceSurfactant.CGsolverNum = 1;
		InterfaceSurfactant.cflCondition = cflCondition;
		InterfaceSurfactant.dt = dt;

		maxIteration = ceil(2.0 / dt); // Up to 2 sec.
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
	MATLAB.Command("SurTube1 = Surfactant0.*(Tube<=1);");
	MATLAB.Command("surf(X,Y,SurTube1),grid on,axis equal,set(gca,'fontsize',20)");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	MATLAB.Command(str.c_str());
	MATLAB.Command("IntSur0 = sum(sum(Surfactant0.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");

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
		InterfaceSurfactant.DimlessNonlinearLangmuirEOS(2);
		InterfaceSurfactant.SurfaceTension.Variable("SurfaceTension");

		//// Step 2 : Navier-Stokes equation
		NSSolver();



		AdvectionMethod2D<double>::LLSPropagatingTVDRK3(levelSet, U, V, dt);
		AdvectionMethod2D<double>::LLSReinitializationTVDRK3(levelSet, dt, reinitialIter);

		InterfaceSurfactant.ConserveSurfactantFactorBeta();

		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);
		levelSet.UpdateInterface();
		levelSet.UpdateLLS();


		MATLAB.Command("subplot(2,1,1)");
		PlotSurfactant();
		MATLAB.Command("subplot(2,1,2)");
		PlotVelocity();
		if (iteration == 1 || iteration % 1 == 0)
		{
			MATLAB.WriteImage("surfactant", iteration, "fig");
			MATLAB.WriteImage("surfactant", iteration, "png");
		}
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;

	}
}

inline void CoalescingDrop::NSSolver()
{
}

inline void CoalescingDrop::EulerMethod()
{
}

inline void CoalescingDrop::EulerMethod1()
{
}

inline void CoalescingDrop::EulerMethod2ndOrder()
{
}

inline void CoalescingDrop::EulerMethod2ndOrder1()
{
}

inline void CoalescingDrop::EulerMethod2ndOrder1stIteration1()
{
}

inline void CoalescingDrop::PlotSurfactant()
{
	string str;
	const char* cmd;

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





