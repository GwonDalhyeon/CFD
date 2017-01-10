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
	int iteration;
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

	bool singularForce = false;
	FD SurfaceForceX;
	FD SurfaceForceY;
	FV SurfGradSurfTension;


	int& ProjectionOrder = Fluid.ProjectionOrder;
	double& reynoldNum = Fluid.reynoldNum;
	double& cflCondition = Fluid.cflCondition;
	double& Ca = Fluid.Ca;
	double& Xi = Fluid.Xi;
	double& El = Fluid.El;
	double& Pe = Fluid.Pe;

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
		reynoldNum = 1000;
		Ca = 0.5;
		Xi = 0.3;
		El = 0.2;
		Pe = 10;
		cflCondition = 0.25;
		dt = cflCondition*min(grid.dx, grid.dy);

		singularForce = true;
		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		// Initialize Surfactant Fields
		InterfaceSurfactant.InitialCondition(7);
		InterfaceSurfactant.cflCondition = cflCondition;
		InterfaceSurfactant.dt = dt;

		double finalT = 5;// Up to 2 sec.
		maxIteration = ceil(finalT / dt);
		totalT = 0;
		writeOutputIteration = 30;
		iteration = 0;
		totalT = 0;
	}
}

inline void InsolubleSurfactant::ContinuumMethodWithSurfactantSolver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;


	InitialCondition(example);

	grid.Variable();
	levelSet.phi.Variable("phi0");
	levelSet.tube.Variable("Tube");
	Surfactant.Variable("Surfactant0");

	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
	Surfactant.Variable("SurTube1");
	PlotVelocity();
	MATLAB.Command("IntSur0 = sum(sum(SurTube1.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");

	//MATLAB.WriteImage("surfactant", 0, "fig");
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
		//InterfaceSurfactant.LSurfactantDiffusion(iteration);
		cout << "Diffusion End" << endl;
		cout << endl;

		//// Step 1-2 : New Surface Tension
		//InterfaceSurfactant.DimlessNonlinearLangmu1rEOS(2);
		//InterfaceSurfactant.SurfaceTension.Variable("SurfaceTension");

		//// Step 2 : Navier-Stokes equation
		NSSolver();
		Pressure.Variable("Pressure");


		AdvectionMethod2D<double>::LLSPropagatingTVDRK3(levelSet, U, V, dt);
		AdvectionMethod2D<double>::LLSReinitializationTVDRK3(levelSet, dt, reinitialIter);

		InterfaceSurfactant.ConserveSurfactantFactorBeta();

		AdvectionMethod2D<double>::LLSQuantityExtension(levelSet, Surfactant, 3, 3, extensionIter);
		levelSet.UpdateInterface();
		levelSet.UpdateLLS();


		//MATLAB.Command("subplot(2,1,1)");
		//PlotSurfactant();
		//MATLAB.Command("subplot(2,1,2)");
		
		if (iteration == 1 || iteration % 1 == 0)
		{
			PlotVelocity();
			//MATLAB.WriteImage("surfactant", iteration, "fig");
			//MATLAB.WriteImage("surfactant", iteration, "png");
		}
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;
		cout << "*******************************************************************" << endl;
	}
}

inline void InsolubleSurfactant::NSSolver()
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
	//// Compute Surface Force
	if (singularForce)
	{
		ComputeSurfaceForce();
	}

	Array2D<double>& K1U = U.K1;
	Array2D<double>& K1V = V.K1;
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
	double oneOverRe = 1 / reynoldNum;
#pragma omp parallel for
	for (int i = K1U.iStart; i <= K1U.iEnd; i++)
	{
		for (int j = K1U.jStart; j <= K1U.jEnd; j++)
		{
			if (UBC(i, j)<0)
			{
				continue;
			}
			K1U(i, j) = dt*(-AdvectionU(i, j) + oneOverRe*DiffusionU(i, j) + SurfaceForceX(i, j));
			U(i, j) = U(i, j) + K1U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = K1V.iStart; i <= K1V.iEnd; i++)
	{
		for (int j = K1V.jStart; j <= K1V.jEnd; j++)
		{
			if (VBC(i, j)<0)
			{
				continue;
			}
			K1V(i, j) = dt*(-AdvectionV(i, j) + oneOverRe*DiffusionV(i, j) + SurfaceForceY(i, j));
			V(i, j) = V(i, j) + K1V(i, j);
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
	GenerateLinearSystemUV(Fluid.Ub, U, Fluid.gradientPx, Fluid.AdvectionU, SurfaceForceX, Fluid.U.innerGrid, 1);
	GenerateLinearSystemUV(Fluid.Vb, V, Fluid.gradientPy, Fluid.AdvectionV, SurfaceForceY, Fluid.V.innerGrid, 1);

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
	FD& meanCurvature = levelSet.meanCurvature;
	FV& unitNormal = levelSet.unitNormal;
	Array2D<VT>& gradient = Surfactant.gradient;

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
			gradient(i, j) = VT(Surfactant.dxPhi(i, j), Surfactant.dyPhi(i, j));
			SurfGradSurfTension(i, j) = Surfactant.gradient(i, j)
				- DotProduct(unitNormal(i, j), gradient(i, j))*unitNormal(i, j);

			SurfaceForceX(i, j) = -2 * meanCurvature(i, j)*SurfaceTension(i, j)*unitNormal(i, j).x + SurfGradSurfTension(i, j).x;
			SurfaceForceX(i, j) *= AdvectionMethod2D<double>::DeltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude() / (reynoldNum*Ca);
			SurfaceForceY(i, j) = -2 * meanCurvature(i, j)*SurfaceTension(i, j)*unitNormal(i, j).y + SurfGradSurfTension(i, j).y;
			SurfaceForceY(i, j) *= -AdvectionMethod2D<double>::DeltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude() / (reynoldNum*Ca);

		}
		else
		{
			SurfGradSurfTension(i, j) = 0;

			SurfaceForceX(i, j) = 0;
			SurfaceForceY(i, j) = 0;
		}
	}
}

inline void InsolubleSurfactant::GenerateLinearSystemUV(VectorND<double>& vectorB, const FD & vel, const FD & gradP, const FD & advec, const FD& force, const Grid2D & ipGrid, const double & scaling)
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
			vectorB(index) = vel(i, j) + dt*(-gradP(i, j) - advec(i, j) + 1. / (2. * reynoldNum)*(vel.dxxPhi(i, j) + vel.dyyPhi(i, j)) + force(i, j));

			if (i == innerIStart)
			{
				vectorB(index) += (2 * vel(i - 1, j) - vel.dataArrayOld(i - 1, j))*ipGrid.oneOverdx2* dt / (2 * reynoldNum);
			}
			if (i == innerIEnd)
			{
				vectorB(index) += (2 * vel(i + 1, j) - vel.dataArrayOld(i + 1, j))*ipGrid.oneOverdx2* dt / (2 * reynoldNum);
			}
			if (j == innerJStart)
			{
				vectorB(index) += (2 * vel(i, j - 1) - vel.dataArrayOld(i, j - 1))*ipGrid.oneOverdy2* dt / (2 * reynoldNum);
			}
			if (j == innerJEnd)
			{
				vectorB(index) += (2 * vel(i, j + 1) - vel.dataArrayOld(i, j + 1))*ipGrid.oneOverdy2* dt / (2 * reynoldNum);
			}

			vectorB(index) *= scaling;
		}
	}
}

inline void InsolubleSurfactant::PlotSurfactant()
{
	string str;

	Surfactant.Variable("Surfactant");
	levelSet.tube.Variable("Tube");
	MATLAB.Command("SurTube1 = Surfactant.*(Tube<=1);");
	//MATLAB.Command("contour(X,Y,Tube,'r'),axis equal,axis([X(1) X(end) Y(1) Y(end)]), hold on,surf(X,Y,SurTube1), h=colorbar,h.Limits=[0 max(max(SurTube1))],hold off;set(gca,'fontsize',20)");
	//MATLAB.Command("surf(X,Y,SurTube1), h=colorbar,h.Limits=[0 max(max(SurTube1))],axis equal,axis([X(1) X(end) Y(1) Y(end)]), set(gca,'fontsize',20)");
	MATLAB.Command("surf(X,Y,Surfactant), h=colorbar,h.Limits=[0 max(max(SurTube1))],axis equal,set(gca,'fontsize',20)");

	//// Measure Surfactant Loss ////
	MATLAB.Command("IntSur = sum(sum(SurTube1.*(Tube==1)))*(Y(2)-Y(1))*(Y(2)-Y(1));");
	MATLAB.Command("loss = (IntSur-IntSur0)/IntSur0*100;");
	MATLAB.Variable("i", iteration);
	MATLAB.Variable("totalT", totalT);
	str = string("title(['iteration : ', num2str(i),', time : ', num2str(totalT), ', error(%)  :',num2str(loss),'%']);");
	MATLAB.Command(str.c_str());
}

inline void InsolubleSurfactant::PlotVelocity()
{
	string str;

	U.Variable("U");
	V.Variable("V");
	levelSet.phi.Variable("phi");

	str = string("quiver(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),set(gca,'fontsize',20);");

	str = str + string("hold on,streamslice(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,'g'),hold off;");
	//str = str + string("hold on,streamline(X,Y,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,-100:0.1:100,-100:0.1:100),hold off;");
	str = str + string("hold on, contour(X,Y,phi,[0 0],'r'), hold off;axis tight,axis equal;");
	MATLAB.Command(str.c_str());
	str = string("title(['iteration : ', num2str(") + to_string(iteration) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	MATLAB.Command(str.c_str());
	MATLAB.Command("divU =U(:,2:end)-U(:,1:end-1),divV =V(2:end,:)-V(1:end-1,:);div=divU+divV;");
}
