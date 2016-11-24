#pragma once

#include "AdvectionMethod2D.h"
#include "EulerianFluidSolver.h"
#include "EulerianMovingInterface.h"

class InsolubleSurfactant
{
public:
	InsolubleSurfactant(EulerianFluidSolver2D & ipFluid, MovingInterface& ipInterfaceSurfactant);
	~InsolubleSurfactant();

	int ExamNum;
	int iteration;
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

inline InsolubleSurfactant::InsolubleSurfactant(EulerianFluidSolver2D & ipFluid, MovingInterface & ipInterfaceSurfactant)
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
		cout << "*************************" << endl;
		cout << "        A level-set continuum method  " << endl;
		cout << " for two-phase flows with insoluble surfactant" << endl;
		cout << "      --JJ Xu, Y Yang, J Lowengrub-- " << endl;
		cout << "               Example 1 " << endl;
		cout << "*************************" << endl;

		int gridSize = 200;
		grid = Grid2D(-5, 5, gridSize + 1, -2, 2, gridSize * 2. / 5. + 1);
		
		// Initialize Velocity Fields
		Fluid.InitialCondition(3);
		Fluid.CGsolverNum = 2;
		accuracyOrder = 2;
		reynoldNum = 10;
		Ca = 0.5;
		Xi = 0.3;
		El = 0.2;
		Pe = 10;
		cflCondition = 1.0 / 8.0;
		dt = cflCondition*min(grid.dx, grid.dy);

		singularForce = true;
		SurfaceForceX = FD(Fluid.gridU);
		SurfaceForceY = FD(Fluid.gridV);
		SurfGradSurfTension = FV(grid);

		// Initialize Surfactant Fields
		InterfaceSurfactant.InitialCondition(7);
		InterfaceSurfactant.CGsolverNum = 1;		
		InterfaceSurfactant.cflCondition = cflCondition;
		InterfaceSurfactant.dt = dt;
		
		maxIteration = ceil(2.0/dt); // Up to 2 sec.
		totalT = 0;
	}
}

inline void InsolubleSurfactant::ContinuumMethodWithSurfactantSolver(const int & example)
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
		if (iteration ==1 || iteration%1==0)
		{
			MATLAB.WriteImage("surfactant", iteration, "fig");
			MATLAB.WriteImage("surfactant", iteration, "png");
		}
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;
	
	}
}

inline void InsolubleSurfactant::NSSolver()
{


	if (iteration == 1)
	{
		if (accuracyOrder==1)
		{
			Fluid.vectorB = VectorND<double>(Fluid.gridPinner.iRes*Fluid.gridPinner.jRes);
			Fluid.tempP = VectorND<double>(Fluid.gridPinner.iRes*Fluid.gridPinner.jRes);
			Fluid.poissonMatrix = Array2D<double>(1, Fluid.gridPinner.iRes*Fluid.gridPinner.jRes, 1, Fluid.gridPinner.iRes*Fluid.gridPinner.jRes);
			
			Fluid.GenerateLinearSystem(Fluid.poissonMatrix, -Fluid.gridP.dx*Fluid.gridP.dx);
			if (Fluid.CGsolverNum == 1)
			{
				//// CG solver 1
				Fluid.poissonCSR = CSR<double>(Fluid.poissonMatrix);
			}
			else
			{
				//// CG solver 2
				CGSolver::SparseA(Fluid.poissonMatrix, Fluid.a, Fluid.row, Fluid.col, Fluid.nonzeroNum);
			}
			
		}
		else
		{
			Fluid.GenerateLinearSystemUV(Fluid.UCNMatrix, Fluid.gridUinner, 1);
			Fluid.GenerateLinearSystemUV(Fluid.VCNMatrix, Fluid.gridVinner, 1);
			Fluid.GenerateLinearSystemPhi(Fluid.PhiCNMatrix, -Fluid.gridP.dx2 / dt);
			if (Fluid.CGsolverNum == 1)
			{
				//// CG solver 1
				Fluid.UCN_CSR = CSR<double>(Fluid.UCNMatrix);
				Fluid.VCN_CSR = CSR<double>(Fluid.VCNMatrix);
				Fluid.PhiCN_CSR = CSR<double>(Fluid.PhiCNMatrix);
			}
			else if (Fluid.CGsolverNum == 2)
			{
				//// CG solver 2
				CGSolver::SparseA(Fluid.UCNMatrix, Fluid.Ua, Fluid.Urow, Fluid.Ucol, Fluid.UnonzeroNum);
				CGSolver::SparseA(Fluid.VCNMatrix, Fluid.Va, Fluid.Vrow, Fluid.Vcol, Fluid.VnonzeroNum);
				CGSolver::SparseA(Fluid.PhiCNMatrix, Fluid.Phia, Fluid.Phirow, Fluid.Phicol, Fluid.PhinonzeroNum);
			}
		}
		
		

	}

	Fluid.originU.dataArray = U.dataArray;
	Fluid.originV.dataArray = V.dataArray;
	
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

inline void InsolubleSurfactant::EulerMethod()
{
	////////////////////////////////////////////////
	////     Projection Method 1 : advection    ////
	////////////////////////////////////////////////

	EulerMethod1();

	////////////////////////////////////////////////
	////     Projection Method 2 : Poisson Eq   ////
	////////////////////////////////////////////////
	Fluid.EulerMethod2();

	//////////////////////////////////////////////
	////     Projection Method 3 : New U,V    ////
	//////////////////////////////////////////////
	Fluid.EulerMethod3();
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
	Fluid.AdvectionTerm(U, V, Fluid.advectionU, Fluid.advectionV);
	Fluid.DiffusionTerm(U, V, Fluid.diffusionU, Fluid.diffusionV);
	//advectionU.Variable("advectionU");
	//diffusionU.Variable("diffusionU");
	//advectionV.Variable("advectionV");
	//diffusionV.Variable("diffusionV");

#pragma omp parallel for
	for (int i = K1U.iStart; i <= K1U.iEnd; i++)
	{
		for (int j = K1U.jStart; j <= K1U.jEnd; j++)
		{
			K1U(i, j) = dt*(-Fluid.advectionU(i, j) + 1. / reynoldNum*Fluid.diffusionU(i, j) + SurfaceForceX(i, j));
			U(i, j) = U(i, j) + K1U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = K1V.iStart; i <= K1V.iEnd; i++)
	{
		for (int j = K1V.jStart; j <= K1V.jEnd; j++)
		{
			K1V(i, j) = dt*(-Fluid.advectionV(i, j) + 1. / reynoldNum*Fluid.diffusionV(i, j) + SurfaceForceY(i, j));
			V(i, j) = V(i, j) + K1V(i, j);
		}
	}
	//K1U.Variable("k1u");
	//K1V.Variable("k1v");
	//// Boundary : Linear extension.
	Fluid.BdryCondVel();


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
	//advectionV.Variable("advectionV");
	//MATLAB.Command("figure(2),subplot(1,2,1),plot(Vnew(51,2:end-1),'-o'),grid on, subplot(1,2,2),plot(advectionV(end,:),'-o'),grid on");
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


	Fluid.AdvectionTerm(U, V, Fluid.advectionU, Fluid.advectionV);
	Fluid.DiffusionTerm(U, V, Fluid.diffusionU, Fluid.diffusionV);

	double viscosity = 1;
	//// 2nd-order Adams-Bashforth formula
	//Fluid.advectionU.Variable("advectionU");
	//Fluid.diffusionU.Variable("diffusionU");

	//Fluid.advectionV.Variable("advectionV");
	//Fluid.diffusionV.Variable("diffusionV");
	
	//// Crank-Nicolson
	GenerateLinearSystemUV(Fluid.Ub, U, Fluid.gradientPx, Fluid.advectionU, SurfaceForceX, Fluid.gridUinner, 1);
	GenerateLinearSystemUV(Fluid.Vb, V, Fluid.gradientPy, Fluid.advectionV, SurfaceForceY, Fluid.gridVinner, 1);

	//Fluid.Ub.Variable("Ub");
	//Fluid.Vb.Variable("Vb");

	if (Fluid.CGsolverNum == 1)
	{
		Fluid.tempU = CGSolver::SolverCSR(Fluid.UCN_CSR, Fluid.Ub, Fluid.gridU.dx*Fluid.gridU.dy);
		Fluid.tempV = CGSolver::SolverCSR(Fluid.VCN_CSR, Fluid.Vb, Fluid.gridV.dx*Fluid.gridV.dy);
	}
	else if (Fluid.CGsolverNum == 2)
	{
		CGSolver::SolverSparse(Fluid.UCNMatrix.iRes, Fluid.Ua, Fluid.Urow, Fluid.Ucol, Fluid.Ub, Fluid.tempU);
		CGSolver::SolverSparse(Fluid.VCNMatrix.iRes, Fluid.Va, Fluid.Vrow, Fluid.Vcol, Fluid.Vb, Fluid.tempV);
	}
	//tempU.Variable("tempU");
	//tempV.Variable("tempV");

	int index;
#pragma omp parallel for private(index)
	for (int i = Fluid.gridUinner.iStart; i <= Fluid.gridUinner.iEnd; i++)
	{
		for (int j = Fluid.gridUinner.jStart; j <= Fluid.gridUinner.jEnd; j++)
		{
			index = (i - Fluid.gridUinner.iStart) + (j - Fluid.gridUinner.jStart)*Fluid.gridUinner.iRes;
			U(i, j) = Fluid.tempU(index);
		}
	}
#pragma omp parallel for private(index)
	for (int i = Fluid.gridVinner.iStart; i <= Fluid.gridVinner.iEnd; i++)
	{
		for (int j = Fluid.gridVinner.jStart; j <= Fluid.gridVinner.jEnd; j++)
		{
			index = (i - Fluid.gridVinner.iStart) + (j - Fluid.gridVinner.jStart)*Fluid.gridVinner.iRes;
			V(i, j) = Fluid.tempV(index);
		}
	}

	// Boundary : Linear extension. (But, Dirichlet로 줄 방법은 없나???)
	Fluid.BdryCondVel();

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
				- dotProduct(unitNormal(i, j), gradient(i, j))*unitNormal(i, j);

			SurfaceForceX(i, j) = - meanCurvature(i, j)*SurfaceTension(i, j)*unitNormal(i, j).x + SurfGradSurfTension(i, j).x;
			SurfaceForceX(i, j) *= AdvectionMethod2D<double>::DeltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude() / (reynoldNum*Ca);
			SurfaceForceY(i, j) = - meanCurvature(i, j)*SurfaceTension(i, j)*unitNormal(i, j).y + SurfGradSurfTension(i, j).y;
			SurfaceForceY(i, j) *= -AdvectionMethod2D<double>::DeltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude() / (reynoldNum*Ca);

		}
		else
		{
			SurfGradSurfTension(i, j) = 0;

			SurfaceForceX(i, j) = 0;
			SurfaceForceY(i, j) = 0;
		}
	}

//#pragma omp parallel for private(i, j)
//	for (int k = 1; k <= levelSet.numTube; k++)
//	{
//		levelSet.TubeIndex(k, i, j);
//		if (levelSet.tube(i, j) == ComputedTubeRange)
//		{
//			SurfaceForceX(i, j) = meanCurvature(i, j)*SurfaceTension(i, j)*unitNormal(i, j).x - SurfGradSurfTension(i, j).x;
//			SurfaceForceX(i, j) *= -AdvectionMethod2D<double>::DeltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude() / (reynoldNum*Ca);
//			SurfaceForceY(i, j) = meanCurvature(i, j)*SurfaceTension(i, j)*unitNormal(i, j).y - SurfGradSurfTension(i, j).y;
//			SurfaceForceY(i, j) *= -AdvectionMethod2D<double>::DeltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude() / (reynoldNum*Ca);
//		}
//		else
//		{
//			SurfaceForceX(i, j) = 0;
//			SurfaceForceY(i, j) = 0;
//		}
//
//	}

}

inline void InsolubleSurfactant::GenerateLinearSystemUV(VectorND<double>& vectorB, const FD & vel, const FD & gradP, const FD & advec, const FD& force,const Grid2D & ipGrid, const double & scaling)
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

			//cout << endl;
			//cout << "(i,j) = (" << i << "," << j << ")" << endl;
			//cout << "index = " << index << endl;
			//cout << -gradP(i, j)<<"  "<< - advec(i, j) << endl;
			//cout << vel.dxxPhi(i, j) + vel.dyyPhi(i, j) << endl;
			//cout << "B " << vectorB(index) << endl;
			//cout << endl;

			vectorB(index) *= scaling;
		}
	}
}

inline void InsolubleSurfactant::PlotSurfactant()
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

inline void InsolubleSurfactant::PlotVelocity()
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
