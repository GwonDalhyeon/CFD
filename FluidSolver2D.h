#pragma once

#include "AdvectionMethod2D.h"
#include "CSR.h"

class FluidSolver2D
{
public:
	FluidSolver2D();
	~FluidSolver2D();

	Grid2D grid;
	Grid2D gridU;
	Grid2D gridV;

	//// Physical Variable
	FD Pressure; // Pressure

	FD U; // x velocity
	FD originU;

	FD V; // y velocity
	FD originV;

	FD Density; // density
	FD DensityU; // density for x velocity
	FD DensityV; // density for y velocity
	FD Mu; // Viscosity
	FD MuU; // Viscosity for x velocity
	FD MuV; // Viscosity for y velocity
	FD Nu;
	double reynoldNum;

	//////////////////////////////////////////////////////////
	////			Variable for						  ////
	////		two-phase flows with insoluble surfactant ////
	//////////////////////////////////////////////////////////
	double Ca;
	double Xi;
	double El;
	double Pe;

	//////////////////////////////////////////////////////////
	////			Variable for						  ////
	////		Coalescing Drop                           ////
	//////////////////////////////////////////////////////////
	double Oh;
	double We;


	bool isENOAdvection = true;;
	FD AdvectionU;
	FD AdvectionV;

	FD DiffusionU;
	FD DiffusionV;

	LS levelSet;

	Array2D<double>PMatrix;
	CSR<double> P_CSR;
	VectorND<double> Pb;
	VectorND<double> tempP;
	CSR<double> M;
	bool isPCG = false;
	//////////////////////////////////////////////////////////
	////			Numerical Variable              //////////
	//////////////////////////////////////////////////////////
	int ProjectionOrder = 1;
	double cflCondition;
	double dt;
	double finalT;
	double totalT;
	int maxIteration;
	int writeOutputIteration;
	int	iteration = 0;
	bool isPlot = true;


	inline void InitialCondition(const int& example);

	inline void Solver(const int & example);

	inline void TVDRK3TimeAdvection();
	inline void EulerMethod();
	inline void EulerMethodStep1();
	inline void EulerMethodStep2();
	inline void EulerMethodStep3();
	inline void GenerateLinearSystemPressure(CSR<double> & ipCSR);
	inline void GenerateLinearSystemPressure(VectorND<double>& vectorB);

	inline void AdvectionTerm(FD& U, FD& V, FD& TermU, FD& TermV);
	inline void DiffusionTerm(const FD& U, const FD& V, FD& TermU, FD& TermV);

	inline void TreatBCAlongXaxis(FD& ipField);
	inline void TreatBCAlongYaxis(FD& ipField);
	inline void TreatVelocityBC(FD& U, FD& V);
	inline void VectorToGrid(const VTN & ipVector, FD& ipField);

	inline void PlotVelocity();
	//////////////////////////////////////////////////////////////////////////
	//
	// "Accurate Projection Methods for the Incompressible Navier-Stokes Equations "
	//                 Brown, Cortez, Minion (2001)
	//     Sec 5.2 Projection Methods with a Lagged Pressure Term
	//
	/////////////////////////////////////////////////////////////////////////

	// Second Order Projection Method
	FD oldU;
	FD oldV;
	FD gradientPx;
	FD gradientPy;
	FD Phi;
	FD Phixxyy;
	FD AdvectionU1;
	FD AdvectionV1;
	FD AdvectionU2;
	FD AdvectionV2;

	// 'CN' means 'Crank-Nicolson'.
	Array2D<double> UCNMatrix;
	CSR<double> UCN_CSR;
	Array2D<int> Ubc;
	VectorND<double> Ub;
	VectorND<double> tempU;

	Array2D<double> VCNMatrix;
	CSR<double> VCN_CSR;
	Array2D<int> Vbc;
	VectorND<double> Vb;
	VectorND<double> tempV;

	Array2D<double> PhiCNMatrix;
	CSR<double> PhiCN_CSR;
	Array2D<int> Phibc;
	VectorND<double> Phib;
	VectorND<double> tempPhi;

	inline void EulerMethod2ndOrder();
	inline void EulerMethod2ndOrder1();
	inline void EulerMethod2ndOrder2();
	inline void EulerMethod2ndOrder3();
	inline void EulerMethod2ndOrder1stIteration1();

	inline void GenerateLinearSystemUV2Order(CSR<double> & ipU_CSR, CSR<double> & ipV_CSR);
	inline void GenerateLinearSystemUV2Order(VectorND<double>& vectorB);
	inline void GenerateLinearSystempPhi2Order(CSR<double> & ipCSR);
	inline void GenerateLinearSystempPhi2Order(VectorND<double>& vectorB);
private:

};

FluidSolver2D::FluidSolver2D()
{
}

FluidSolver2D::~FluidSolver2D()
{
}

inline void FluidSolver2D::InitialCondition(const int & example)
{
	if (example == 1)
	{
		cout << "*************************" << endl;
		cout << "    Navier-Stokes equation" << endl;
		cout << "    Chorin's Projection Method" << endl;
		cout << "    Cavity Flow" << endl;
		cout << "*************************" << endl;
		cout << endl;
		
		double domainLength = 1;
		int domainRes = 41;
		grid = Grid2D(-domainLength, domainLength, domainRes, -domainLength, domainLength, domainRes);
		gridU = Grid2D(grid.xMin - grid.dx / 2, grid.xMax + grid.dx / 2, grid.iRes + 1,
			grid.yMin, grid.yMax, grid.jRes);

		gridV = Grid2D(grid.xMin, grid.xMax, grid.iRes,
			grid.yMin - grid.dy / 2, grid.yMax + grid.dy / 2, grid.jRes + 1);

		Pressure = FD(grid);
		int tempBC = 0;
		//// Boundary Condition
		for (int j = Pressure.jStart; j <= Pressure.jEnd; j++)
		{
			for (int i = Pressure.iStart; i <= Pressure.iEnd; i++)
			{
				Pressure.BC(i, j) = tempBC++;
			}
		}

		U = FD(gridU);
		tempBC = 0;
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			for (int i = U.iStart; i <= U.iEnd; i++)
			{
				if (i == U.iStart || i == U.iEnd)
				{
					U.BC(i, j) = BC_REFLECTION;
					continue;
				}
				if (j == U.jStart)
				{
					U.BC(i, j) = BC_DIR;
					U(i, j) = 0;
					continue;
				}
				if (j == U.jEnd)
				{
					U.BC(i, j) = BC_DIR;
					U(i, j) = 1;
					continue;
				}
				U.BC(i, j) = tempBC++;
			}
		}
		TreatBCAlongYaxis(U);
		originU = U;

		V = FD(gridV);

		tempBC = 0;
		for (int j = V.jStart; j <= V.jEnd; j++)
		{
			for (int i = V.iStart; i <= V.iEnd; i++)
			{
				if (i == V.iStart || i == V.iEnd)
				{
					V.BC(i, j) = BC_DIR;
					continue;
				}
				if (j == V.jStart || j == V.jEnd)
				{
					V.BC(i, j) = BC_REFLECTION;
					continue;
				}
				V.BC(i, j) = tempBC++;
			}
		}
		originV = V;

		Density = FD(grid);
		DensityU = FD(gridU);
		DensityV = FD(gridV);
		Mu = FD(grid);
		MuU = FD(gridU);
		MuV = FD(gridV);
		Nu = FD(grid);

		reynoldNum = 1000;

		AdvectionU = FD(gridU);
		AdvectionV = FD(gridV);

		DiffusionU = FD(gridU);
		DiffusionV = FD(gridV);

		levelSet = LS(grid);

		Phi = FD(grid);
		Phixxyy = FD(grid);
		gradientPx = FD(U.innerGrid);
		gradientPy = FD(V.innerGrid);

		//// Projection Accuracy Order
		ProjectionOrder = 1; 
		//isPCG = false;
		if (ProjectionOrder == 1)
		{
			GenerateLinearSystemPressure(P_CSR);
		}
		else if (ProjectionOrder == 2)
		{
			GenerateLinearSystemUV2Order(UCN_CSR, VCN_CSR);
			GenerateLinearSystempPhi2Order(PhiCN_CSR);
		}

		cflCondition = 0.4;
		dt = cflCondition*grid.dx;
		finalT = 30;
		maxIteration = ceil(finalT / dt);
		totalT = 0;
		writeOutputIteration = 30;
		iteration = 0;
	}

	if (example == 2)
	{
		cout << "*************************" << endl;
		cout << "    Navier-Stokes equation" << endl;
		cout << "    Chorin's Projection Method" << endl;
		cout << "    Tube" << endl;
		cout << "*************************" << endl;
		cout << endl;

		grid = Grid2D(-0.5, 0.5, 41, -0.5, 0.5, 41);
		gridU = Grid2D(grid.xMin - grid.dx / 2, grid.xMax + grid.dx / 2, grid.iRes + 1,
			grid.yMin, grid.yMax, grid.jRes);

		gridV = Grid2D(grid.xMin, grid.xMax, grid.iRes,
			grid.yMin - grid.dy / 2, grid.yMax + grid.dy / 2, grid.jRes + 1);

		Pressure = FD(grid);
		int tempBC = 0;
		//// Boundary Condition
		for (int j = Pressure.jStart; j <= Pressure.jEnd; j++)
		{
			for (int i = Pressure.iStart; i <= Pressure.iEnd; i++)
			{
				Pressure.BC(i, j) = tempBC++;
			}
		}

		U = FD(gridU);

#pragma omp parallel for
		for (int i = gridU.iStart; i <= gridU.iEnd; i++)
		{
			U(i, U.jStart) = 0;
			U.BC(i, U.jStart) = BC_DIR;
			U(i, U.jEnd) = 0;
			U.BC(i, U.jEnd) = BC_DIR;

		}
		tempBC = 0;
		for (int j = U.jStartI; j <= U.jEndI; j++)
		{
			for (int i = U.iStart; i <= U.iEnd; i++)
			{
				U(i, j) = 1;
				if (i == U.iStart)
				{
					U.BC(i, j) = BC_DIR;
					continue;
				}
				if (i == U.iEnd)
				{
					U.BC(i, j) = BC_NEUM;
					continue;
				}
				U.BC(i, j) = tempBC++;
			}
		}
		originU = U;

		V = FD(gridV);

		tempBC = 0;
		for (int j = V.jStart; j <= V.jEnd; j++)
		{
			for (int i = V.iStart; i <= V.iEnd; i++)
			{
				if (i == V.iStart || i == V.iEnd)
				{
					V.BC(i, j) = BC_DIR;
					continue;
				}
				if (j == V.jStart || j == V.jEnd)
				{
					V.BC(i, j) = BC_REFLECTION;
					continue;
				}
				V.BC(i, j) = tempBC++;
			}
		}
		originV = V;

		Density = FD(grid);
		DensityU = FD(gridU);
		DensityV = FD(gridV);
		Mu = FD(grid);
		MuU = FD(gridU);
		MuV = FD(gridV);
		Nu = FD(grid);

		reynoldNum = 1000;

		AdvectionU = FD(gridU);
		AdvectionV = FD(gridV);

		DiffusionU = FD(gridU);
		DiffusionV = FD(gridV);

		levelSet = LS(grid);

		Phi = FD(grid);
		Phixxyy = FD(grid);
		gradientPx = FD(U.innerGrid);
		gradientPy = FD(V.innerGrid);

		//// Projection Accuracy Order
		ProjectionOrder = 1;
		if (ProjectionOrder == 1)
		{
			GenerateLinearSystemPressure(P_CSR);
		}
		else if (ProjectionOrder == 2)
		{
			GenerateLinearSystemUV2Order(UCN_CSR, VCN_CSR);
			GenerateLinearSystempPhi2Order(PhiCN_CSR);
		}

		cflCondition = 0.1;
		dt = cflCondition*grid.dx;
		finalT = 10;
		maxIteration = int(finalT / dt);
		totalT = 0;
		writeOutputIteration = 10;
		iteration = 0;
	}

	if (example == 3)
	{
		/////////////////////////////////////////////////////////
		/////  A level-set continuum method
		/////  for two-phase flows with insoluble surfactant
		/////  --JJ Xu, Y Yang, J Lowengrub--
		/////  Example 1
		/////////////////////////////////////////////////////////
		int numP = grid.iRes;
		double ddx = grid.dx;

		gridU = Grid2D(grid.xMin - grid.dx / 2, grid.xMax + grid.dx / 2, grid.iRes + 1,
			grid.yMin, grid.yMax, grid.jRes);
		gridV = Grid2D(grid.xMin, grid.xMax, grid.iRes,
			grid.yMin - grid.dy / 2, grid.yMax + grid.dy / 2, grid.jRes + 1);

		Pressure = FD(grid);
		int tempBC = 0;
		//// Boundary Condition
		for (int j = Pressure.jStart; j <= Pressure.jEnd; j++)
		{
			for (int i = Pressure.iStart; i <= Pressure.iEnd; i++)
			{
				if (i == Pressure.iStart && j == Pressure.jStartI)
				{
					Pressure.BC(i, j) = BC_DIR;
					Pressure(i, j) = 1;
					continue;
				}
				if (i == Pressure.iStart || i == Pressure.iEnd)
				{
					Pressure.BC(i, j) = BC_NEUM;
					continue;
				}
				if (j == Pressure.jStart || j == Pressure.jEnd)
				{
					Pressure.BC(i, j) = BC_NEUM;
					continue;
				}
				Pressure.BC(i, j) = tempBC++;;
			}
		}

		U = FD(gridU);
		tempBC = 0;
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			for (int i = U.iStart; i <= U.iEnd; i++)
			{
				U(i, j) = gridU(i, j).y;

				if (i == U.iStart || i == U.iEnd)
				{
					U.BC(i, j) = BC_DIR;

					continue;
				}
				if (j == U.jStart || j == U.jEnd)
				{
					U.BC(i, j) = BC_DIR;
					continue;
				}
				U.BC(i, j) = tempBC++;
			}
		}
		originU = U;

		V = FD(gridV);
		tempBC = 0;
		for (int j = V.jStart; j <= V.jEnd; j++)
		{
			for (int i = V.iStart; i <= V.iEnd; i++)
			{
				if (i == V.iStart || i == V.iEnd)
				{
					V.BC(i, j) = BC_NEUM;
					continue;
				}
				if (j == V.jStart || j == V.jEnd)
				{
					V.BC(i, j) = BC_NEUM;
					continue;
				}
				V.BC(i, j) = tempBC++;
			}
		}
		originV = V;

		Density = FD(grid);
		DensityU = FD(gridU);
		DensityV = FD(gridV);
		Mu = FD(grid);
		MuU = FD(gridU);
		MuV = FD(gridV);
		Nu = FD(grid);

		reynoldNum = 1000;

		AdvectionU = FD(gridU);
		AdvectionV = FD(gridV);

		DiffusionU = FD(gridU);
		DiffusionV = FD(gridV);

		levelSet = LS(grid);

		Phi = FD(grid);
		Phixxyy = FD(grid);
		gradientPx = FD(U.innerGrid);
		gradientPy = FD(V.innerGrid);

		//// Projection Accuracy Order
		ProjectionOrder = 1;
		if (ProjectionOrder == 1)
		{
			GenerateLinearSystemPressure(P_CSR);
			//P_CSR.indPrt.Variable("PindPrt2");
			//P_CSR.values.Variable("Pvalues2");
			//P_CSR.columns.Variable("Pcolumns2");

		}
		else if (ProjectionOrder == 2)
		{
			GenerateLinearSystemUV2Order(UCN_CSR, VCN_CSR);
			GenerateLinearSystempPhi2Order(PhiCN_CSR);
		}

		cflCondition = 0.1;
		dt = cflCondition*grid.dx;
		finalT = 10;
		maxIteration = int(finalT / dt);
		totalT = 0;
		writeOutputIteration = 10;
		iteration = 0;

		isPlot = false;
	}

	if (example == 5)
	{
		/////////////////////////////////////////////////////////
		/////  Simulations of surfactant effects
		/////  on the dynamics of coalescing drops and bubbles
		/////  --DW Martin and F Blanchette--
		/////  Example 1
		/////////////////////////////////////////////////////////
		int numP = grid.iRes;
		double ddx = grid.dx;

		gridU = Grid2D(grid.xMin - grid.dx / 2, grid.xMax + grid.dx / 2, grid.iRes + 1,
			grid.yMin, grid.yMax, grid.jRes);
		gridV = Grid2D(grid.xMin, grid.xMax, grid.iRes,
			grid.yMin - grid.dy / 2, grid.yMax + grid.dy / 2, grid.jRes + 1);

		Pressure = FD(grid);
		int tempBC = 0;
		for (int j = Pressure.jStart; j <= Pressure.jEnd; j++)
		{
			for (int i = Pressure.iStart; i <= Pressure.iEnd; i++)
			{
				if (i == Pressure.iStart && j == Pressure.jStartI)
				{
					Pressure.BC(i, j) = BC_DIR;
					Pressure(i, j) = 1;
					continue;
				}
				if (i == Pressure.iStart || i == Pressure.iEnd || j == Pressure.jStart || j == Pressure.jEnd)
				{
					Pressure.BC(i, j) = BC_NEUM;
					continue;
				}
				Pressure.BC(i, j) = tempBC++;
			}
		}

		U = FD(gridU);
		tempBC = 0;
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			for (int i = U.iStart; i <= U.iEnd; i++)
			{
				if (i == U.iStart || i == U.iEnd || j == U.jStart || j == U.jEnd)
				{
					U.BC(i, j) = BC_NEUM;
					continue;
				}
				U.BC(i, j) = tempBC++;
			}
		}
		originU = U;

		V = FD(gridV);
		tempBC = 0;
		for (int j = V.jStart; j <= V.jEnd; j++)
		{
			for (int i = V.iStart; i <= V.iEnd; i++)
			{
				if (i == V.iStart || i == V.iEnd || j == V.jStart || j == V.jEnd)
				{
					V.BC(i, j) = BC_NEUM;
					continue;
				}
				V.BC(i, j) = tempBC++;
			}
		}
		originV = V;



		Density = FD(grid);
		DensityU = FD(gridU);
		DensityV = FD(gridV);
		Mu = FD(grid);
		MuU = FD(gridU);
		MuV = FD(gridV);
		Nu = FD(grid);

		reynoldNum = 1000;

		AdvectionU = FD(gridU);
		AdvectionV = FD(gridV);

		DiffusionU = FD(gridU);
		DiffusionV = FD(gridV);

		levelSet = LS(grid);

		Phi = FD(grid);
		Phixxyy = FD(grid);
		gradientPx = FD(U.innerGrid);
		gradientPy = FD(V.innerGrid);

		isPlot = false;
		ProjectionOrder = 1;

		if (ProjectionOrder == 1)
		{
			//GenerateLinearSystem(PMatrix, -grid.dx*grid.dx);
			GenerateLinearSystemPressure(P_CSR);
		}
		else
		{
			GenerateLinearSystemUV2Order(UCN_CSR, VCN_CSR);
			GenerateLinearSystempPhi2Order(PhiCN_CSR);

			//GenerateLinearSystemUV(Fluid.UCNMatrix, Fluid.U.innerGrid, 1, Fluid.UCN_CSR);
			//GenerateLinearSystemUV(Fluid.VCNMatrix, Fluid.V.innerGrid, 1, Fluid.VCN_CSR);

		}

	}
}

inline void FluidSolver2D::Solver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;

	InitialCondition(example);
	grid.Variable("Xp", "Yp");
	gridU.Variable("Xu", "Yu");
	gridV.Variable("Xv", "Yv");

	//isPlot = false;
	if (isPlot)
	{
		MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
		PlotVelocity();
		MATLAB.WriteImage("fluid", iteration, "png");
	}

	for (iteration = 1; iteration <= maxIteration; iteration++)
	{

		cout << endl;
		cout << "********************************" << endl;
		cout << "       Iteration " << to_string(iteration) << " : Start" << endl;
		TVDRK3TimeAdvection();
		totalT += dt;
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;
		cout << "********************************" << endl;

		if (isPlot)
		{
			if (iteration == 1 || iteration % 1 == 0)
			{
				//MATLAB.Command("subplot(1,2,1)");
				//PlotVelocity();
				//MATLAB.Command("subplot(1,2,2)");
				//Pressure.Variable("P");
				//MATLAB.Command("surf(P)");
				//MATLAB.Command("divU =U(:,2:end)-U(:,1:end-1),divV =V(2:end,:)-V(1:end-1,:);div=divU+divV;");
				if ((iteration == 1 || iteration % 10 == 0) && true)
				{
					MATLAB.Command("subplot(1,2,1)");
					PlotVelocity();
					MATLAB.Command("subplot(1,2,2)");
					Pressure.Variable("P");
					MATLAB.Command("surf(P)");
					//MATLAB.WriteImage("fluid", iteration, "fig");
					MATLAB.WriteImage("fluid", iteration, "png");
				}
			}
		}


		if (writeFile && iteration%writeOutputIteration == 0)
		{
			fileName = "pressure" + to_string(iteration);
			Pressure.WriteFile(fileName);
			fileName = "xVelocity" + to_string(iteration);
			U.WriteFile(fileName);
			fileName = "yVelocity" + to_string(iteration);
			V.WriteFile(fileName);
		}
	}


}

inline void FluidSolver2D::TVDRK3TimeAdvection()
{
	originU.dataArray = U.dataArray;
	originV.dataArray = V.dataArray;


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

inline void FluidSolver2D::EulerMethod()
{
	////////////////////////////////////////////////
	////     Projection Method 1 : advection    ////
	////////////////////////////////////////////////
	EulerMethodStep1();
	
	TreatBCAlongYaxis(U);
	TreatBCAlongXaxis(V);

	//U.Variable("Ustar");
	//V.Variable("Vstar");
	//MATLAB.Command("divUstar =Ustar(:,2:end)-Ustar(:,1:end-1),divVstar =Vstar(2:end,:)-Vstar(1:end-1,:);divstar=divUstar+divVstar;");
	//MATLAB.Command("quiver(Xp,Yp,Ustar(:,1:end-1),Vstar(1:end-1,:))");
	
	////////////////////////////////////////////////
	////     Projection Method 2 : Poisson Eq   ////
	////////////////////////////////////////////////
	EulerMethodStep2();


	//////////////////////////////////////////////
	////     Projection Method 3 : New U,V    ////
	//////////////////////////////////////////////
	EulerMethodStep3();
	
	//U.Variable("Unew");
	//V.Variable("Vnew");
	//MATLAB.Command("divUnew =Unew(:,2:end)-Unew(:,1:end-1),divVnew =Vnew(2:end,:)-Vnew(1:end-1,:);divnew=divUnew+divVnew;");
	//MATLAB.Command("quiver(Xp,Yp,Unew(1:end-1,:),Vnew(:,1:end-1))");

	TreatBCAlongYaxis(U);
	TreatBCAlongXaxis(V);
}

inline void FluidSolver2D::EulerMethodStep1()
{
	Array2D<double>& K1U = U.K1;
	Array2D<double>& K1V = V.K1;

	AdvectionTerm(U, V, AdvectionU, AdvectionV);
	DiffusionTerm(U, V, DiffusionU, DiffusionV);
	double oneOverRe = 1 / reynoldNum;
	//AdvectionU.Variable("advectionU");
	//DiffusionU.Variable("diffusionU");
	//AdvectionV.Variable("advectionV");
	//DiffusionV.Variable("diffusionV");

#pragma omp parallel for
	for (int i = U.iStart; i <= U.iEnd; i++)
	{
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			if (U.BC(i, j) < 0) continue;

			K1U(i, j) = dt*(-AdvectionU(i, j) + oneOverRe*DiffusionU(i, j));
			U(i, j) += K1U(i, j);
		}
	}
#pragma omp parallel for
	for (int i = V.iStart; i <= V.iEnd; i++)
	{
		for (int j = V.jStart; j <= V.jEnd; j++)
		{
			if (V.BC(i, j) < 0) continue;

			K1V(i, j) = dt*(-AdvectionV(i, j) + oneOverRe*DiffusionV(i, j));
			V(i, j) += K1V(i, j);
		}
	}
	//K1U.Variable("k1u");
	//K1V.Variable("k1v");
}

inline void FluidSolver2D::EulerMethodStep2()
{
	GenerateLinearSystemPressure(Pb);
	Array2D<int> BC = Pressure.BC;
#pragma omp parallel for
	for (int i = Pressure.iStart; i <= Pressure.iEnd; i++)
	{
		for (int j = Pressure.jStart; j <= Pressure.jEnd; j++)
		{
			if (BC(i, j)<0) continue;

			tempP(BC(i, j)) = Pressure(i, j);
		}
	}

	if (isPCG)
	{
		PCGSolver::Solver(P_CSR, M, Pb, tempP, BC);
	}
	else
	{
		CGSolver::Solver(P_CSR, Pb, tempP);
	}

	VectorToGrid(tempP, Pressure);
}

inline void FluidSolver2D::EulerMethodStep3()
{
	double oneOverdx = Pressure.oneOverdx;
#pragma omp parallel for
	for (int i = U.iStart; i <= U.iEnd; i++)
	{
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			if (U.BC(i, j) < 0) continue;

			U(i, j) -= dt * (Pressure(i, j) - Pressure(i - 1, j))*oneOverdx;
		}
	}

	double oneOverdy = Pressure.oneOverdy;
#pragma omp parallel for
	for (int i = V.iStart; i <= V.iEnd; i++)
	{
		for (int j = V.jStart; j <= V.jEnd; j++)
		{
			if (V.BC(i, j) < 0) continue;

			V(i, j) -= dt * (Pressure(i, j) - Pressure(i, j - 1))*oneOverdy;
		}
	}
}


inline void FluidSolver2D::GenerateLinearSystemPressure(CSR<double>& ipCSR)
{
	Pressure.CountNonZero();
	Pb = VectorND<double>(Pressure.num_all_full_cells);
	tempP = VectorND<double>(Pressure.num_all_full_cells);
	P_CSR = CSR<double>(Pressure.num_all_full_cells, Pressure.nnz);

	int iStart = Pressure.iStart, iEnd = Pressure.iEnd, jStart = Pressure.jStart, jEnd = Pressure.jEnd;

	double dx = Pressure.dx, dy = Pressure.dy;
	double dx2 = dx*dx, dy2 = dy*dy, dxdy = dx*dy;
	double oneOverdx = 1 / dx, oneOverdx2 = 1 / dx2;
	double oneOverdy = 1 / dy, oneOverdy2 = 1 / dy2;

	Array2D<int>& BC = Pressure.BC;

	double coefIJ;
	for (int j = jStart; j <= jEnd; j++)
	{
		for (int i = iStart; i <= iEnd; i++)
		{
			coefIJ = 0;

			if (BC(i, j) < 0) continue;

			//// If neighbor is full cell
			if (i>iStart)
			{
				if (BC(i - 1, j) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i - 1, j), -oneOverdx2);
				}
			}
			if (i<iEnd)
			{
				if (BC(i + 1, j) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i + 1, j), -oneOverdx2);
				}
			}
			if (j>jStart)
			{
				if (BC(i, j - 1) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i, j - 1), -oneOverdy2);
				}
			}
			if (j<jEnd)
			{
				if (BC(i, j + 1) > -1)
				{
					coefIJ++;
					ipCSR.AssignValue(BC(i, j), BC(i, j + 1), -oneOverdy2);
				}
			}

			////// Dirichlet Boundary Condition
			//if (i>iStart)	if (BC(i - 1, j) == BC_DIR)	coefIJ++;
			//if (i<iEnd)		if (BC(i + 1, j) == BC_DIR)	coefIJ++;
			//if (j>jStart)	if (BC(i, j - 1) == BC_DIR)	coefIJ++;
			//if (j<jEnd)		if (BC(i, j + 1) == BC_DIR)	coefIJ++;

			//if (i == Pressure.iStartI && j == Pressure.jStartI)
			//{
			//	if (BC(i - 1, j) == BC_NEUM) coefIJ++;
			//	if (BC(i + 1, j) == BC_NEUM) coefIJ++;
			//	if (BC(i, j - 1) == BC_NEUM) coefIJ++;
			//	if (BC(i, j + 1) == BC_NEUM) coefIJ++;
			//}
			//else
			//{
			//	if (i>iStart)	if (BC(i - 1, j) == BC_NEUM) coefIJ += 0;
			//	if (i<iEnd)		if (BC(i + 1, j) == BC_NEUM) coefIJ += 0;
			//	if (j>jStart)	if (BC(i, j - 1) == BC_NEUM) coefIJ += 0;
			//	if (j<jEnd)		if (BC(i, j + 1) == BC_NEUM) coefIJ += 0;
			//}

			//if (coefIJ == 0) coefIJ = 1;
			
			ipCSR.AssignValue(BC(i, j), BC(i, j), oneOverdx2*coefIJ);
		}
	}

	if (isPCG)
	{
		int N = Pressure.num_all_full_cells;
		int nz = P_CSR.valueNum;

		M = CSR<double>(Pressure.num_all_full_cells, (nz - N) / 2 + N);
		PCGSolver::IncompleteCholeskyDecomposition(P_CSR, M, Pressure.BC);
	}
}

inline void FluidSolver2D::GenerateLinearSystemPressure(VectorND<double>& vectorB)
{
	double tempSum = 0;
	Array2D<int>& BC = Pressure.BC;
	double oneOverdt = 1. / dt;
	double oneOverdx = U.oneOverdx;
	double oneOverdy = U.oneOverdy;
	double oneOverdx2 = Pressure.oneOverdx2;
	double oneOverdy2 = Pressure.oneOverdy2;
//#pragma omp parallel for reduction(+ : tempSum)
	for (int j = Pressure.jStart; j <= Pressure.jEnd; j++)
	{
		for (int i = Pressure.iStart; i <= Pressure.iEnd; i++)
		{
			if (BC(i, j) < 0) continue;

			vectorB(BC(i, j)) = -oneOverdt * ((U(i + 1, j) - U(i, j))*oneOverdx + (V(i, j + 1) - V(i, j))*oneOverdy);
			
			if (i > Pressure.iStart)	if (BC(i - 1, j) == BC_DIR) vectorB(BC(i, j)) += Pressure(i - 1, j)*oneOverdx2;
			if (i < Pressure.iEnd)		if (BC(i + 1, j) == BC_DIR) vectorB(BC(i, j)) += Pressure(i + 1, j)*oneOverdx2;
			if (j > Pressure.jStart)	if (BC(i, j - 1) == BC_DIR) vectorB(BC(i, j)) += Pressure(i, j - 1)*oneOverdy2;
			if (j < Pressure.jEnd)		if (BC(i, j + 1) == BC_DIR) vectorB(BC(i, j)) += Pressure(i, j + 1)*oneOverdy2;

			tempSum += vectorB(BC(i, j));
		}
	}
	
	double ave = tempSum / (double) vectorB.iLength;
#pragma omp parallel for
	for (int i = 0; i < vectorB.iLength; i++)
	{
		vectorB(i) -= ave;
	}
}

inline void FluidSolver2D::AdvectionTerm(FD & U, FD & V, FD & TermU, FD & TermV)
{
	Array2D<double>& dUdxM = U.dfdxM;
	Array2D<double>& dUdxP = U.dfdxP;
	Array2D<double>& dUdyM = U.dfdyM;
	Array2D<double>& dUdyP = U.dfdyP;

	Array2D<double>& dVdxM = V.dfdxM;
	Array2D<double>& dVdxP = V.dfdxP;
	Array2D<double>& dVdyM = V.dfdyM;
	Array2D<double>& dVdyP = V.dfdyP;
	//isENOAdvection = false;
	if (isENOAdvection)
	{
		AdvectionMethod2D<double>::ENO3rdDerivation(U, dUdxM, dUdxP, dUdyM, dUdyP);
		AdvectionMethod2D<double>::ENO3rdDerivation(V, dVdxM, dVdxP, dVdyM, dVdyP);
	}
	else
	{
		AdvectionMethod2D<double>::WENO3rdDerivation(U, dUdxM, dUdxP, dUdyM, dUdyP);
		AdvectionMethod2D<double>::WENO3rdDerivation(V, dVdxM, dVdxP, dVdyM, dVdyP);
	}
	
	double Ux, Uy;
	double aveV;
#pragma omp parallel for private(aveV, Ux, Uy)
	for (int i = TermU.iStart; i <= TermU.iEnd; i++)
	{
		for (int j = TermU.jStart; j <= TermU.jEnd; j++)
		{
			if (U.BC(i, j) < 0) continue;

			aveV = (V(i - 1, j) + V(i, j) + V(i - 1, j + 1) + V(i, j + 1)) / 4.;
			
			if (U(i, j) > 0)
			{
				Ux = dUdxM(i, j);
			}
			else
			{
				Ux = dUdxP(i, j);
			}

			if (aveV > 0)
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
			if (V.BC(i, j) < 0) continue;

			aveU = (U(i, j - 1) + U(i + 1, j - 1) + U(i, j) + U(i + 1, j)) / 4.;
			if (aveU > 0)
			{
				Vx = dVdxM(i, j);
			}
			else
			{
				Vx = dVdxP(i, j);
			}

			if (V(i, j) > 0)
			{
				Vy = dVdyM(i, j);
			}
			else
			{
				Vy = dVdyP(i, j);
			}
			TermV(i, j) = aveU*Vx + V(i, j)*Vy;
		}
	}
}

inline void FluidSolver2D::DiffusionTerm(const FD & U, const FD & V, FD & TermU, FD & TermV)
{
#pragma omp parallel for
	for (int i = U.iStart; i <= U.iEnd; i++)
	{
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			if (U.BC(i, j) < 0) continue;

			TermU(i, j) = 2 * ((U(i + 1, j) - U(i, j))*U.oneOverdx - (U(i, j) - U(i - 1, j))*U.oneOverdx)*U.oneOverdx
				+ ((U(i, j + 1) - U(i, j))*U.oneOverdy + (V(i, j + 1) - V(i - 1, j + 1))*V.oneOverdx
					- ((U(i, j) - U(i, j - 1))*U.oneOverdy + (V(i, j) - V(i - 1, j))*V.oneOverdx))*U.oneOverdy;
		}
	}

#pragma omp parallel for 
	for (int i = TermV.iStart; i <= TermV.iEnd; i++)
	{
		for (int j = TermV.jStart; j <= TermV.jEnd; j++)
		{
			if (V.BC(i, j) < 0) continue;

			TermV(i, j) = 2 * ((V(i, j + 1) - V(i, j))*V.oneOverdy - (V(i, j) - V(i, j - 1))*V.oneOverdy)*V.oneOverdy
				+ ((U(i + 1, j) - U(i + 1, j - 1))*U.oneOverdy + (V(i + 1, j) - V(i, j))*V.oneOverdx
					- ((U(i, j) - U(i, j - 1))*U.oneOverdy + (V(i, j) - V(i - 1, j))*V.oneOverdx))*V.oneOverdx;
		}
	}
}

inline void FluidSolver2D::TreatBCAlongXaxis(FD & ipField)
{
#pragma omp parallel for 
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		if (ipField.BC(i, ipField.jStart) == BC_NEUM)		ipField(i, ipField.jStart) = ipField(i, ipField.jStartI);
		if (ipField.BC(i, ipField.jStart) == BC_REFLECTION) ipField(i, ipField.jStart) = -ipField(i, ipField.jStartI);
		if (ipField.BC(i, ipField.jEnd) == BC_NEUM)			ipField(i, ipField.jEnd)   = ipField(i, ipField.jEndI);
		if (ipField.BC(i, ipField.jEnd) == BC_REFLECTION)	ipField(i, ipField.jEnd)   = -ipField(i, ipField.jEndI);
	}
}

inline void FluidSolver2D::TreatBCAlongYaxis(FD & ipField)
{
#pragma omp parallel for 
	for (int j = ipField.jStart; j <= ipField.jEnd; j++)
	{
		if (ipField.BC(ipField.iStart, j) == BC_NEUM)		ipField(ipField.iStart, j) = ipField(ipField.iStartI, j);
		if (ipField.BC(ipField.iStart, j) == BC_REFLECTION) ipField(ipField.iStart, j) = -ipField(ipField.iStartI, j);
		if (ipField.BC(ipField.iEnd, j) == BC_NEUM)			ipField(ipField.iEnd, j)   = ipField(ipField.iEndI, j);
		if (ipField.BC(ipField.iEnd, j) == BC_REFLECTION)	ipField(ipField.iEnd, j)   = -ipField(ipField.iEndI, j);
	}
}


inline void FluidSolver2D::TreatVelocityBC(FD & U, FD & V)
{
#pragma omp parallel for 
	for (int i = U.iStart; i <= U.iEnd; i++)
	{
		if (U.BC(i, U.jStart) == BC_NEUM)		U(i, U.jStart) =  U(i, U.jStartI);
		if (U.BC(i, U.jEnd)   == BC_NEUM)		U(i, U.jEnd)   =  U(i, U.jEndI);
		if (U.BC(i, U.jStart) == BC_REFLECTION) U(i, U.jStart) = -U(i, U.jStartI);
		if (U.BC(i, U.jEnd)   == BC_REFLECTION)	U(i, U.jEnd)   = -U(i, U.jEndI);
	}
#pragma omp parallel for 
	for (int j = U.jStart; j <= U.jEnd; j++)
	{
		if (U.BC(U.iStart, j) == BC_NEUM)		U(U.iStart, j) =  U(U.iStartI, j); 
		if (U.BC(U.iEnd, j)   == BC_NEUM)		U(U.iEnd, j)   =  U(U.iEndI,   j);
		if (U.BC(U.iStart, j) == BC_REFLECTION) U(U.iStart, j) = -U(U.iStartI, j);
		if (U.BC(U.iEnd, j)   == BC_REFLECTION) U(U.iEnd, j)   = -U(U.iEndI,   j);
	}

#pragma omp parallel for 
	for (int j = V.jStart; j <= V.jEnd; j++)
	{
		if (V.BC(V.iStart, j) == BC_NEUM)		V(V.iStart, j) =  V(V.iStartI, j);
		if (V.BC(V.iEnd, j)   == BC_NEUM)		V(V.iEnd, j)   =  V(V.iEndI,   j);
		if (V.BC(V.iStart, j) == BC_REFLECTION) V(V.iStart, j) = -V(V.iStartI, j);
		if (V.BC(V.iEnd, j)   == BC_REFLECTION) V(V.iEnd, j)   = -V(V.iEndI,   j);
	}
#pragma omp parallel for 
	for (int i = V.iStart; i <= V.iEnd; i++)
	{
		if (V.BC(i, V.jStart) == BC_NEUM)		V(i, V.jStart) =  V(i, V.jStartI);
		if (V.BC(i, V.jEnd)   == BC_NEUM)		V(i, V.jEnd)   =  V(i, V.jEndI);
		if (V.BC(i, V.jStart) == BC_REFLECTION) V(i, V.jStart) = -V(i, V.jStartI);
		if (V.BC(i, V.jEnd)   == BC_REFLECTION) V(i, V.jEnd)   = -V(i, V.jEndI);
	}
}

inline void FluidSolver2D::VectorToGrid(const VTN & ipVector, FD & ipField)
{
#pragma omp parallel for 
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
		{
			if (ipField.BC(i, j)<0) continue;

			ipField(i, j) = ipVector(ipField.BC(i, j));
		}
	}
}

inline void FluidSolver2D::PlotVelocity()
{
	string str;

	U.Variable("U");
	V.Variable("V");
	
	str = string("quiver(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),axis([Xp(1)-(Xp(end)-Xp(1))/10 Xp(end)+(Xp(end)-Xp(1))/10 Yp(1)-(Yp(end)-Yp(1))/10 Yp(end)+(Yp(end)-Yp(1))/10]);");
	//str = str + string("hold on,streamline(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,Xp(1:30:end),Yp(1:30:end)), hold off;axis equal; axis([Xp(1) Xp(end) Yp(1) Yp(end)]);");
	str = str + string("hold on,streamslice(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2), hold off;axis equal; axis([Xp(1) Xp(end) Yp(1) Yp(end)]);");
	MATLAB.Command(str.c_str());
	str = string("title(['iteration : ', num2str(") + to_string(iteration) + string("),', time : ', num2str(") + to_string(totalT) + string(")]);");
	MATLAB.Command(str.c_str());
}

inline void FluidSolver2D::EulerMethod2ndOrder()
{
}

inline void FluidSolver2D::EulerMethod2ndOrder1()
{
}

inline void FluidSolver2D::EulerMethod2ndOrder2()
{
}

inline void FluidSolver2D::EulerMethod2ndOrder3()
{
}

inline void FluidSolver2D::EulerMethod2ndOrder1stIteration1()
{
}

inline void FluidSolver2D::GenerateLinearSystemUV2Order(CSR<double>& ipU_CSR, CSR<double>& ipV_CSR)
{
	U.CountNonZero();
	Ub = VectorND<double>(U.num_all_full_cells);
	tempU = VectorND<double>(U.num_all_full_cells);

	V.CountNonZero();
	Vb = VectorND<double>(V.num_all_full_cells);
	tempV = VectorND<double>(V.num_all_full_cells);
}

inline void FluidSolver2D::GenerateLinearSystemUV2Order(VectorND<double>& vectorB)
{

}

inline void FluidSolver2D::GenerateLinearSystempPhi2Order(CSR<double>& ipCSR)
{
	Phi.CountNonZero();
	Phib = VectorND<double>(Phi.num_all_full_cells);
	tempPhi = VectorND<double>(Phi.num_all_full_cells);
}

inline void FluidSolver2D::GenerateLinearSystempPhi2Order(VectorND<double>& vectorB)
{

}