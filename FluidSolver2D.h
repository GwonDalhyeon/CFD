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
	FD Viscosity; // Viscosity
	FD ViscosityU; // Viscosity for x velocity
	FD ViscosityV; // Viscosity for y velocity
	FD Nu;
	double Re = 1000;


	//////////////////////////////////////////////////////////
	////			Variable for						  ////
	////		two-phase flows                           ////
	////    with boundary condition capturing method      ////
	//////////////////////////////////////////////////////////
	bool dimensionlessForm = true;
	bool isMultiPhase = false;
	bool isGravity = false;
	bool isCSFmodel = true;
	const double gravity = -9.8; // -9.8 m/s^2
	// Water : 1000kg/m^3, Air : 1.226kg/m^3
	const double densityWater = 1000;
	const double densityAir = 1.226;
	double densityI = 1; 
	double densityE = 1;
	// Water : 1.137*E-3kg/ms, Air : 1.78*E-5kg/ms
	const double viscosityWater = 1.137*pow(10, -3);
	const double viscosityAir = 1.78*pow(10, -5);
	double viscosityI = 1; //
	double viscosityE = 1;
	double gamma0 = 1; // 0.7825kg/s^2

	
	FD SurfaceForce;
	FD SurfaceForceX;
	FD SurfaceForceY;

	//////////////////////////////////////////////////////////
	////			Variable for						  ////
	////		two-phase flows with insoluble surfactant ////
	//////////////////////////////////////////////////////////
	double Ca = 1;
	double Xi = 0.3;
	double El = 1;
	double Pe = 1;

	//////////////////////////////////////////////////////////
	////			Variable for						  ////
	////		Coalescing Drop                           ////
	//////////////////////////////////////////////////////////
	double Oh = 1;
	double We = 1;
	double Bo = 1; // Bond number, the buoyancy force of the interfior fluid with respect to the exterior fluid.


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
	int examNum;
	int ProjectionOrder = 1;
	int spatialOrder = 3;
	int temporalOrder = 3;
	double cflCondition = 0.4;
	double dt;
	double finalT = 10;
	double totalT = 0;
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

	inline void DetermineViscosity();
	inline void DetermineDensity();
	inline void ComputeJJJJ(Array2D<double>& J11, Array2D<double>& J12, Array2D<double>& J21, Array2D<double>& J22);

	template <class TT>
	inline TT InterpolationGridtoU(const Array2D<TT>& ipData, const int& ui, const int& uj);
	template <class TT>
	inline TT InterpolationGridtoV(const Array2D<TT>& ipData, const int& vi, const int& vj);

	inline void ComputeSurfaceForce();
	inline void ComputeSurfaceForceUV();

	inline double AdaptiveTimeStep();

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

	inline void GenerateLinearSysteviscosityV2Order(CSR<double> & ipU_CSR, CSR<double> & ipV_CSR);
	inline void GenerateLinearSysteviscosityV2Order(VectorND<double>& vectorB);
	inline void GenerateLinearSystempPhi2Order(CSR<double> & ipCSR);
	inline void GenerateLinearSystempPhi2Order(VectorND<double>& vectorB);

	inline void SetLinearSystem(const int& iter);
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
	examNum = example;
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
		Density.dataArray = 1;
		DensityU = FD(gridU);
		DensityV = FD(gridV);
		Viscosity = FD(grid);
		Viscosity.dataArray = 1;
		ViscosityU = FD(gridU);
		ViscosityV = FD(gridV);
		Nu = FD(grid);

		Re = 1000;

		AdvectionU = FD(gridU);
		AdvectionV = FD(gridV);

		DiffusionU = FD(gridU);
		DiffusionV = FD(gridV);

		SurfaceForce = FD(grid);
		SurfaceForceX = FD(gridU);
		SurfaceForceY = FD(gridV);

		levelSet = LS(grid);
		levelSet.phi.dataArray = 1;

		ProjectionOrder = 1; 

		dt = 0;
		finalT = 30;
		totalT = 0;
		writeOutputIteration = 30;
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
		tempBC = 0;
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			for (int i = U.iStart; i <= U.iEnd; i++)
			{
				if (j == U.jStart || j == U.jEnd)
				{
					U.BC(i, j) = BC_DIR;
					continue;
				}
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
				if (j == V.jStart || j == V.jEnd)
				{
					V.BC(i, j) = BC_REFLECTION;
					continue;
				}
				if (i == V.iStart || i == V.iEnd)
				{
					V.BC(i, j) = BC_NEUM;
					continue;
				}
				V.BC(i, j) = tempBC++;
			}
		}
		originV = V;

		Density = FD(grid);
		Density.dataArray = 1;
		DensityU = FD(gridU);
		DensityV = FD(gridV);
		Viscosity = FD(grid);
		Viscosity.dataArray = 1;
		ViscosityU = FD(gridU);
		ViscosityV = FD(gridV);
		Nu = FD(grid);

		Re = 1000;

		AdvectionU = FD(gridU);
		AdvectionV = FD(gridV);

		DiffusionU = FD(gridU);
		DiffusionV = FD(gridV);

		SurfaceForce = FD(grid);
		SurfaceForceX = FD(gridU);
		SurfaceForceY = FD(gridV);

		levelSet = LS(grid);
		levelSet.phi.dataArray = 1;

		ProjectionOrder = 1;

		dt = 0;
		finalT = 4;
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
		/////  Example 1 or Example 4.3
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
				Pressure.BC(i, j) = tempBC++;
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
					U.BC(i, j) = BC_NEUM;
					continue;
				}
				if (j == U.jStart || j == U.jEnd)
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
		Density.dataArray = 1;
		DensityU = FD(gridU);
		DensityV = FD(gridV);
		Viscosity = FD(grid);
		Viscosity.dataArray = 1;
		ViscosityU = FD(gridU);
		ViscosityV = FD(gridV);
		Nu = FD(grid);


		AdvectionU = FD(gridU);
		AdvectionV = FD(gridV);

		DiffusionU = FD(gridU);
		DiffusionV = FD(gridV);

		SurfaceForce = FD(grid);
		SurfaceForceX = FD(gridU);
		SurfaceForceY = FD(gridV);

		levelSet = LS(grid);
		levelSet.phi.dataArray = 1;

		ProjectionOrder = 1;

		//cflCondition = 0.1;
		//dt = cflCondition*grid.dx;
		//finalT = 10;
		//maxIteration = int(finalT / dt);
		//totalT = 0;
		//writeOutputIteration = 10;
		//iteration = 0;

		isPlot = false;
	}

	if (example == 4)
	{
		cout << "*************************************************" << endl;
		cout << "    --- Navier-Stokes equation ---" << endl;
		cout << "    A Boundary Conditio nCapturing Method for" << endl;
		cout << "        Multiphase Incompressible Flow" << endl;
		cout << "           Example : Air Bubble" << endl;
		cout << "*************************************************" << endl;
		cout << endl;
		//////////////////////////////////
		// Large bubble or Small Bubble //
		//////////////////////////////////
		bool isSmallBubble = true;

		////////////////////////////////////////////////////////////////////////////
		// CSF model with delta function 
		// Boundary Condition Capturing method without delta function 
		///////////////////////////////////////////////////////////////////////////
		isCSFmodel = false;

		dimensionlessForm = false;
		isMultiPhase = true;
		isGravity = true;
		isPCG = true;

		double scaling = 1;
		if (isSmallBubble) scaling = 0.01;

		double xLength = 1 * scaling;
		double yLength = 2 * scaling;
		int domainRes = 80;
		grid = Grid2D(-xLength, xLength, domainRes + 1, -xLength, yLength, 1.5*domainRes + 1);
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
					//U.BC(i, j) = BC_NEUM;
					U.BC(i, j) = BC_DIR;
					continue;
				}
				if (j == U.jEnd)
				{
					//U.BC(i, j) = BC_NEUM;
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
					V.BC(i, j) = BC_DIR;
					continue;
				}
				if (j == V.jStart)
				{
					//V.BC(i, j) = BC_NEUM;
					V.BC(i, j) = BC_REFLECTION;
					continue;
				}
				if (j == V.jEnd)
				{
					//V.BC(i, j) = BC_NEUM;
					V.BC(i, j) = BC_REFLECTION;
					continue;
				}
				V.BC(i, j) = tempBC++;
			}
		}
		originV = V;

		levelSet = LS(grid);
		double x, y;
		double radius = 1.0 / 3.0 * scaling;
#pragma omp parallel for private(x, y)
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				x = grid(i, j).x, y = grid(i, j).y;
				levelSet(i, j) = sqrt(x*x + y*y) - radius;
			}
		}
		levelSet.InitialTube();
		levelSet.LComputeMeanCurvature();

		// Negative Level Set : Inside
		// Positive Level Set : Outside
		densityE = densityWater;
		densityI = densityAir;
		viscosityE = viscosityWater;
		viscosityI = viscosityAir;
		gamma0 = 0.7825;

		Density = FD(grid);
		Density.dataArray = 1;
		DensityU = FD(gridU);
		DensityV = FD(gridV);
		Viscosity = FD(grid);
		Viscosity.dataArray = 1;
		ViscosityU = FD(gridU);
		ViscosityV = FD(gridV);
		Nu = FD(grid);

		AdvectionU = FD(gridU);
		AdvectionV = FD(gridV);

		DiffusionU = FD(gridU);
		DiffusionV = FD(gridV);

		SurfaceForce = FD(grid);
		SurfaceForceX = FD(gridU);
		SurfaceForceY = FD(gridV);

		ProjectionOrder = 1;

		dt = 0;
		finalT = 1 * scaling;
		totalT = 0;
		writeOutputIteration = 30;
		iteration = 0;
	}

	if (example == 5)
	{
		/////////////////////////////////////////////////////////
		/////  Siviscositylations of surfactant effects
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
		Density.dataArray = 1;
		DensityU = FD(gridU);
		DensityV = FD(gridV);
		Viscosity = FD(grid);
		Viscosity.dataArray = 1;
		ViscosityU = FD(gridU);
		ViscosityV = FD(gridV);
		Nu = FD(grid);

		//Re = 1000;

		AdvectionU = FD(gridU);
		AdvectionV = FD(gridV);

		DiffusionU = FD(gridU);
		DiffusionV = FD(gridV);

		levelSet = LS(grid);
		levelSet.phi.dataArray = 1;

		

		isPlot = false;
		ProjectionOrder = 1;

	}


	if (ProjectionOrder == 1)
	{
		Pressure.CountNonZero();
		Pb = VectorND<double>(Pressure.num_all_full_cells);
		tempP = VectorND<double>(Pressure.num_all_full_cells);
		P_CSR = CSR<double>(Pressure.num_all_full_cells, Pressure.nnz);
	}
	else
	{
		Phi = FD(grid);
		Phixxyy = FD(grid);
		gradientPx = FD(U.innerGrid);
		gradientPy = FD(V.innerGrid);

		U.CountNonZero();
		Ub = VectorND<double>(U.num_all_full_cells);
		tempU = VectorND<double>(U.num_all_full_cells);

		V.CountNonZero();
		Vb = VectorND<double>(V.num_all_full_cells);
		tempV = VectorND<double>(V.num_all_full_cells);

		Phi.CountNonZero();
		Phib = VectorND<double>(Phi.num_all_full_cells);
		tempPhi = VectorND<double>(Phi.num_all_full_cells);
	}
	AdvectionMethod2D<double>::alpha = 1.5*grid.dx;
}

inline void FluidSolver2D::Solver(const int & example)
{
	bool writeFile = false;
	string fileName;
	string str;
	clock_t startTime = clock();
	double  before = 0, after = 0, timeCheck = 0;

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

	
	while (totalT <= finalT)
	{
		iteration++;
		before = clock();
		cout << endl;
		cout << "********************************" << endl;
		cout << "       Iteration " << to_string(iteration) << " : Start" << endl;
		totalT += (dt = AdaptiveTimeStep());

		TVDRK3TimeAdvection();

		if (isMultiPhase)
		{
			int reinitialIter = int(levelSet.gamma1 / min(levelSet.phi.dx, levelSet.phi.dy)) * 2;
			AdvectionMethod2D<double>::LLSPropagatingTVDRK3MACGrid(levelSet, U, V, dt, 3);
			AdvectionMethod2D<double>::LLSReinitializationTVDRK3(levelSet, dt, reinitialIter, 3);
			levelSet.UpdateInterface();
			levelSet.UpdateLLS();
		}
		cout << "       Iteration " << to_string(iteration) << " : End" << endl;
		cout << "********************************" << endl;

		if (iteration % 10 == 0 && isPlot)
		{
			MATLAB.Command("subplot(1,2,1)");
			PlotVelocity();
			MATLAB.Command("subplot(1,2,2)");
			Pressure.Variable("P");
			MATLAB.Command("surf(Xp, Yp, P), axis tight;");
			//MATLAB.WriteImage("fluid", iteration, "fig");
			MATLAB.WriteImage("fluid", iteration, "png");
		}


		timeCheck += ((after = clock()) - before) / CLOCKS_PER_SEC;
		cout << "Consuming Time : " + to_string((after - before) / CLOCKS_PER_SEC) + " / " + to_string(timeCheck) << endl;
	}


	if (writeFile)
	{
		fileName = "pressure" + to_string(iteration);
		Pressure.WriteFile(fileName);
		fileName = "xVelocity" + to_string(iteration);
		U.WriteFile(fileName);
		fileName = "yVelocity" + to_string(iteration);
		V.WriteFile(fileName);
	}
}

inline void FluidSolver2D::TVDRK3TimeAdvection()
{
	originU.dataArray = U.dataArray;
	originV.dataArray = V.dataArray;

	double* uVal(U.dataArray.values);
	double* vVal(V.dataArray.values);
	double* uOriginVal(originU.dataArray.values);
	double* vOriginVal(originV.dataArray.values);
	int uRes = U.dataArray.ijRes;
	int vRes = V.dataArray.ijRes;

	DetermineViscosity();
	DetermineDensity();
	SetLinearSystem(iteration);
	
	//// Compute Surface Force
	if (isMultiPhase)
	{
		if (isCSFmodel) ComputeSurfaceForceUV();
		else ComputeSurfaceForce();
	}
	//SurfaceForce.Variable("SurfaceForce");
	//MATLAB.Command("surf(Xp,Yp,SurfaceForce), axis equal tight;");
	/////////////////
	//// Step 1  ////
	/////////////////
	if (ProjectionOrder == 1)		 EulerMethod();
	else if (ProjectionOrder == 2)	 EulerMethod2ndOrder();

	/////////////////
	//// Step 2  ////
	/////////////////
	if (ProjectionOrder == 1)		 EulerMethod();
	else if (ProjectionOrder == 2)	 EulerMethod2ndOrder();

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
	if (ProjectionOrder == 1)		 EulerMethod();
	else if (ProjectionOrder == 2)	 EulerMethod2ndOrder();

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

inline void FluidSolver2D::EulerMethod()
{
	////////////////////////////////////////////////
	////     Projection Method 1 : advection    ////
	////////////////////////////////////////////////
	EulerMethodStep1();
	if (examNum == 1)
	{
		TreatBCAlongYaxis(U);
		TreatBCAlongXaxis(V);
	}
	else if (examNum == 2)
	{
		TreatBCAlongYaxis(U);
		TreatBCAlongYaxis(V);
		TreatBCAlongXaxis(V);
	}
	else if (examNum == 4)
	{
		TreatBCAlongXaxis(U);
		TreatBCAlongYaxis(U);
		TreatBCAlongYaxis(V);
		TreatBCAlongXaxis(V);
	}
	

	//U.Variable("Ustar");
	//V.Variable("Vstar");
	//MATLAB.Command("divUstar =Ustar(:,2:end)-Ustar(:,1:end-1),divVstar =Vstar(2:end,:)-Vstar(1:end-1,:);divstar=divUstar+divVstar;");
	//MATLAB.Command("quiver(Xp,Yp,Ustar(:,1:end-1),Vstar(1:end-1,:)), axis equal tight;");
	
	////////////////////////////////////////////////
	////     Projection Method 2 : Poisson Eq   ////
	////////////////////////////////////////////////
	EulerMethodStep2();
	//Pressure.Variable("Pressure");

	//////////////////////////////////////////////
	////     Projection Method 3 : New U,V    ////
	//////////////////////////////////////////////
	EulerMethodStep3();

	if (examNum == 1)
	{
		TreatBCAlongYaxis(U);
		TreatBCAlongXaxis(V);
	}
	else if (examNum == 2)
	{
		TreatBCAlongYaxis(U);
		TreatBCAlongYaxis(V);
		TreatBCAlongXaxis(V);
	}
	else if (examNum == 4)
	{
		TreatBCAlongXaxis(U);
		TreatBCAlongYaxis(U);
		TreatBCAlongYaxis(V);
		TreatBCAlongXaxis(V);
	}

	//U.Variable("Unew");
	//V.Variable("Vnew");
	//MATLAB.Command("divUnew =Unew(:,2:end)-Unew(:,1:end-1),divVnew =Vnew(2:end,:)-Vnew(1:end-1,:);divnew=divUnew+divVnew;");
	//MATLAB.Command("quiver(Xp,Yp,Unew(1:end-1,:),Vnew(:,1:end-1)), axis equal tight;");
	
}

inline void FluidSolver2D::EulerMethodStep1()
{
	//Array2D<double>& K1U = U.K1;
	//Array2D<double>& K1V = V.K1;


	AdvectionTerm(U, V, AdvectionU, AdvectionV);
	DiffusionTerm(U, V, DiffusionU, DiffusionV);
	//AdvectionU.Variable("advectionU");
	//DiffusionU.Variable("diffusionU");
	//AdvectionV.Variable("advectionV");
	//DiffusionV.Variable("diffusionV");
	int iStart = grid.iStart, jStart = grid.jStart;
	double oneOverRe = 1 / Re;
	double oneOverDensity;
	double ls;
	int ii, jj;
#pragma omp parallel for private(oneOverDensity, ii, ls)
	for (int i = U.iStart; i <= U.iEnd; i++)
	{
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			if (U.BC(i, j) < 0) continue;
			//K1U(i, j) = dt*(-AdvectionU(i, j) + Oh * oneOverRe*DiffusionU(i, j));
			if (dimensionlessForm)
			{
				U(i, j) += dt*(-AdvectionU(i, j) + Oh * oneOverRe*DiffusionU(i, j) + SurfaceForceX(i, j));
			}
			else
			{
				oneOverDensity = 1. / DensityU(i, j);
				U(i, j) += dt*(-AdvectionU(i, j) + oneOverDensity*(DiffusionU(i, j) + SurfaceForceX(i, j)));
			}
		}
	}
#pragma omp parallel for private(oneOverDensity, jj, ls)
	for (int i = V.iStart; i <= V.iEnd; i++)
	{
		for (int j = V.jStart; j <= V.jEnd; j++)
		{
			if (V.BC(i, j) < 0) continue;
			//K1V(i, j) = dt*(-AdvectionV(i, j) + Oh *oneOverRe*DiffusionV(i, j));
			if (dimensionlessForm)
			{
				V(i, j) += dt*(-AdvectionV(i, j) + Oh *oneOverRe*DiffusionV(i, j) + SurfaceForceY(i, j));
			}
			else
			{
				oneOverDensity = 1. / DensityV(i, j);
				V(i, j) += dt*(-AdvectionV(i, j) + oneOverDensity*(DiffusionV(i, j) + SurfaceForceY(i, j)));
			}
			if (isGravity) V(i, j) += dt*gravity;

		}
	}
	//K1U.Variable("k1u");
	//K1V.Variable("k1v");
}

inline void FluidSolver2D::EulerMethodStep2()
{
	GenerateLinearSystemPressure(Pb);
	Array2D<int> BC = Pressure.BC;
	double* pVal(tempP.values);
	int index;
#pragma omp parallel for private(index)
	for (int i = Pressure.iStart; i <= Pressure.iEnd; i++)
	{
		for (int j = Pressure.jStart; j <= Pressure.jEnd; j++)
		{
			index = BC(i, j);
			if (index < 0) continue;

			pVal[index] = Pressure(i, j);
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
	//Pressure.Variable("P");
	//MATLAB.Command("surf(Xp,Yp,P)");
}

inline void FluidSolver2D::EulerMethodStep3()
{
	double oneOverdx = Pressure.oneOverdx;
	double oneOverdy = Pressure.oneOverdy;
	double dx = Pressure.dx;
	double dy = Pressure.dy;
	int iL, iR, jB, jT, iStart = grid.iStart, iEnd = grid.iEnd, jStart = grid.jStart, jEnd = grid.jEnd;
	double lsL, lsR,  lsB, lsT;
	double JPL, JPR, JPB, JPT, JPx, JPy, JP;
	double oneOverDensity;
#pragma omp parallel for private(iL, iR, lsL, lsR, JPL, JPR, JP, oneOverDensity)
	for (int i = U.iStart; i <= U.iEnd; i++)
	{
		for (int j = U.jStart; j <= U.jEnd; j++)
		{
			if (U.BC(i, j) < 0) continue;

			iL = max(i - 1, iStart);
			iR = i;
			lsL = levelSet(iL, j);
			lsR = levelSet(iR, j);

			if (dimensionlessForm) U(i, j) -= dt * (Pressure(i, j) - Pressure(i - 1, j))*oneOverdx;
			else
			{
				if (isCSFmodel) JP = 0;
				else
				{
					if (lsL*lsR < 0)
					{
						JPL = SurfaceForce(iL, j);
						JPR = SurfaceForce(iR, j);
						JP = (JPL*abs(lsR) + JPR*abs(lsL)) / (abs(lsR) + abs(lsL));
						if (lsR > 0) JP *= -1;
					}
				}
				oneOverDensity = 1. / DensityU(i, j);
				U(i, j) -= dt * ((Pressure(i, j) - Pressure(iL, j) + JP) * oneOverdx) *oneOverDensity;
			}
		}
	}

#pragma omp parallel for private(jB, jT, lsB, lsT, JPB, JPT, JP, oneOverDensity)
	for (int i = V.iStart; i <= V.iEnd; i++)
	{
		for (int j = V.jStart; j <= V.jEnd; j++)
		{
			if (V.BC(i, j) < 0) continue;

			jB = max(j - 1, jStart);
			jT = j;
			lsB = levelSet(i, jB);
			lsT = levelSet(i, jT);

			if (dimensionlessForm) V(i, j) -= dt * (Pressure(i, j) - Pressure(i, j - 1))*oneOverdy;
			else
			{
				if (isCSFmodel) JP = 0;
				else
				{
					if (lsB*lsT < 0)
					{
						JPB = SurfaceForce(i, jB);
						JPT = SurfaceForce(i, jT);
						JP = (JPB*abs(lsT) + JPT*abs(lsB)) / (abs(lsT) + abs(lsB));
						if (lsT > 0) JP *= -1;
					}
				}
				oneOverDensity = 1. / DensityV(i, j);
				V(i, j) -= dt * ((Pressure(i, jT) - Pressure(i, jB) + JP)*oneOverdy) *oneOverDensity;
			}
		}
	}
}


inline void FluidSolver2D::GenerateLinearSystemPressure(CSR<double>& ipCSR)
{
	int iStart = Pressure.iStart, iEnd = Pressure.iEnd, jStart = Pressure.jStart, jEnd = Pressure.jEnd;

	double dx = Pressure.dx, dy = Pressure.dy;
	double dx2 = dx*dx, dy2 = dy*dy, dxdy = dx*dy;
	double oneOverdx = 1 / dx, oneOverdx2 = 1 / dx2;
	double oneOverdy = 1 / dy, oneOverdy2 = 1 / dy2;

	Array2D<int>& BC = Pressure.BC;

	double coefIJ;
	double betaL = 1, betaR = 1, betaB = 1, betaT = 1;
	double lsC = 1, lsL = 1, lsR = 1, lsB = 1, lsT = 1;
	double dC, dL, dR, dB, dT;
	for (int j = jStart; j <= jEnd; j++)
	{
		for (int i = iStart; i <= iEnd; i++)
		{
			coefIJ = 0;
			if (BC(i, j) < 0) continue;
		
			//// Density jump condition
			lsC = levelSet(i, j);
			if (i > iStart) lsL = levelSet(i - 1, j);
			else			lsL = lsC;
			if (i < iEnd)	lsR = levelSet(i + 1, j);
			else			lsR = lsC;
			if (j > jStart) lsB = levelSet(i, j - 1);
			else			lsB = lsC;
			if (j < jEnd)	lsT = levelSet(i, j + 1);
			else			lsT = lsC;

			dC = 1 / Density(i, j);
			if (i > iStart) dL = 1 / Density(i - 1, j);
			else			dL = dC;
			if (i < iEnd)	dR = 1 / Density(i + 1, j);
			else			dR = dC;
			if (j > jStart) dB = 1 / Density(i, j - 1);
			else			dB = dC;
			if (j < jEnd)	dT = 1 / Density(i, j + 1);
			else			dT = dC;

			betaL = dC*dL * (abs(lsL) + abs(lsC)) / (dC*abs(lsL) + dL*abs(lsC));
			betaR = dC*dR * (abs(lsR) + abs(lsC)) / (dC*abs(lsR) + dR*abs(lsC));
			betaB = dC*dB * (abs(lsB) + abs(lsC)) / (dC*abs(lsB) + dB*abs(lsC));
			betaT = dC*dT * (abs(lsT) + abs(lsC)) / (dC*abs(lsT) + dT*abs(lsC));


			//// If neighbor is full cell
			if (i>iStart)
			{
				if (BC(i - 1, j) > -1)
				{
					coefIJ += betaL;
					ipCSR.AssignValue(BC(i, j), BC(i - 1, j), - betaL * oneOverdx2);
				}
			}
			if (i<iEnd)
			{
				if (BC(i + 1, j) > -1)
				{
					coefIJ += betaR;
					ipCSR.AssignValue(BC(i, j), BC(i + 1, j), - betaR * oneOverdx2);
				}
			}
			if (j>jStart)
			{
				if (BC(i, j - 1) > -1)
				{
					coefIJ += betaB;
					ipCSR.AssignValue(BC(i, j), BC(i, j - 1), - betaB * oneOverdy2);
				}
			}
			if (j<jEnd)
			{
				if (BC(i, j + 1) > -1)
				{
					coefIJ += betaT;
					ipCSR.AssignValue(BC(i, j), BC(i, j + 1), - betaT * oneOverdy2);
				}
			}

			////// Dirichlet Boundary Condition
			//if (i>iStart)	if (BC(i - 1, j) == BC_DIR)	coefIJ += betaL;
			//if (i<iEnd)		if (BC(i + 1, j) == BC_DIR)	coefIJ += betaR;
			//if (j>jStart)	if (BC(i, j - 1) == BC_DIR)	coefIJ += betaB;
			//if (j<jEnd)		if (BC(i, j + 1) == BC_DIR)	coefIJ += betaT;

			//if (i == Pressure.iStartI && j == Pressure.jStartI)
			//{
			//	if (BC(i - 1, j) == BC_NEUM) coefIJ += betaL;
			//	if (BC(i + 1, j) == BC_NEUM) coefIJ += betaR;
			//	if (BC(i, j - 1) == BC_NEUM) coefIJ += betaB;
			//	if (BC(i, j + 1) == BC_NEUM) coefIJ += betaT;
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
	int iStart = Pressure.iStart, iEnd = Pressure.iEnd, jStart = Pressure.jStart, jEnd = Pressure.jEnd;
	double tempSum = 0;
	Array2D<int>& BC = Pressure.BC;
	double oneOverdt = 1. / dt;
	double oneOverdx = U.oneOverdx;
	double oneOverdy = U.oneOverdy;
	double oneOverdx2 = Pressure.oneOverdx2;
	double oneOverdy2 = Pressure.oneOverdy2;

	double* bVal(vectorB.values);
	int tempIndex;
	double betaL, betaR, betaB, betaT;
	double lsC, lsL, lsR, lsB, lsT; // Level Set
	double dC, dL, dR, dB, dT; // Density
	double jC, jL, jR, jB, jT; // Jump condition
	double fL, fR, fB, fT; // Force
	double aGamma, bGamma, theta;

	for (int j = jStart; j <= jEnd; j++)
	{
		for (int i = iStart; i <= iEnd; i++)
		{
			if (j==27)
			{
				if (i==35)
				{
					fL = 0;
				}
			}
			tempIndex = BC(i, j);
			if (tempIndex < 0) continue;
			
			bVal[tempIndex] = 0;
			
			fL = 0, fR = 0, fB = 0, fT = 0;

			//// Density jump condition
			lsC = levelSet(i, j);
			if (i > iStart) lsL = levelSet(i - 1, j);
			else			lsL = lsC;
			if (i < iEnd)	lsR = levelSet(i + 1, j);
			else			lsR = lsC;
			if (j > jStart) lsB = levelSet(i, j - 1);
			else			lsB = lsC;
			if (j < jEnd)	lsT = levelSet(i, j + 1);
			else			lsT = lsC;

			dC = 1. / Density(i, j);
			if (i > iStart) dL = 1. / Density(i - 1, j);
			else			dL = dC;
			if (i < iEnd)	dR = 1. / Density(i + 1, j);
			else			dR = dC;
			if (j > jStart) dB = 1. / Density(i, j - 1);
			else			dB = dC;
			if (j < jEnd)	dT = 1. / Density(i, j + 1);
			else			dT = dC;

			jC = SurfaceForce(i, j);
			if (i > iStart) jL = SurfaceForce(i - 1, j);
			else jL = jC;
			if (i < iEnd)	jR = SurfaceForce(i + 1, j);
			else jR = jC;
			if (j > jStart)	jB = SurfaceForce(i, j - 1);
			else jB = jC;
			if (j < jEnd)	jT = SurfaceForce(i, j + 1);
			else jT = jC;

			betaL = dC*dL * (abs(lsL) + abs(lsC)) / (dC*abs(lsL) + dL*abs(lsC));
			betaR = dC*dR * (abs(lsR) + abs(lsC)) / (dC*abs(lsR) + dR*abs(lsC));
			betaB = dC*dB * (abs(lsB) + abs(lsC)) / (dC*abs(lsB) + dB*abs(lsC));
			betaT = dC*dT * (abs(lsT) + abs(lsC)) / (dC*abs(lsT) + dT*abs(lsC));

			theta = abs(lsL) / (abs(lsL) + abs(lsC));
			aGamma = (jC*abs(lsL) + jL*abs(lsC)) / (abs(lsC) + abs(lsL));
			if (lsL > 0 && lsC <= 0)			 fL = oneOverdx2 * aGamma * betaL;
			else if (lsL < 0 && lsC >= 0)	 fL = -oneOverdx2 * aGamma * betaL;

			theta = abs(lsR) / (abs(lsR) + abs(lsC));
			aGamma = (jC*abs(lsR) + jR*abs(lsC)) / (abs(lsC) + abs(lsR));
			if (lsR > 0 && lsC <= 0)		 fR = oneOverdx2 * aGamma * betaR;
			else if (lsR < 0 && lsC >= 0)	 fR = -oneOverdx2 * aGamma * betaR;

			theta = abs(lsB) / (abs(lsB) + abs(lsC));
			aGamma = (jC*abs(lsB) + jB*abs(lsC)) / (abs(lsC) + abs(lsB));
			if (lsB > 0 && lsC <= 0)		 fB = oneOverdy2 * aGamma * betaB;
			else if (lsB < 0 && lsC >= 0)	 fB = -oneOverdy2 * aGamma * betaB;

			theta = abs(lsT) / (abs(lsT) + abs(lsC));
			aGamma = (jC*abs(lsT) + jT*abs(lsC)) / (abs(lsC) + abs(lsT));
			if (lsT > 0 && lsC <= 0)		 fT = oneOverdy2 * aGamma * betaT;
			else if (lsT < 0 && lsC >= 0)	 fT = -oneOverdy2 * aGamma * betaT;

			bVal[tempIndex] = -(fL + fR + fB + fT);

			bVal[tempIndex] += -oneOverdt * ((U(i + 1, j) - U(i, j))*oneOverdx + (V(i, j + 1) - V(i, j))*oneOverdy);

			if (i > iStart)	if (BC(i - 1, j) == BC_DIR)		bVal[tempIndex] += betaL * Pressure(i - 1, j)*oneOverdx2;
			if (i < iEnd)		if (BC(i + 1, j) == BC_DIR) bVal[tempIndex] += betaR * Pressure(i + 1, j)*oneOverdx2;
			if (j > jStart)	if (BC(i, j - 1) == BC_DIR)		bVal[tempIndex] += betaB * Pressure(i, j - 1)*oneOverdy2;
			if (j < jEnd)		if (BC(i, j + 1) == BC_DIR) bVal[tempIndex] += betaT * Pressure(i, j + 1)*oneOverdy2;

			tempSum += bVal[tempIndex];

			tempP.values[tempIndex] = Pressure(i, j);
		}
	}

	int bLength = vectorB.iLength;
	double ave = tempSum / (double)bLength;
#pragma omp parallel for
	for (int i = 0; i < bLength; i++)
	{
		bVal[i] -= ave;
	}
	//vectorB.Variable("b");
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
	double u;
#pragma omp parallel for private(aveV, Ux, Uy,u)
	for (int i = TermU.iStart; i <= TermU.iEnd; i++)
	{
		for (int j = TermU.jStart; j <= TermU.jEnd; j++)
		{
			if (U.BC(i, j) < 0) continue;

			aveV = 0.25 * (V(i - 1, j) + V(i, j) + V(i - 1, j + 1) + V(i, j + 1));
			u = U(i, j);
			if (u > 0)
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
			TermU(i, j) = u*Ux + aveV*Uy;
		}
	}

	double aveU;
	double Vx, Vy;
	double v;
#pragma omp parallel for private(aveU, Vx, Vy, v)
	for (int i = TermV.iStart; i <= TermV.iEnd; i++)
	{
		for (int j = TermV.jStart; j <= TermV.jEnd; j++)
		{
			if (V.BC(i, j) < 0) continue;

			aveU = 0.25 * (U(i, j - 1) + U(i + 1, j - 1) + U(i, j) + U(i + 1, j));
			if (aveU > 0)
			{
				Vx = dVdxM(i, j);
			}
			else
			{
				Vx = dVdxP(i, j);
			}
			v = V(i, j);
			if (v > 0)
			{
				Vy = dVdyM(i, j);
			}
			else
			{
				Vy = dVdyP(i, j);
			}
			TermV(i, j) = aveU*Vx + v*Vy;
		}
	}
}

inline void FluidSolver2D::DiffusionTerm(const FD & U, const FD & V, FD & TermU, FD & TermV)
{
	double dx = U.dx, dy = U.dy;
	double oneOverdx = U.oneOverdx, oneOverdy = U.oneOverdy;
	double lsL, lsR, lsB, lsT;
	double muL, muR, muB, muT;
	if (isCSFmodel)
	{
		int iStart = grid.iStart, iEnd = grid.iEnd, jStart = grid.jStart, jEnd = grid.jEnd;
#pragma omp parallel for private(lsL, lsR, lsB, lsT, muL, muR, muB, muT)
		for (int i = U.iStart; i <= U.iEnd; i++)
		{
			for (int j = U.jStart; j <= U.jEnd; j++)
			{
				if (U.BC(i, j) < 0) continue;
				lsL = levelSet(max(i - 1,iStart), j);
				lsR = levelSet(min(i, iEnd), j);
				lsB = 0.25*(lsL + lsR + levelSet(max(i - 1, iStart), max(j - 1, jStart)) + levelSet(i, max(j - 1, jStart)));
				lsT = 0.25*(lsL + lsR + levelSet(max(i - 1, iStart), min(j + 1, jEnd)) + levelSet(i, min(j + 1, jEnd)));
				muL = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(lsL);
				muR = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(lsR);
				muB = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(lsB);
				muT = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(lsT);

				TermU(i, j) = 2 * (muR*(U(i + 1, j) - U(i, j))*oneOverdx - muL*(U(i, j) - U(i - 1, j))*oneOverdx)*oneOverdx
					+ (muT * ((U(i, j + 1) - U(i, j))*oneOverdy + (V(i, j + 1) - V(i - 1, j + 1))*oneOverdx)
						- muB * ((U(i, j) - U(i, j - 1))*oneOverdy + (V(i, j) - V(i - 1, j))*oneOverdx))*oneOverdy;
			}
		}
#pragma omp parallel for private(lsL, lsR, lsB, lsT, muL, muR, muB, muT)
		for (int i = V.iStart; i <= V.iEnd; i++)
		{
			for (int j = V.jStart; j <= V.jEnd; j++)
			{
				if (V.BC(i, j) < 0) continue;
				lsB = levelSet(i, max(j - 1, jStart));
				lsT = levelSet(i, min(j, jEnd));
				lsL = 0.25 * (lsB + lsT + levelSet(max(i - 1, iStart), j) + levelSet(max(i - 1, iStart), max(j - 1, jStart)));
				lsR = 0.25 * (lsB + lsT + levelSet(min(i + 1,iEnd), j) + levelSet(min(i + 1, iEnd), max(j - 1, jStart)));
				muL = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(lsL);
				muR = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(lsR);
				muB = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(lsB);
				muT = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(lsT);

				TermV(i, j) = 2 * (muT * (V(i, j + 1) - V(i, j))*oneOverdy - muB * (V(i, j) - V(i, j - 1))*oneOverdy)*oneOverdy
					+ ( muR * ((U(i + 1, j) - U(i + 1, j - 1))*oneOverdy + (V(i + 1, j) - V(i, j))*oneOverdx)
						- muL * ((U(i, j) - U(i, j - 1))*oneOverdy + (V(i, j) - V(i - 1, j))*oneOverdx))*oneOverdx;
			}
		}
	}
	else
	{
		Array2D<double> J11(gridU), J12(gridU), J21(gridV), J22(gridV);
		ComputeJJJJ(J11, J12, J21, J22);

		int iStartL = grid.iStart, iEndL = grid.iEnd, jStartL = grid.jStart, jEndL = grid.jEnd;
		int iStart = U.iStart, iEnd = U.iEnd, jStart = U.jStart, jEnd = U.jEnd;

		double lsL, lsC, lsR, lsB, lsT;
		double viscoL, viscoC, viscoR, viscoB, viscoT;
		double uL, uC, uR, uI, uB, uT;
		double theta;
		double JL, JC, JR, JI, JB, JT;
		double densityJ = viscosityE - viscosityI;
		double mudUdxL, mudUdxR, mudUdxB, mudUdxT;
#pragma omp parallel for private(lsL, lsC, lsR, lsB, lsT, viscoL, viscoC, viscoR, viscoB, viscoT, uL, uC, uR, uI, uB, uT, theta, JL, JC, JR, JI, JB, JT, mudUdxL, mudUdxR, mudUdxB, mudUdxT)
		for (int i = U.iStart; i <= U.iEnd; i++)
		{
			for (int j = U.jStart; j <= U.jEnd; j++)
			{
				if (U.BC(i, j) < 0) continue;

				TermU(i, j) = 0;

				lsL = InterpolationGridtoU(levelSet.phi.dataArray, i - 1, j);
				lsC = InterpolationGridtoU(levelSet.phi.dataArray, i, j);
				lsR = InterpolationGridtoU(levelSet.phi.dataArray, i + 1, j);
				lsB = InterpolationGridtoU(levelSet.phi.dataArray, i, j - 1);
				lsT = InterpolationGridtoU(levelSet.phi.dataArray, i, j + 1);

				viscoL = ViscosityU(i - 1, j);
				viscoC = ViscosityU(i, j);
				viscoR = ViscosityU(i + 1, j);
				viscoB = ViscosityU(i, j - 1);
				viscoT = ViscosityU(i, j + 1);

				uL = U(i - 1, j);
				uC = U(i, j);
				uR = U(i + 1, j);
				uB = U(i, j - 1);
				uT = U(i, j + 1);

				//////////////////////
				//// Compute Uxx
				//////////////////////
				JL = J11(i - 1, j);
				JC = J11(i, j);
				JR = J11(i + 1, j);

				if (lsL *lsC >= 0)
				{
					mudUdxL = viscoC*(uC - uL)*oneOverdx;
				}
				else
				{
					theta = abs(lsL) / (abs(lsL) + abs(lsC));
					JI = theta*JC + (1 - theta)*JL;
					uI = (viscoC*uC*theta + viscoL*uL*(1 - theta) - JI*theta*(1 - theta)*dx) / (viscoC*theta + viscoL*(1 - theta));
					mudUdxL = viscoC*(uC - uI) / ((1 - theta)*dx + DBL_EPSILON);
				}

				if (lsC *lsR >= 0)
				{
					mudUdxR = viscoC*(uR - uC)*oneOverdx;
				}
				else
				{
					theta = abs(lsR) / (abs(lsR) + abs(lsC));
					JI = theta*JC + (1 - theta)*JR;
					uI = (viscoC*uC*theta + viscoR*uR*(1 - theta) - JI*theta*(1 - theta)*dx) / (viscoC*theta + viscoR*(1 - theta));
					mudUdxR = viscoC*(uI - uC) / ((1 - theta)*dx + DBL_EPSILON);
				}
				TermU(i, j) += (mudUdxR - mudUdxL)*oneOverdx;

				//////////////////////
				//// Compute Uyy
				/////////////////////
				JC = J12(i, j);
				JB = J12(i, j - 1);
				JT = J12(i, j + 1);

				if (lsB *lsC >= 0)
				{
					mudUdxB = viscoC*(uC - uB)*oneOverdy;
				}
				else
				{
					theta = abs(lsB) / (abs(lsB) + abs(lsC));
					JI = theta*JC + (1 - theta)*JB;
					uI = (viscoC*uC*theta + viscoB*uB*(1 - theta) - JI*theta*(1 - theta)*dx) / (viscoC*theta + viscoB*(1 - theta));
					mudUdxB = viscoC*(uC - uI) / ((1 - theta)*dy + DBL_EPSILON);
				}

				if (lsC *lsT >= 0)
				{
					mudUdxT = viscoC*(uT - uC)*oneOverdy;
				}
				else
				{
					theta = abs(lsT) / (abs(lsT) + abs(lsC));
					JI = theta*JC + (1 - theta)*JT;
					uI = (viscoC*uC*theta + viscoT*uT*(1 - theta) - JI*theta*(1 - theta)*dx) / (viscoC*theta + viscoT*(1 - theta));
					mudUdxT = viscoC*(uI - uC) / ((1 - theta)*dy + DBL_EPSILON);
				}
				TermU(i, j) += (mudUdxT - mudUdxB)*oneOverdy;

			}
		}

		iStart = V.iStart, iEnd = V.iEnd, jStart = V.jStart, jEnd = V.jEnd;
		double vL, vC, vR, vI, vB, vT;
		double mudVdxL, mudVdxR, mudVdxB, mudVdxT;
#pragma omp parallel for private(lsL, lsC, lsR, lsB, lsT, viscoL, viscoC, viscoR, viscoB, viscoT, vL, vC, vR, vI, vB, vT, theta, JL, JC, JR, JI, JB, JT, mudVdxL, mudVdxR, mudVdxB, mudVdxT)
		for (int i = iStart; i <= iEnd; i++)
		{
			for (int j = jStart; j <= jEnd; j++)
			{
				if (V.BC(i, j) < 0) continue;

				TermV(i, j) = 0;

				lsL = InterpolationGridtoV(levelSet.phi.dataArray, i - 1, j);
				lsC = InterpolationGridtoV(levelSet.phi.dataArray, i, j);
				lsR = InterpolationGridtoV(levelSet.phi.dataArray, i + 1, j);
				lsB = InterpolationGridtoV(levelSet.phi.dataArray, i, j - 1);
				lsT = InterpolationGridtoV(levelSet.phi.dataArray, i, j + 1);

				viscoL = ViscosityV(i - 1, j);
				viscoC = ViscosityV(i, j);
				viscoR = ViscosityV(i + 1, j);
				viscoB = ViscosityV(i, j - 1);
				viscoT = ViscosityV(i, j + 1);

				vL = V(i - 1, j);
				vC = V(i, j);
				vR = V(i + 1, j);
				vB = V(i, j - 1);
				vT = V(i, j + 1);

				//////////////////////
				//// Compute Vxx
				//////////////////////
				JL = J21(i - 1, j);
				JC = J21(i, j);
				JR = J21(i + 1, j);

				if (lsL *lsC >= 0)
				{
					mudVdxL = viscoC*(vC - vL)*oneOverdx;
				}
				else
				{
					theta = abs(lsL) / (abs(lsL) + abs(lsC));
					JI = theta*JC + (1 - theta)*JL;
					vI = (viscoC*vC*theta + viscoL*vL*(1 - theta) - JI*theta*(1 - theta)*dx) / (viscoC*theta + viscoL*(1 - theta));
					mudVdxL = viscoC*(vC - vI) / ((1 - theta)*dx + DBL_EPSILON);
				}

				if (lsC *lsR >= 0)
				{
					mudVdxR = viscoC*(vR - vC)*oneOverdx;
				}
				else
				{
					theta = abs(lsR) / (abs(lsR) + abs(lsC));
					JI = theta*JC + (1 - theta)*JR;
					vI = (viscoC*vC*theta + viscoR*vR*(1 - theta) - JI*theta*(1 - theta)*dx) / (viscoC*theta + viscoR*(1 - theta));
					mudVdxR = viscoC*(vI - vC) / ((1 - theta)*dx + DBL_EPSILON);
				}
				TermV(i, j) += (mudVdxR - mudVdxL)*oneOverdx;

				//////////////////////
				//// Compute Vyy
				/////////////////////
				JC = J22(i, j);
				JB = J22(i, j - 1);
				JT = J22(i, j + 1);

				if (lsB *lsC >= 0)
				{
					mudVdxB = viscoC*(vC - vB)*oneOverdy;
				}
				else
				{
					theta = abs(lsB) / (abs(lsB) + abs(lsC));
					JI = theta*JC + (1 - theta)*JB;
					vI = (viscoC*vC*theta + viscoB*vB*(1 - theta) - JI*theta*(1 - theta)*dx) / (viscoC*theta + viscoB*(1 - theta));
					mudVdxB = viscoC*(vC - vI) / ((1 - theta)*dy + DBL_EPSILON);
				}

				if (lsC *lsT >= 0)
				{
					mudVdxT = viscoC*(vT - vC)*oneOverdy;
				}
				else
				{
					theta = abs(lsT) / (abs(lsT) + abs(lsC));
					JI = theta*JC + (1 - theta)*JT;
					vI = (viscoC*vC*theta + viscoT*vT*(1 - theta) - JI*theta*(1 - theta)*dx) / (viscoC*theta + viscoT*(1 - theta));
					mudVdxT = viscoC*(vI - vC) / ((1 - theta)*dy + DBL_EPSILON);
				}
				TermV(i, j) += (mudVdxT - mudVdxB)*oneOverdy;
			}
		}
	}

}

// Interpolation Viscosity using Level Set and Heaviside function.
inline void FluidSolver2D::DetermineViscosity()
{
	int iStart = grid.iStart, iEnd = grid.iEnd, jStart = grid.jStart, jEnd = grid.jEnd;
	int UiStart = U.iStart, UiEnd = U.iEnd, UjStart = U.jStart, UjEnd = U.jEnd;
	int ViStart = V.iStart, ViEnd = V.iEnd, VjStart = V.jStart, VjEnd = V.jEnd;
	if (isMultiPhase)
	{
		// Level Set Viscosity.
#pragma omp parallel for
		for (int i = iStart; i <= iEnd; i++)
		{
			for (int j = jStart; j <= jEnd; j++)
			{
				if (levelSet(i, j) < 0) Viscosity(i, j) = viscosityI;
				else					Viscosity(i, j) = viscosityE;
			}
		}

		int iL, iR;
		double ls;
		// U Viscosity
#pragma omp parallel for private (iL, iR, ls)
		for (int i = UiStart; i <= UiEnd; i++)
		{
			for (int j = UjStart; j <= UjEnd; j++)
			{
				iL = max(i - 1, iStart);
				iR = min(i, iEnd);
				ls = 0.5*(levelSet(iL, j) + levelSet(iR, j));
				ViscosityU(i, j) = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(ls);
			}
		}
		
		int jB, jT;
		// V Viscosity
#pragma omp parallel for private(jB, jT, ls)
		for (int i = ViStart; i <= ViEnd; i++)
		{
			for (int j = VjStart; j <= VjEnd; j++)
			{
				jB = max(j - 1, jStart);
				jT = min(j, jEnd);
				ls = 0.5*(levelSet(i, jB) + levelSet(i, jT));
				ViscosityV(i, j) = viscosityI + (viscosityE - viscosityI)*AdvectionMethod2D<double>::Heaviside(ls);
			}
		}
	}
	else
	{
		Viscosity.dataArray = viscosityI;
		ViscosityU.dataArray = viscosityI;
		ViscosityV.dataArray = viscosityI;
	}
}

// Interpolation Density using Level Set and Heaviside function.
inline void FluidSolver2D::DetermineDensity()
{
	int iStart = grid.iStart, iEnd = grid.iEnd, jStart = grid.jStart, jEnd = grid.jEnd;
	int UiStart = U.iStart, UiEnd = U.iEnd, UjStart = U.jStart, UjEnd = U.jEnd;
	int ViStart = V.iStart, ViEnd = V.iEnd, VjStart = V.jStart, VjEnd = V.jEnd;
	if (isMultiPhase)
	{
		// Level Set Density.
#pragma omp parallel for 
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (levelSet(i, j) < 0) Density(i, j) = densityI;
				else					Density(i, j) = densityE;
			}
		}

		int iL, iR;
		double ls;
		// U Density
#pragma omp parallel for private (iL, iR, ls)
		for (int i = UiStart; i <= UiEnd ; i++)
		{
			for (int j = UjStart; j <= UjEnd; j++)
			{
				iL = max(i - 1, iStart);
				iR = min(i, iEnd);
				ls = 0.5*(levelSet(iL, j) + levelSet(iR, j));
				DensityU(i, j) = densityI + (densityE - densityI)*AdvectionMethod2D<double>::Heaviside(ls);
			}
		}

		int jB, jT;
		// V Density
#pragma omp parallel for private(jB, jT, ls)
		for (int i = ViStart; i <= ViEnd; i++)
		{
			for (int j = VjStart; j <= VjEnd; j++)
			{
				jB = max(j - 1, jStart);
				jT = min(j, jEnd);
				ls = 0.5*(levelSet(i, jB) + levelSet(i, jT));
				DensityV(i, j) = densityI + (densityE - densityI)*AdvectionMethod2D<double>::Heaviside(ls);
			}
		}
	}
	else
	{
		Density.dataArray = densityI;
		DensityU.dataArray = densityI;
		DensityV.dataArray = densityI;
	}
}

// Compute Viscosity Jump Condition.
inline void FluidSolver2D::ComputeJJJJ(Array2D<double> & J11, Array2D<double> & J12, Array2D<double> & J21, Array2D<double> & J22)
{
	levelSet.ComputeUnitNormal();
	Array2D<VT>& UnitNormal = levelSet.unitNormal.dataArray;
	Array2D<VT> UnitTangent(grid);

#pragma omp parallel for
	for (int i = 0; i < grid.iRes*grid.jRes; i++)
	{
		UnitTangent.values[i].x = UnitNormal.values[i].y;
		UnitTangent.values[i].y = -UnitNormal.values[i].x;
	}
	
	const double viscosityJump = viscosityE - viscosityI;
	VT  normal, tangent;
	double oneOverdx = U.oneOverdx, oneOverdy = U.oneOverdy;
	double dx = U.dx, dy = U.dy;

	double n1, n2, t1, t2;
	double ux, uy, vx, vy;
	int iStart = U.iStart, iEnd = U.iEnd, jStart = U.jStart, jEnd = U.jEnd;
#pragma omp parallel for private(normal, tangent, n1, n2, t1, t2, ux, uy, vx, vy)
	for (int j = jStart; j <= jEnd; j++)
	{
		for (int i = iStart; i <= iEnd; i++)
		{
			if (U.BC(i, j) < 0) continue;
			J11(i, j) = 0;
			J12(i, j) = 0;

			normal = InterpolationGridtoU(UnitNormal, i, j);
			normal.normalize();
			n1 = normal.x, n2 = normal.y;
			tangent = InterpolationGridtoU(UnitTangent, i, j);
			t1 = tangent.x, t2 = tangent.y;

			ux = U.dxPhi(i, j);
			uy = U.dyPhi(i, j);
			vx = 0.25*(V.dxPhi(i - 1, j) + V.dxPhi(i, j) + V.dxPhi(i - 1, j + 1) + V.dxPhi(i, j + 1));
			vy = 0.25*(V.dyPhi(i - 1, j) + V.dyPhi(i, j) + V.dyPhi(i - 1, j + 1) + V.dyPhi(i, j + 1));

			J11(i, j) = t1*t1*ux + t1*t2*uy + (n1*n1*ux + n1*n2*vx)*n1*n1 + (n1*n1*uy + n1*n2*vy)*n1*n2
				- (t1*t1*ux + t1*t2*uy)*n1*n1 - (t1*t1*vx + t1*t2*vy)*n1*n2;
			J11(i, j) *= viscosityJump;

			J12(i, j) = t1*t2*ux + t2*t2*uy + (n1*n1*ux + n1*n2*vx)*n1*n2 + (n1*n1*uy + n1*n2*vy)*n2*n2
				- (t1*t1*ux + t1*t2*uy)*n1*n2 - (t1*t1*vx + t1*t2*vy)*n2*n2;
			J12(i, j) *= viscosityJump;
		}
	}

	iStart = V.iStart, iEnd = V.iEnd, jStart = V.jStart, jEnd = V.jEnd;
#pragma omp parallel for private(normal, tangent, n1, n2, t1, t2, ux, uy, vx, vy)
	for (int j = jStart; j <= jEnd; j++)
	{
		for (int i = iStart; i <= iEnd; i++)
		{
			if (V.BC(i, j) < 0) continue;
			J21(i, j) = 0;
			J22(i, j) = 0;

			normal = InterpolationGridtoV(UnitNormal, i, j);
			normal.normalize();
			n1 = normal.x, n2 = normal.y;
			tangent = InterpolationGridtoV(UnitTangent, i, j);
			t1 = tangent.x, t2 = tangent.y;

			ux = 0.25*(U.dxPhi(i, j - 1) + U.dxPhi(i + 1, j - 1) + U.dxPhi(i, j) + U.dxPhi(i + 1, j));
			uy = 0.25*(U.dyPhi(i, j - 1) + U.dyPhi(i + 1, j - 1) + U.dyPhi(i, j) + U.dyPhi(i + 1, j));
			vx = V.dxPhi(i, j);
			vy = V.dyPhi(i, j);

			J21(i, j) = t1*t1*vx + t1*t2*vy + (n1*n2*ux + n2*n2*vx)*n1*n1 + (n1*n2*uy + n2*n2*vy)*n1*n2
				- (t1*t2*ux + t2*t2*uy)*n1*n1 - (t1*t2*vx + t2*t2*vy)*n1*n2;
			J21(i, j) *= viscosityJump;
			J22(i, j) = t1*t2*vx + t2*t2*vy + (n1*n2*ux + n2*n2*vx)*n1*n2 + (n1*n2*uy + n2*n2*vy)*n2*n2
				-(t1*t2*ux + t2*t2*uy)*n1*n2 - (t1*t2*vx + t2*t2*vy)*n2*n2;
			J22(i, j) *= viscosityJump;
		}
	}


}

// Interpolation Data using Level Set and Internally Dividign Point.
template <class TT>
inline TT FluidSolver2D::InterpolationGridtoU(const Array2D<TT>& ipData, const int & ui, const int & uj)
{
	int iStart = grid.iStart, iEnd = grid.iEnd, jStart = grid.jStart, jEnd = grid.jEnd;
	int tempJ;

	if (uj >= jStart && uj <= jEnd)
	{
		tempJ = uj;
	}
	if (uj < jStart)
	{
		tempJ = jStart;
	}
	if (uj > jEnd)
	{
		tempJ = jEnd;
	}
	if (ui >= iStart + 1 && ui <= iEnd)
	{
		double theta = abs(levelSet(ui - 1, tempJ)) / (abs(levelSet(ui - 1, tempJ)) + abs(levelSet(ui, tempJ)));
		return (1 - theta)*ipData(ui - 1, tempJ) + theta*ipData(ui, tempJ);
	}
	else if (ui <= iStart)  return ipData(iStart, tempJ);
	else					return ipData(iEnd, tempJ);
}

// Interpolation Data using Level Set and Internally Dividign Point.
template <class TT>
inline TT FluidSolver2D::InterpolationGridtoV(const Array2D<TT>& ipData, const int & vi, const int & vj)
{
	int iStart = grid.iStart, iEnd = grid.iEnd, jStart = grid.jStart, jEnd = grid.jEnd;
	int tempI;

	if (vi >= iStart && vi <= iEnd)
	{
		tempI = vi;
	}
	if (vi < iStart)
	{
		tempI = iStart;
	}
	if (vi > iEnd)
	{
		tempI = iEnd;
	}

	if (vj >= jStart + 1 && vj <= jEnd)
	{
		double theta = abs(levelSet(tempI, vj - 1)) / (abs(levelSet(tempI, vj - 1)) + abs(levelSet(tempI, vj)));
		return (1 - theta)*ipData(tempI, vj - 1) + theta*ipData(tempI, vj);
	}
	else if (vj <= jStart)	return ipData(tempI, jStart);
	else 					return ipData(tempI, jEnd);

}

inline void FluidSolver2D::ComputeSurfaceForce()
{
	//////////////////////////
	// L : Level Set
	// S : Surfactant
	// G : Gradient
	// U : Unit
	// N : Normal
	//////////////////////////
	Array2D<double>& meanCurvature = levelSet.meanCurvature.dataArray;
	Array2D<VT>& LunitNormal = levelSet.unitNormal.dataArray;
	Array2D<int>& tube = levelSet.tube;


	int ComputedTubeRange = 1;

	levelSet.LComputeMeanCurvature(ComputedTubeRange);
	levelSet.LComputeUnitNormal(ComputedTubeRange);
	U.Gradient();
	V.Gradient();
	
	double oneOverReCa = 1 / (Re*Ca);
	double oneOverrho;
	int numTube = levelSet.numTube;
	const double viscosityJump = viscosityE - viscosityI;
	int iEndu = grid.iEnd, jEndv = grid.jEnd;
	int i, ii, j, jj;
	double curvature;
	VT normal, GU, GV;
#pragma omp parallel for private(i, ii, j, jj, curvature, normal, GU, GV)
	for (int k = 1; k <= numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		if (tube(i, j) <= ComputedTubeRange)
		{
			curvature = -meanCurvature(i, j);
			if (dimensionlessForm) SurfaceForce(i, j) = gamma0 * curvature * oneOverReCa;
			else
			{
				normal = LunitNormal(i, j);
				ii = min(i + 1, iEndu);
				jj = min(j + 1, jEndv);

				GU = 0.5 * (U.gradient(i, j) + U.gradient(ii, j));
				GV = 0.5 * (V.gradient(i, j) + V.gradient(i, jj));
				
				SurfaceForce(i, j) = gamma0 * curvature;
				SurfaceForce(i, j) += 2*viscosityJump * ((GU.x*normal.x + GU.y*normal.y)*normal.x+ (GV.x*normal.x + GV.y*normal.y)*normal.y);
			}
		}
		else
		{
			SurfaceForce(i, j) = 0;
		}
	}
}

inline void FluidSolver2D::ComputeSurfaceForceUV()
{
	//////////////////////////
	// L : Level Set
	// S : Surfactant
	// G : Gradient
	// U : Unit
	// N : Normal
	//////////////////////////
	Array2D<double>& meanCurvature = levelSet.meanCurvature.dataArray;
	Array2D<VT>& LunitNormal = levelSet.unitNormal.dataArray;
	Array2D<int>& tube = levelSet.tube;

	int ComputedTubeRange = 1;

	levelSet.LComputeMeanCurvature(ComputedTubeRange);
	levelSet.LComputeUnitNormal(ComputedTubeRange);

	double oneOverReCa = 1; 1 / (Re*Ca);
	int numTube = levelSet.numTube;
	int iStart = grid.iStart, jStart = grid.jStart;
	double oneOverdx = grid.oneOverdx, oneOverdy = grid.oneOverdy;
	int i, j, ii, jj;
	VT LUN, LG;
	double LGMag, deltaL, curvature, ST;
	double Lval, LvalLeft, LvalBottom;
	ST = gamma0;
#pragma omp parallel for private(i, j, ii, jj, LUN, LG, LGMag, deltaL, curvature, Lval, LvalLeft, LvalBottom)
	for (int k = 1; k <= numTube; k++)
	{
		levelSet.TubeIndex(k, i, j);
		ii = max(i - 1, iStart);
		jj = max(j - 1, jStart);
		Lval = levelSet(i, j), LvalLeft = levelSet(ii, j), LvalBottom = levelSet(i, jj);
		if (tube(i, j) == ComputedTubeRange || tube(ii, j) == ComputedTubeRange)
		{
			/////////////////////////////////
			// SurfaceForceX on MAC grid.  //
			/////////////////////////////////
			LG.x = (Lval - LvalLeft) * oneOverdx;
			LG.y = 0.5 * (levelSet.dyPhi(ii, j) + levelSet.dyPhi(i, j));
			LGMag = LG.magnitude();
			LUN = LG / LGMag;
			deltaL = AdvectionMethod2D<double>::DeltaFt(0.5 * (Lval + LvalLeft));
			curvature = 0.5 * (- meanCurvature(i, j) - meanCurvature(ii, j));

			SurfaceForceX(i, j) = curvature * ST * LUN.x;
			if (dimensionlessForm) SurfaceForceX(i, j) *= deltaL * LGMag * oneOverReCa;
			else SurfaceForceX(i, j) *= deltaL;
		}
		if (tube(i, j) == ComputedTubeRange || tube(i, jj) == ComputedTubeRange)
		{
			/////////////////////////////////
			// SurfaceForceY on MAC grid.  //
			/////////////////////////////////
			LG.x = 0.5 * (levelSet.dxPhi(i, j) + levelSet.dxPhi(i, jj));
			LG.y = (Lval - LvalBottom) * oneOverdy;
			LGMag = LG.magnitude();
			LUN = LG / LGMag;
			deltaL = AdvectionMethod2D<double>::DeltaFt(0.5 * (Lval + LvalBottom));
			curvature = 0.5 * (- meanCurvature(i, j) - meanCurvature(i, jj));

			SurfaceForceY(i, j) = curvature * ST * LUN.y;
			if (dimensionlessForm) SurfaceForceY(i, j) *= deltaL * LGMag * oneOverReCa;
			else SurfaceForceY(i, j) *= deltaL;
		}
		else
		{
			SurfaceForceX(i, j) = 0;
			SurfaceForceY(i, j) = 0;
		}
	}
}

inline double FluidSolver2D::AdaptiveTimeStep()
{
	int iStart = grid.iStart, iEnd = grid.iEnd, jStart = grid.jStart, jEnd = grid.jEnd;
	int UiStart = U.iStart, UiEnd = U.iEnd, UjStart = U.jStart, UjEnd = U.jEnd;
	int ViStart = V.iStart, ViEnd = V.iEnd, VjStart = V.jStart, VjEnd = V.jEnd;
	double dx = grid.dx, dy = grid.dy;
	double Ccfl = 1, Vcfl = 1, Gcfl = 0, Scfl = 1;
	double uMax = 0, vMax = 0, curvatureMax = 0;

	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j < jEnd; j++)
		{
			if (abs(U(i, j)) > uMax) uMax = abs(U(i, j));
			if (abs(V(i, j)) > vMax) vMax = abs(V(i, j));
			if (abs(levelSet.meanCurvature(i, j)) > curvatureMax) curvatureMax = abs(levelSet.meanCurvature(i, j));
		}
	}
	
	Ccfl = (uMax / dx + vMax / dy);

	Vcfl = max(viscosityE / densityE, viscosityI / densityI)*sqrt(2 / (dx*dx) + 2 / (dx*dx));

	if (isGravity) Gcfl = sqrt(abs(gravity) / dy);

	Scfl = sqrt(gamma0*curvatureMax / (min(densityE, densityI)*min(dx*dx, dy*dy)));

	double finalCFL = ((Ccfl + Vcfl) + sqrt((Ccfl + Vcfl)*(Ccfl + Vcfl) + 4 * Gcfl*Gcfl + 4 * Scfl*Scfl)) / 2;

	return 0.5/finalCFL;
}

inline void FluidSolver2D::TreatBCAlongXaxis(FD & ipField)
{
	Array2D<int>& BC = ipField.BC;
	Array2D<double>& fieldData = ipField.dataArray;
#pragma omp parallel for 
	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
	{
		if (BC(i, ipField.jStart) == BC_NEUM)		fieldData(i, ipField.jStart) = fieldData(i, ipField.jStartI);
		if (BC(i, ipField.jStart) == BC_REFLECTION) fieldData(i, ipField.jStart) = -fieldData(i, ipField.jStartI);
		if (BC(i, ipField.jEnd) == BC_NEUM)			fieldData(i, ipField.jEnd)   = fieldData(i, ipField.jEndI);
		if (BC(i, ipField.jEnd) == BC_REFLECTION)	fieldData(i, ipField.jEnd)   = -fieldData(i, ipField.jEndI);
	}
}

inline void FluidSolver2D::TreatBCAlongYaxis(FD & ipField)
{
	Array2D<int>& BC = ipField.BC;
	Array2D<double>& fieldData = ipField.dataArray;
#pragma omp parallel for 
	for (int j = ipField.jStart; j <= ipField.jEnd; j++)
	{
		if (BC(ipField.iStart, j) == BC_NEUM)		fieldData(ipField.iStart, j) = fieldData(ipField.iStartI, j);
		if (BC(ipField.iStart, j) == BC_REFLECTION) fieldData(ipField.iStart, j) = -fieldData(ipField.iStartI, j);
		if (BC(ipField.iEnd, j) == BC_NEUM)			fieldData(ipField.iEnd, j)   = fieldData(ipField.iEndI, j);
		if (BC(ipField.iEnd, j) == BC_REFLECTION)	fieldData(ipField.iEnd, j)   = -fieldData(ipField.iEndI, j);
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
	double* val(ipVector.values);
	Array2D<int>& BC = ipField.BC;
	int iStart = ipField.iStart, iEnd = ipField.iEnd, jStart = ipField.jStart, jEnd = ipField.jEnd;
	int tempIndex = 0;
	Array2D<double>& dataArray = ipField.dataArray;

#pragma omp parallel for private(tempIndex)
	for (int i = iStart; i <= iEnd; i++)
	{
		for (int j = jStart; j <= jEnd; j++)
		{
			tempIndex = BC(i, j);

			if (tempIndex < 0) continue;

			dataArray(i, j) = val[tempIndex];
		}
	}
}

inline void FluidSolver2D::PlotVelocity()
{
	string str;

	U.Variable("U");
	V.Variable("V");
	
	str = string("quiver(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,2),axis equal,axis([Xp(1) Xp(end) Yp(1) Yp(end)]);");
	//str = str + string("hold on,streamline(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2,Xp(1:30:end),Yp(1:30:end)), hold off;axis equal tight;");
	str = str + string("hold on,streamslice(Xp,Yp,U(:,1:end-1)/2+U(:,2:end)/2,V(1:end-1,:)/2+V(2:end,:)/2), hold off;axis equal,axis([Xp(1) Xp(end) Yp(1) Yp(end)]);");
	if (isMultiPhase)
	{
		levelSet.phi.Variable("phi");
		str = str + string("hold on, contour(Xp,Yp,phi,[0 0],'r'), grid on,hold off;axis equal,axis([Xp(1) Xp(end) Yp(1) Yp(end)]);");
	}
	MATLAB.Command(str.c_str());
	str = string("title(['iteration : ', num2str(") + to_string(iteration) + string("),', time : ', num2str(") + to_string(totalT) + string("),', dt : ', num2str(") + to_string(dt) + string(")]);");
	MATLAB.Command(str.c_str());
	//MATLAB.Command("divU =U(:,2:end)-U(:,1:end-1),divV =V(2:end,:)-V(1:end-1,:);div=divU+divV;");
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

inline void FluidSolver2D::GenerateLinearSysteviscosityV2Order(CSR<double>& ipU_CSR, CSR<double>& ipV_CSR)
{

}

inline void FluidSolver2D::GenerateLinearSysteviscosityV2Order(VectorND<double>& vectorB)
{

}

inline void FluidSolver2D::GenerateLinearSystempPhi2Order(CSR<double>& ipCSR)
{

}

inline void FluidSolver2D::GenerateLinearSystempPhi2Order(VectorND<double>& vectorB)
{

}

inline void FluidSolver2D::SetLinearSystem(const int& iter)
{
	if (isMultiPhase || iter == 1)
	{
		if (ProjectionOrder == 1)
		{
			P_CSR = CSR<double>(Pressure.num_all_full_cells, Pressure.nnz);
			GenerateLinearSystemPressure(P_CSR);
			//P_CSR.RecoverCSR("Pmat");
		}
		else
		{
			GenerateLinearSysteviscosityV2Order(UCN_CSR, VCN_CSR);
			GenerateLinearSystempPhi2Order(PhiCN_CSR);

			//GenerateLinearSysteviscosityV(Fluid.UCNMatrix, Fluid.U.innerGrid, 1, Fluid.UCN_CSR);
			//GenerateLinearSysteviscosityV(Fluid.VCNMatrix, Fluid.V.innerGrid, 1, Fluid.VCN_CSR);

		}
	}

}