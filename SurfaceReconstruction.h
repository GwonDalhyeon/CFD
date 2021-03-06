#pragma once
#include "CommonDef.h"
#include "Grid2D.h"
#include "Field2D.h"
#include "CombineStructure.h"
#include "LevelSet2D.h"
#include "BregmanMethod.h"

double AdvectionMethod2D<double>::alpha;

template <class TT>
class SurfaceReconst
{
public:
	Grid2D grid;
	LS levelSet;
	FD distance;
	FD velocity;

	// Zero level set point.
	VectorND<VT> givenPoint;
	int givenPointNum;
	double dt;
	double cflCondition;

	int maxIteration;
	int reconstMaxIteration;
	int reinitialIter;
	int writeIter;
	
	double LpNorm;
	double distanceThreshold;
	double curvatureThreshold;

	double lambda;
	double mu;

	SurfaceReconst();
	~SurfaceReconst();

	void InitialCondition(int example);

	// surface Reconstruction using Variational Level Set Method.
	void SurfaceReconstructionSolver(int example);


	double Distance2Data(const int& i, const int& j);
	void ExactDistance();
	void SweepingDistance();



	void ComputeVelocity();
	double ComputeIntegralTerm();
	void LSPropagatingTVDRK3();


	bool StoppingCriterion();

	// Adaptive time step functions.
	double AdaptiveTimeStep();
	double AdaptiveTimeStep(const FD& velocity1);
	double AdaptiveTimeStep(const FD& velocity1, const FD& velocity2);



	//////////////////////////////////////////////////////////////////////////////
	//// Surface reconstruction : Split Bregman Method
	//// Notation : Geodesic Application of the Split Bregman Method ... Osher.

	int innerIStart;
	int innerJStart;
	int innerIEnd;
	int innerJEnd;
	int innerIRes;
	int innerJRes;

	void SurfaceReconstructionSplitBregman(const int & example);

	void InitialF(const int & example, FD& f);
	void GenerateLinearSystem(Array2D<double>& matrixA, const double & scaling);
	void GenerateLinearSystem(const FD& u, const FD& f, const FV& d, const FV& b, VectorND<double>& vectorB, const double & scaling);
	void OptimalU(const FD& f, const FV& d, const FV& b, const CSR<double>& csrA, FD& u);
	void OptimalD(const FD& ipField, const FV& b, FV& d);
	void OptimalB(const FD& ipField, const FV& d, FV& b);

private:

};

template <class TT>
SurfaceReconst<TT>::SurfaceReconst()
{
}

template <class TT>
SurfaceReconst<TT>::~SurfaceReconst()
{
}

template<class TT>
inline void SurfaceReconst<TT>::InitialCondition(int example)
{
	cflCondition = 0.5;
	reinitialIter = 20;

	if (example == 1)
	{

		cout << "******************************************************" << endl;
		cout << "         Surface reconstruction : One circles." << endl;
		cout << "              Good Initial Level Set." << endl;
		cout << "******************************************************" << endl;

		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LS(grid);
		distance = FD(grid);
		//distance.dataArray = 100;
		velocity = FD(grid);
		givenPointNum = 100;
		givenPoint = VectorND<VT>(givenPointNum);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		//dt = 5 * grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 500;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		lambda = 0.5;
		mu = 10e-5;

		for (int i = 0; i < givenPointNum; i++)
		{
			givenPoint(i) = 0.25*VT(cos(2 * PI*i / givenPointNum), sin(2 * PI*i / givenPointNum)) + 0.5 + grid.dx / 2;
		}

		//int i0, i1, j0, j1;
		//for (int k = 0; k < givenPointNum; k++)
		//{
		//	i0 = floor((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j0 = floor((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);
		//	i1 = ceil((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j1 = ceil((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);

		//	distance2Data(i0, j0);
		//	distance2Data(i0, j1);
		//	distance2Data(i1, j0);
		//	distance2Data(i1, j1);
		//}

		//distance2Data(); // distance initialization
		ExactDistance();

		// level set initialization
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if ((grid(i, j) - 0.5 - grid.dx / 2).magnitude() < 0.25)
				{
					levelSet.phi(i, j) = -distance(i, j) - grid.dx * 3;
				}
				else
				{
					levelSet.phi(i, j) = distance(i, j) - grid.dx * 3;
				}
			}
		}
		//for (int i = grid.iStart; i <= grid.iEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		levelSet.phi(i, j) = -(abs((grid(i, j) - 0.5).x) + abs((grid(i, j) - 0.5).y) - 0.28*sqrt(2));
		//	}
		//}

	}
	else if (example == 2)
	{
		cout << "******************************************************" << endl;
		cout << "         Surface reconstruction : One circles." << endl;
		cout << "              Good Initial Level Set." << endl;
		cout << "******************************************************" << endl;

		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LS(grid);
		distance = FD(grid);
		//distance.dataArray = 100;
		velocity = FD(grid);
		givenPointNum = 100;
		givenPoint = VectorND<VT>(givenPointNum);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		//dt = 5 * grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 500;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		lambda = 0.5;
		mu = 10e-5;
#pragma omp parallel for
		for (int i = 0; i < givenPointNum; i++)
		{
			givenPoint(i) = 0.25*VT(cos(2 * PI*i / givenPointNum), sin(2 * PI*i / givenPointNum)) + 0.5 + grid.dx / 2;
		}

		//int i0, i1, j0, j1;
		//for (int k = 0; k < givenPointNum; k++)
		//{
		//	i0 = floor((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j0 = floor((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);
		//	i1 = ceil((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j1 = ceil((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);

		//	distance2Data(i0, j0);
		//	distance2Data(i0, j1);
		//	distance2Data(i1, j0);
		//	distance2Data(i1, j1);
		//}

		//distance2Data(); // distance initialization
		ExactDistance();
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet.phi(i, j) = -(grid(i, j) - 0.5).magnitude() + 0.28;
			}
		}

		//for (int i = grid.iStart; i <= grid.iEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		levelSet.phi(i, j) = -(abs((grid(i, j) - 0.5).x) + abs((grid(i, j) - 0.5).y) - 0.28*sqrt(2));
		//	}
		//}
	}
	else if (example == 3)
	{
		cout << "******************************************************" << endl;
		cout << "         Surface reconstruction : Two circles." << endl;
		cout << "       Initial Level Set is far from real shape." << endl;
		cout << "******************************************************" << endl;

		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LS(grid);
		distance = FD(grid);
		//distance.dataArray = 100;
		velocity = FD(grid);
		givenPointNum = 400;
		givenPoint = VectorND<VT>(givenPointNum);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		//dt = 5 * grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 5000;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		lambda = 0.5;
		mu = 10e-5;

		VT point1(0.6, 0.4);
		VT point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*VT(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*VT(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		//int i0, i1, j0, j1;
		//for (int k = 0; k < givenPointNum; k++)
		//{
		//	i0 = floor((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j0 = floor((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);
		//	i1 = ceil((givenPoint(k)(0) - grid.xMin)*grid.oneOverdx);
		//	j1 = ceil((givenPoint(k)(1) - grid.yMin)*grid.oneOverdy);

		//	distance2Data(i0, j0);
		//	distance2Data(i0, j1);
		//	distance2Data(i1, j0);
		//	distance2Data(i1, j1);
		//}

		//distance2Data(); // distance initialization
		ExactDistance();

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet.phi(i, j) = -(grid(i, j) - 0.5).magnitude() + 0.28;
			}
		}

		//for (int i = grid.iStart; i <= grid.iEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		levelSet.phi(i, j) = -(abs((grid(i, j) - 0.5).x) + abs((grid(i, j) - 0.5).y) - 0.28*sqrt(2));
		//	}
		//}
	}
	else if (example == 4)
	{
		cout << "******************************************************" << endl;
		cout << "   Surface reconstruction : Two circles with outlier." << endl;
		cout << "       Initial Level Set is far from real shape." << endl;
		cout << "******************************************************" << endl;


		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LS(grid);
		distance = FD(grid);
		velocity = FD(grid);
		givenPointNum = 400;
		int outlier = 20;
		givenPoint = VectorND<VT>(givenPointNum + outlier);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		reconstMaxIteration = 2000;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		lambda = 0.5;
		mu = 10e-5;

		VT point1(0.6, 0.4);
		VT point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*VT(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*VT(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		srand(time(NULL));
		for (int i = 0; i < outlier; i++)
		{
			VT tempVector(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
			if (((tempVector - 0.5) / 3).magnitude()<0.28)
			{
				givenPoint(i + givenPointNum) = (tempVector - 0.5) / 3 + 0.5;
			}
			else
			{
				i--;
			}

		}
		givenPointNum = givenPointNum + outlier;

		ExactDistance();

#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet.phi(i, j) = -(grid(i, j) - 0.5).magnitude() + 0.28;
			}
		}
	}
	else if (example == 5)
	{
		cout << "******************************************************" << endl;
		cout << "  Surface reconstruction : Two circles with outlier." << endl;
		cout << "              Good Initial Level Set." << endl;
		cout << "******************************************************" << endl;


		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LS(grid);
		distance = FD(grid);
		velocity = FD(grid);
		givenPointNum = 400;
		int outlier = 20;
		givenPoint = VectorND<VT>(givenPointNum + outlier);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		reconstMaxIteration = 2000;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		lambda = 0.5;
		mu = 10e-5;

		VT point1(0.6, 0.4);
		VT point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*VT(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*VT(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		srand(time(NULL));
		for (int i = 0; i < outlier; i++)
		{
			VT tempVector(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
			if (((tempVector - 0.5) / 3).magnitude()<0.28)
			{
				givenPoint(i + givenPointNum) = (tempVector - 0.5) / 3 + 0.5;
			}
			else
			{
				i--;
			}

		}
		givenPointNum = givenPointNum + outlier;

		ExactDistance();

		double initialEpsilon = 0.02;
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if ((grid(i, j) - point1 - grid.dx / 2).magnitude()<0.1 || (grid(i, j) - point2 - grid.dx / 2).magnitude()<0.1)
				{
					levelSet(i, j) = distance(i, j) + initialEpsilon;
				}
				else
				{
					levelSet(i, j) = -distance(i, j) + initialEpsilon;
				}
			}
		}
	}
	else if (example == 6)
	{
		cout << "******************************************************" << endl;
		cout << "   Surface reconstruction : Two circles with outlier." << endl;
		cout << "       Initial Level Set is far from real shape." << endl;
		cout << "******************************************************" << endl;

		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LS(grid);
		distance = FD(grid);
		velocity = FD(grid);
		givenPointNum = 400;
		int outlier = 20;
		givenPoint = VectorND<VT>(givenPointNum + outlier);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		reconstMaxIteration = 2000;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		lambda = 0.5;
		mu = 10e-5;

		VT point1(0.6, 0.4);
		VT point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*VT(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*VT(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		srand(time(NULL));
		for (int i = 0; i < outlier; i++)
		{
			VT tempVector(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
			if (((tempVector - 0.5) / 3).magnitude()<0.28)
			{
				givenPoint(i + givenPointNum) = (tempVector - 0.5) / 3 + 0.5;
			}
			else
			{
				i--;
			}

		}
		givenPointNum = givenPointNum + outlier;

		ExactDistance();
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet.phi(i, j) = -(grid(i, j) - 0.5).magnitude() + 0.28;
			}
		}
	}


}

template<class TT>
inline void SurfaceReconst<TT>::SurfaceReconstructionSolver(int example)
{
	bool writeFile = false;
	string str;
	const char* cmd;

	InitialCondition(example);
	
	if (writeFile)
	{
		str = "distance";
		distance.WriteFile("distance");
		str = "pointData";
		givenPoint.WriteFile("pointData");
	}
	
	grid.Variable();
	distance.Variable("distance");
	VecND2DVariable("pointData", givenPoint);
	levelSet.phi.Variable("phi0");

	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1/2 1])");

	// MATLAB command : start
	velocity.Variable("velocity");
	//MATLAB.Command("subplot(1, 2, 1)");
	MATLAB.Command("plot(pointData(:,1), pointData(:,2),'ro');axis([0 1 0 1]);grid on, axis equal");
	MATLAB.Command("hold on");
	MATLAB.Command("contour(X,Y,phi0,[0 0],'color','r');");
	MATLAB.Command("contour(X,Y,phi,[0 0]);");
	MATLAB.Command("axis([0 1 0 1]);");
	MATLAB.Command("hold off");
	str = string("title(['iteration : ', num2str(") + to_string(0) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);
	MATLAB.Command("drawnow;");
	// MATLAB command : end

	//// MATLAB command : start
	//levelSet.phi.Variable("phi");
	//MATLAB.Command("subplot(1, 2, 2)");
	//MATLAB.Command("surf(X,Y,phi), axis equal");
	//str = string("title(['velocity : ', num2str(") + to_string(0) + string(")]);");
	//cmd = str.c_str();
	//MATLAB.Command(cmd);
	//// MATLAB command : end

	int i = 0;
	//while (!StoppingCriterion() && i < reconstMaxIteration)
	while (i < reconstMaxIteration)
	{
		i++;
		cout << "Surface reconstruction : " << i << endl;

		ComputeVelocity();
		velocity.Variable("velocity");

		// MATLAB command : start
		//MATLAB.Command("subplot(1, 2, 1)");
		levelSet.phi.Variable("phi");
		MATLAB.Command("plot(pointData(:,1), pointData(:,2),'ro');axis([0 1 0 1]);grid on, axis equal");
		MATLAB.Command("hold on");
		MATLAB.Command("contour(X,Y,phi0,[0 0],'color','r');");
		MATLAB.Command("contour(X,Y,phi,[0 0]);");
		MATLAB.Command("axis([0 1 0 1]);");
		MATLAB.Command("hold off");
		str = string("title(['iteration : ', num2str(") + to_string(i) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
		MATLAB.Command("drawnow;");
		// MATLAB command : end

		dt = AdaptiveTimeStep(velocity);
		LSPropagatingTVDRK3();

		//// MATLAB command : start
		//MATLAB.Command("subplot(1, 2, 2)");
		//MATLAB.Command("surf(X,Y,phi), axis equal");
		//str = string("title(['velocity : ', num2str(") + to_string(i) + string(")]);");
		//cmd = str.c_str();
		//MATLAB.Command(cmd);
		//// MATLAB command : end

		if (writeFile)
		{
			str = "velocity" + to_string(i);
			velocity.WriteFile(str);
		}
		
		dt = AdaptiveTimeStep();
		for (int j = 0; j < reinitialIter; j++)
		{
			cout << "Reinitialization : " << i << "-" << j + 1 << endl;
			AdvectionMethod2D<double>::LSPropagatingTVDRK3(levelSet, dt);
		}

		if (i%writeIter && writeFile)
		{
			str = "phi" + to_string(i);
			levelSet.phi.WriteFile(str);
		}
	}
}


template<class TT>
inline double SurfaceReconst<TT>::Distance2Data(const int & i, const int & j)
{
	double distance = 100;
	double tempDist;

	for (int k = 0; k < givenPointNum; k++)
	{
		tempDist = (givenPoint(k) - grid(i, j)).magnitude();
		if (distance > tempDist)
		{
			distance = tempDist;
		}
	}

	return distance;
}

template<class TT>
inline void SurfaceReconst<TT>::ExactDistance()
{
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			distance(i, j) = Distance2Data(i, j);
		}
	}
}

template<class TT>
inline void SurfaceReconst<TT>::SweepingDistance()
{
	double xMin, yMin;
	double h = min(grid.dx, grid.dy);

#pragma omp parallel for private(xMin, yMin)
	for (int i = distance.iStart; i <= distance.iEnd; i++)
	{
		for (int j = distance.jStart; j <= distance.jEnd; j++)
		{
			xMin = min(distance(max(i - 1, distance.iStart), j), distance(min(i + 1, distance.iEnd), j));
			yMin = min(distance(i, max(j - 1, distance.jStart)), distance(i, min(j + 1, distance.jEnd)));
			//#pragma omp parallel for
			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				distance(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}

#pragma omp parallel for private(xMin, yMin)
	for (int i = distance.iStart; i <= distance.iEnd; i++)
	{
		for (int j = distance.jEnd; j >= distance.jStart; j--)
		{
			xMin = min(distance(max(i - 1, distance.iStart), j), distance(min(i + 1, distance.iEnd), j));
			yMin = min(distance(i, max(j - 1, distance.jStart)), distance(i, min(j + 1, distance.jEnd)));

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				distance(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}


#pragma omp parallel for private(xMin, yMin)
	for (int i = distance.iEnd; i >= distance.iStart; i--)
	{
		for (int j = distance.jStart; j <= distance.jEnd; j++)
		{
			xMin = min(distance(max(i - 1, distance.iStart), j), distance(min(i + 1, distance.iEnd), j));
			yMin = min(distance(i, max(j - 1, distance.jStart)), distance(i, min(j + 1, distance.jEnd)));

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				distance(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}

#pragma omp parallel for private(xMin, yMin)
	for (int i = distance.iEnd; i >= distance.iStart; i--)
	{
		for (int j = distance.jEnd; j >= distance.jStart; j--)
		{
			xMin = min(distance(max(i - 1, distance.iStart), j), distance(min(i + 1, distance.iEnd), j));
			yMin = min(distance(i, max(j - 1, distance.jStart)), distance(i, min(j + 1, distance.jEnd)));

			if (abs(xMin - yMin) >= h)
			{
				distance(i, j) = min(xMin, yMin) + h;
			}
			else
			{
				distance(i, j) = (xMin + yMin + sqrt(2 * h*h - (xMin - yMin)*(xMin - yMin))) / 2.0;
			}
		}
	}
}



template<class TT>
inline void SurfaceReconst<TT>::ComputeVelocity()
{
	double integralTerm = ComputeIntegralTerm();
	double lastTerm;
	Vector2D <double> gradPhi;

	double tempDist;
	double tempCurvature;

	levelSet.ComputeMeanCurvature();

#pragma omp parallel for private(lastTerm, tempDist,tempCurvature, gradPhi)
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			velocity(i, j) = 0;
			gradPhi = levelSet.gradient(i, j);
			if (gradPhi.magnitude() < 3 * DBL_EPSILON)
			{
				lastTerm = 0;
			}
			else
			{
				//lastTerm = DotProduct(distance.gradient(i, j), gradPhi / (gradPhi.magnitude() + DBL_EPSILON)) + 1.0 / LpNorm*distance(i, j)*levelSet.meanCurvature(i, j);
				tempDist = min(distance(i, j), distanceThreshold);
				tempCurvature = AdvectionMethod2D<double>::sign(levelSet.meanCurvature(i, j))* min(abs(levelSet.meanCurvature(i, j)), curvatureThreshold);
				lastTerm = DotProduct(distance.Gradient(i, j), gradPhi / (gradPhi.magnitude() + DBL_EPSILON)) + 1.0 / LpNorm*tempDist*tempCurvature;
			}

			//velocity(i, j) = gradPhi.magnitude()*integralTerm*pow(distance(i, j), LpNorm - 1)*lastTerm;
			//velocity(i, j) = gradPhi.magnitude()*integralTerm*pow(min(distance(i, j), distanceThreshold), LpNorm - 1)*lastTerm;
			velocity(i, j) = gradPhi.magnitude()*integralTerm*lastTerm;
		}
	}
}

template<class TT>
inline double SurfaceReconst<TT>::ComputeIntegralTerm()
{
	double sum = 0.0;
#pragma omp parallel for reduction(+:sum)
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			sum += pow(distance(i, j), LpNorm)*AdvectionMethod2D<double>::DeltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude()*grid.dx*grid.dy;
		}
	}
	return pow(sum, 1.0 / LpNorm - 1.0);
}

template<class TT>
inline void SurfaceReconst<TT>::LSPropagatingTVDRK3()
{
	LS originLevelSet = levelSet;

	Array2D<double>& k1 = levelSet.phi.K1;
	Array2D<double>& k2 = levelSet.phi.K2;
	Array2D<double>& k3 = levelSet.phi.K3;

	Array2D<double>& wenoXMinus = levelSet.phi.dfdxM;
	Array2D<double>& wenoXPlus = levelSet.phi.dfdxP;
	Array2D<double>& wenoYMinus = levelSet.phi.dfdyM;
	Array2D<double>& wenoYPlus = levelSet.phi.dfdyP;

	AdvectionMethod2D<double>::WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = dt*velocity(i, j);
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}
	AdvectionMethod2D<double>::WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = dt*velocity(i, j);
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}
	AdvectionMethod2D<double>::WENO5thDerivation(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k3(i, j) = dt*velocity(i, j);
			levelSet(i, j) = 1.0 / 3.0*originLevelSet(i, j) + 2.0 / 3.0*levelSet(i, j) + 2.0 / 3.0*k3(i, j);
		}
	}
}





template<class TT>
inline bool SurfaceReconst<TT>::StoppingCriterion()
{
	bool criterion = true;

	for (int i = 0; i < givenPointNum; i++)
	{
		if (abs(levelSet.interpolation(givenPoint(i)))>grid.dx + grid.dy)
		{
			criterion = false;
			return criterion;
		}
	}

	return criterion;
}

template<class TT>
inline double SurfaceReconst<TT>::AdaptiveTimeStep()
{
	return cflCondition*max(grid.dx, grid.dy);
}

template<class TT>
inline double SurfaceReconst<TT>::AdaptiveTimeStep(const FD& velocity1)
{
	double maxVel1 = 0;

	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			if (abs(velocity1(i, j)) > maxVel1)
			{
				maxVel1 = abs(velocity1(i, j));
			}
		}
	}
	return cflCondition*max(grid.dx, grid.dy) / maxVel1;
}

template<class TT>
inline double SurfaceReconst<TT>::AdaptiveTimeStep(const FD& velocity1, const FD& velocity2)
{
	double maxVel1 = 0;
	double maxVel2 = 0;

	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			if (abs(velocity1(i, j)) > maxVel1)
			{
				maxVel1 = abs(velocity1(i, j));
			}
			if (abs(velocity2(i, j)) > maxVel2)
			{
				maxVel2 = abs(velocity2(i, j));
			}
		}
	}
	return cflCondition*(grid.dx / maxVel1 + grid.dy / maxVel2);
}

template<class TT>
inline void SurfaceReconst<TT>::SurfaceReconstructionSplitBregman(const int & example)
{

	bool writeFile = false;
	string str;
	const char* cmd;

	InitialCondition(example);

	FV b(grid);
	FV d(grid);
	FD f(grid);

	InitialF(example, f);


	f.Variable("aaaa");

	if (writeFile)
	{
		str = "distance";
		distance.WriteFile("distance");
		str = "pointData";
		givenPoint.WriteFile("pointData");
	}


	grid.Variable();
	distance.Variable("distance");
	VecND2DVariable("pointData", givenPoint);
	levelSet.phi.Variable("phi0");

	if (writeFile)
	{
		givenPoint.WriteFile("pointData");
		distance.WriteFile("distance");
		levelSet.phi.WriteFile("phi0");
	}

	innerIStart = grid.iStart + 1;
	innerJStart = grid.jStart + 1;
	innerIEnd = grid.iEnd - 1;
	innerJEnd = grid.jEnd - 1;
	innerIRes = innerIEnd - innerIStart + 1;
	innerJRes = innerJEnd - innerJStart + 1;;

	

	Array2D<double> matrixA(innerIStart, innerIRes*innerJRes, innerJStart, innerIRes*innerJRes);
	GenerateLinearSystem(matrixA, 1);
	matrixA.Variable("A");

	cout << "Start CSR." << endl;
	CSR<double> csrA(matrixA);
	cout << "End CSR." << endl;

	
	int maxIter = 20;
	for (int i = 1; i < 2; i++)
	{
		OptimalU(distance, d, b, csrA, levelSet.phi);
		levelSet.phi.Variable("u");
		levelSet.phi.Gradient();
		ArrayVec2DVariable("grad", levelSet.phi.gradient);

		OptimalD(levelSet.phi, b, d);

		OptimalB(levelSet.phi, d, b);

		if (WriteFile)
		{
			str = "phi" + to_string(i);
			levelSet.phi.WriteFile(str);
		}
	}

}

template<class TT>
inline void SurfaceReconst<TT>::InitialF(const int & example, FD& f)
{
	if (example == 1)
	{
#pragma omp parallel for
		for (int i = f.iStart; i <= f.iEnd; i++)
		{
			for (int j = f.jStart; j <= f.jEnd; j++)
			{
				if ((f.grid(i,j)- 0.5 - grid.dx / 2).magnitude()<0.25)
				{
					f(i, j) = 1;
				}
				else
				{
					f(i, j) = 0;
				}
			}
		}
	}
}

template<class TT>
inline void SurfaceReconst<TT>::GenerateLinearSystem(Array2D<double>& matrixA, const double & scaling)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int index, leftIndex, rightIndex, bottomIndex, topIndex;

#pragma omp parallel for private(index, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			leftIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart - 1) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			rightIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart + 1) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart)*innerIRes;
			bottomIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart - 1)*innerIRes;
			topIndex = (i - innerIStart)*innerIRes*innerJRes + (i - innerIStart) + (j - innerJStart)*innerIRes*innerIRes*innerJRes + (j - innerJStart + 1)*innerIRes;

			matrixA(index) = scaling * (mu + 4 * lambda * grid.oneOverdx2);

			if (i>innerIStart)
			{
				matrixA(leftIndex) = scaling * (-lambda * grid.oneOverdx2);
			}
			if (i<innerIEnd)
			{
				matrixA(rightIndex) = scaling * (-lambda * grid.oneOverdx2);
			}
			if (j>innerJStart)
			{
				matrixA(bottomIndex) = scaling * (-lambda * grid.oneOverdy2);
			}
			if (j<innerJEnd)
			{
				matrixA(topIndex) = scaling * (-lambda * grid.oneOverdy2);
			}
		}
	}
	cout << "End Generate Linear System : matrix A" << endl;
}

template<class TT>
inline void SurfaceReconst<TT>::GenerateLinearSystem(const FD& u, const FD& f, const FV& d, const FV& b, VectorND<double>& vectorB, const double & scaling)
{
	int index;

	FV temp(grid);
	temp.dataArray = b.dataArray - d.dataArray;
	FD div = FD::Divegence(temp);

#pragma omp parallel for private(index)
	for (int i = innerIStart; i <= innerIEnd; i++)
	{
		for (int j = innerJStart; j <= innerJEnd; j++)
		{
			index = (i - innerIStart) + (j - innerJStart)*innerIRes;

			vectorB(index) = scaling * (mu*f(i, j) + lambda*div(i, j));

			if (i == innerIStart)
			{
				vectorB(index) += scaling * lambda* u(i - 1, j) * grid.oneOverdx2;
			}
			if (i == innerIEnd)
			{
				vectorB(index) += scaling * lambda* u(i + 1, j) * grid.oneOverdx2;
			}
			if (j == innerJStart)
			{
				vectorB(index) += scaling * lambda* u(i, j - 1) * grid.oneOverdy2;
			}
			if (j == innerJEnd)
			{
				vectorB(index) += scaling * lambda* u(i, j + 1) * grid.oneOverdy2;
			}
		}
	}
}

template<class TT>
inline void SurfaceReconst<TT>::OptimalU(const FD& f, const FV& d, const FV& b, const CSR<double>& csrA, FD& u)
{
	cout << "Start Optimal U" << endl;


	VectorND<double> vectorB(innerIRes*innerJRes);
	GenerateLinearSystem(u, f, d, b, vectorB, 1);
	vectorB.Variable("b");

	VectorND<double> uStar(innerIRes*innerJRes);
	uStar = CGSolver::SolverCSR(csrA, vectorB, grid.dx*grid.dy);
	uStar.Variable("uStar");

	int idx;

#pragma omp parallel for private(idx)
	for (int i = u.iStart + 1; i <= u.iEnd - 1; i++)
	{
		for (int j = u.jStart + 1; j <= u.jEnd - 1; j++)
		{
			idx = (i - innerIStart) + (j - innerJStart)*innerIRes;
			u(i, j) = uStar(idx);
		}
	}

	cout << "End Optimal U" << endl;
	cout << endl;
}

template<class TT>
inline void SurfaceReconst<TT>::OptimalD(const FD & ipField, const FV & b, FV & d)
{
#pragma omp parallel for
	for (int i = d.iStart; i <= d.iEnd; i++)
	{
		for (int j = d.jStart; j <= d.jEnd; j++)
		{
			d(i, j) = BregmanMethod<double>::Shrink(ipField.gradient(i, j) + b(i, j), lambda);
		}
	}
}

template<class TT>
inline void SurfaceReconst<TT>::OptimalB(const FD & ipField, const FV & d, FV & b)
{
#pragma omp parallel for
	for (int i = b.iStart; i <= b.iEnd; i++)
	{
		for (int j = b.jStart; j <= b.jEnd; j++)
		{
			b(i, j) = b(i, j) + ipField.gradient(i, j) - d(i, j);
		}
	}
}
