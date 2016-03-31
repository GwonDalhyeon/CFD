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
	LevelSet2D levelSet;
	Field2D<double> distance;
	Field2D<double> velocity;

	// Zero level set point.
	VectorND<Vector2D<double>> givenPoint;
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

	

	// surface Reconstruction using Variational Level Set Method.
	void InitialCondition(int example);
	void SurfaceReconstructionSolver(int example);


	void Distance2Data();
	double Distance2Data(const int& i, const int& j);
	void ExactDistance();
	void SweepingDistance();



	void ComputeVelocity();
	double ComputeIntegralTerm();
	void LevelSetPropagatingTVDRK3();


	bool StoppingCriterion();

	// Adaptive time step functions.
	double AdaptiveTimeStep();
	double AdaptiveTimeStep(const Field2D<double>& velocity1);
	double AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2);

	//////////////////////////////////////////////////////////////////////////////
	//// Surface reconstruction : Split Bregman Method
	//// Notation : Geodesic Application of the Split Bregman Method ... Osher.
	void SurfaceReconstructionSplitBregman(const int & example, const bool & propa, const bool & reinitial, const bool & surfReconst);
	void GenerateLinearSystem(Array2D<double>& matrixA);
	void GenerateLinearSystem(const Field2D<double>& u, const Field2D<double>& f, const Field2D<Vector2D<double>>& d, const Field2D<Vector2D<double>>& b, VectorND<double>& vectorB);
	void OptimalU(const Field2D<double>& f, const Field2D<Vector2D<double>>& d, const Field2D<Vector2D<double>>& b, const CSR<double>& csrA, Field2D<double>& u);
	void OptimalD(const Field2D<Vector2D<double>>& gradientU, const Field2D<Vector2D<double>>& b, Field2D<Vector2D<double>>& d);
	void OptimalB(const Field2D<Vector2D<double>>& gradientU, const Field2D<Vector2D<double>>& d, Field2D<Vector2D<double>>& b);

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
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		//distance.dataArray = 100;
		velocity = Field2D<double>(grid);
		givenPointNum = 100;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		//dt = 5 * grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 500;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		for (int i = 0; i < givenPointNum; i++)
		{
			givenPoint(i) = 0.25*Vector2D<double>(cos(2 * PI*i / givenPointNum), sin(2 * PI*i / givenPointNum)) + 0.5 + grid.dx / 2;
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
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		//distance.dataArray = 100;
		velocity = Field2D<double>(grid);
		givenPointNum = 100;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		//dt = 5 * grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 500;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		for (int i = 0; i < givenPointNum; i++)
		{
			givenPoint(i) = 0.25*Vector2D<double>(cos(2 * PI*i / givenPointNum), sin(2 * PI*i / givenPointNum)) + 0.5 + grid.dx / 2;
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
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		//distance.dataArray = 100;
		velocity = Field2D<double>(grid);
		givenPointNum = 400;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		//dt = 5 * grid.dx*grid.dy / 2.0;
		reconstMaxIteration = 5000;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		Vector2D<double> point1(0.6, 0.4);
		Vector2D<double> point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
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
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		velocity = Field2D<double>(grid);
		givenPointNum = 400;
		int outlier = 20;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum + outlier);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		reconstMaxIteration = 2000;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		Vector2D<double> point1(0.6, 0.4);
		Vector2D<double> point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		srand(time(NULL));
		for (int i = 0; i < outlier; i++)
		{
			Vector2D<double> tempVector(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
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
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		velocity = Field2D<double>(grid);
		givenPointNum = 400;
		int outlier = 20;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum + outlier);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		reconstMaxIteration = 2000;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		Vector2D<double> point1(0.6, 0.4);
		Vector2D<double> point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		srand(time(NULL));
		for (int i = 0; i < outlier; i++)
		{
			Vector2D<double> tempVector(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
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
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		velocity = Field2D<double>(grid);
		givenPointNum = 400;
		int outlier = 20;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum + outlier);
		AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy);
		reconstMaxIteration = 2000;
		writeIter = 10;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		lambda = 0.5;
		mu = 10e-5;

		Vector2D<double> point1(0.6, 0.4);
		Vector2D<double> point2(0.4, 0.6);

#pragma omp parallel for
		for (int i = 0; i < givenPointNum / 2; i++)
		{
			givenPoint(i) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point1 + grid.dx / 2;
			givenPoint(i + givenPointNum / 2) = 0.1*Vector2D<double>(cos(2 * PI*i / givenPointNum * 2), sin(2 * PI*i / givenPointNum * 2)) + point2 + grid.dx / 2;
		}

		srand(time(NULL));
		for (int i = 0; i < outlier; i++)
		{
			Vector2D<double> tempVector(double(rand()) / double(RAND_MAX), double(rand()) / double(RAND_MAX));
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
		distance.WriteFile(str.c_str);
		str = "pointData";
		givenPoint.WriteFile(str.c_str);
	}
	
	
	grid.Variable();
	distance.Variable("distance");
	VecND2DVariable("pointData", givenPoint);
	levelSet.phi.Variable("phi0");


	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");

	// MATLAB command : start
	velocity.Variable("velocity");
	MATLAB.Command("subplot(1, 2, 1)");
	MATLAB.Command("plot(pointData(:,1), pointData(:,2),'ro');axis([0 1 0 1]);grid on");
	MATLAB.Command("hold on");
	MATLAB.Command("contour(X,Y,phi0,[0 0],'color','r');");
	MATLAB.Command("contour(X,Y,phi,[0 0]);");
	MATLAB.Command("axis([0 1 0 1]);");
	MATLAB.Command("hold off");
	str = string("title(['velocity : ', num2str(") + to_string(0) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);
	MATLAB.Command("drawnow;");
	// MATLAB command : end

	dt = AdaptiveTimeStep(velocity);
	LevelSetPropagatingTVDRK3();

	// MATLAB command : start
	levelSet.phi.Variable("phi");
	MATLAB.Command("subplot(1, 2, 2)");
	MATLAB.Command("surf(X,Y,phi)");
	str = string("title(['velocity : ', num2str(") + to_string(0) + string(")]);");
	cmd = str.c_str();
	MATLAB.Command(cmd);
	// MATLAB command : end


	int i = 0;
	//while (i < reconstMaxIteration)
	while (!StoppingCriterion() && i < reconstMaxIteration)
	{
		i++;
		cout << "Surface reconstruction : " << i << endl;

		ComputeVelocity();
		velocity.Variable("velocity");

		// MATLAB command : start
		MATLAB.Command("subplot(1, 2, 1)");
		MATLAB.Command("plot(pointData(:,1), pointData(:,2),'ro');axis([0 1 0 1]);grid on");
		MATLAB.Command("hold on");
		MATLAB.Command("contour(X,Y,phi0,[0 0],'color','r');");
		MATLAB.Command("contour(X,Y,phi,[0 0]);");
		MATLAB.Command("axis([0 1 0 1]);");
		MATLAB.Command("hold off");
		str = string("title(['velocity : ', num2str(") + to_string(i) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
		MATLAB.Command("drawnow;");
		// MATLAB command : end

		dt = AdaptiveTimeStep(velocity);
		LevelSetPropagatingTVDRK3();

		// MATLAB command : start
		levelSet.phi.Variable("phi");
		MATLAB.Command("subplot(1, 2, 2)");
		MATLAB.Command("surf(X,Y,phi)");
		str = string("title(['velocity : ', num2str(") + to_string(i) + string(")]);");
		cmd = str.c_str();
		MATLAB.Command(cmd);
		// MATLAB command : end

		if (writeFile)
		{
			str = "velocity" + to_string(iter);
			reconstructionVelocity.WriteFile(str.c_str);
		}
		
		dt = AdaptiveTimeStep();
		for (int j = 0; j < reinitialIter; j++)
		{
			cout << "Reinitialization : " << i << "-" << j + 1 << endl;
			AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
		}

		if (i%writeIter && writeFile)
		{
			str = "phi" + to_string(i);
			levelSet.phi.WriteFile(str.c_str);
		}
		

		cout << endl;
	}
	OutputResult(i);
}

template<class TT>
inline void SurfaceReconst<TT>::Distance2Data()
{

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

	levelSet.computeMeanCurvature();

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
				//lastTerm = dotProduct(distance.gradient(i, j), gradPhi / (gradPhi.magnitude() + DBL_EPSILON)) + 1.0 / LpNorm*distance(i, j)*levelSet.meanCurvature(i, j);
				tempDist = min(distance(i, j), distanceThreshold);
				tempCurvature = AdvectionMethod2D<double>::sign(levelSet.meanCurvature(i, j))* min(abs(levelSet.meanCurvature(i, j)), curvatureThreshold);
				lastTerm = dotProduct(distance.gradient(i, j), gradPhi / (gradPhi.magnitude() + DBL_EPSILON)) + 1.0 / LpNorm*tempDist*tempCurvature;
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
			sum += pow(distance(i, j), LpNorm)*AdvectionMethod2D<double>::deltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude()*grid.dx*grid.dy;
		}
	}
	return pow(sum, 1.0 / LpNorm - 1.0);
}

template<class TT>
inline void SurfaceReconst<TT>::LevelSetPropagatingTVDRK3()
{
	LevelSet2D originLevelSet = levelSet;
	LevelSet2D tempLevelSet(originLevelSet.grid);

	Field2D<double> k1(levelSet.grid);
	Field2D<double> k2(levelSet.grid);
	Field2D<double> k3(levelSet.grid);

	Field2D<double> wenoXMinus(levelSet.grid);
	Field2D<double> wenoXPlus(levelSet.grid);
	Field2D<double> wenoYMinus(levelSet.grid);
	Field2D<double> wenoYPlus(levelSet.grid);

	AdvectionMethod2D<double>::WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k1(i, j) = dt*velocity(i, j);
			levelSet(i, j) = originLevelSet(i, j) + k1(i, j);
		}
	}

	AdvectionMethod2D<double>::WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
#pragma omp parallel for
	for (int i = levelSet.grid.iStart; i <= levelSet.grid.iEnd; i++)
	{
		for (int j = levelSet.grid.jStart; j <= levelSet.grid.jEnd; j++)
		{
			k2(i, j) = dt*velocity(i, j);
			levelSet(i, j) = 3.0 / 4.0*originLevelSet(i, j) + 1.0 / 4.0*levelSet(i, j) + 1.0 / 4.0*k2(i, j);
		}
	}

	AdvectionMethod2D<double>::WENO5th(levelSet.phi, wenoXMinus, wenoXPlus, wenoYMinus, wenoYPlus);
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
inline double SurfaceReconst<TT>::AdaptiveTimeStep(const Field2D<double>& velocity1)
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
inline double SurfaceReconst<TT>::AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2)
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
inline void SurfaceReconst<TT>::SurfaceReconstructionSplitBregman(const int & example, const bool & propa, const bool & reinitial, const bool & surfReconst)
{
	InitialCondition(example);

	Field2D<Vector2D<double>> b = Field2D<Vector2D<double>>(grid);
	Field2D<Vector2D<double>> d = Field2D<Vector2D<double>>(grid);
	Field2D<Vector2D<double>> gradientU = Field2D<Vector2D<double>>(grid);;
	Field2D<double> f = Field2D<double>(grid);

	distance = distance;
#pragma omp parallel for
	for (int i = distance.iStart; i <= distance.iEnd; i++)
	{
		for (int j = distance.jStart; j <= distance.jEnd; j++)
		{
			distance(i, j) = -distance(i, j);
		}
	}

	Array2D<double> matrixA = Array2D<double>(1, (grid.iRes - 2)*(grid.jRes - 2), 1, (grid.iRes - 2)*(grid.jRes - 2));
	GenerateLinearSystem(matrixA);
	cout << "Start CSR." << endl;
	CSR<double> csrA(matrixA);
	cout << "End CSR." << endl;

	string fileName;

	fileName = "pointData" + to_string(0);
	givenPoint.WriteFile(fileName);
	fileName = "distance" + to_string(0);
	distance.WriteFile(fileName);
	fileName = "phi" + to_string(0);
	levelSet.phi.WriteFile(fileName);



	int maxIter = 20;
	for (int i = 1; i < maxIter; i++)
	{
		OptimalU(distance, d, b, csrA, levelSet.phi);
		gradientU = Field2D<Vector2D<double>>::Gradient(levelSet.phi);
		OptimalD(gradientU, b, d);
		OptimalB(gradientU, d, b);

		fileName = "phi" + to_string(i);
		levelSet.phi.WriteFile(fileName);

	}

}

template<class TT>
inline void SurfaceReconst<TT>::GenerateLinearSystem(Array2D<double>& matrixA)
{
	cout << "Start Generate Linear System : matrix A" << endl;
	int index, leftIndex, rightIndex, bottomIndex, topIndex;
	int innerIRes = grid.iRes - 2;
	int innerJRes = grid.jRes - 2;

#pragma omp parallel for private(index, leftIndex, rightIndex, bottomIndex, topIndex)
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			index = (i - 1)*innerIRes*innerJRes + (i - 1) + (j - 1)*innerIRes*(innerIRes*innerJRes + 1);
			leftIndex = (i - 1)*innerIRes*innerJRes + (i - 1 - 1) + (j - 1)*innerIRes*(innerIRes*innerJRes + 1);
			rightIndex = (i - 1)*innerIRes*innerJRes + i + (j - 1)*innerIRes*(innerIRes*innerJRes + 1);
			bottomIndex = (i - 1)*innerIRes*innerJRes + (i - 1) + (j - 1)*innerIRes*innerIRes*innerJRes + (j - 1 - 1)*innerIRes;
			topIndex = (i - 1)*innerIRes*innerJRes + (i - 1) + (j - 1)*innerIRes*innerIRes*innerJRes + (j)*innerIRes;

			matrixA(index) = -(mu - 4 * lambda / grid.dx2);

			if (i>grid.iStart + 1)
			{
				matrixA(leftIndex) = lambda / grid.dx2;
			}
			if (i<grid.iEnd - 1)
			{
				matrixA(rightIndex) = lambda / grid.dx2;
			}
			if (j>grid.jStart + 1)
			{
				matrixA(bottomIndex) = lambda / grid.dx2;
			}
			if (j<grid.jEnd - 1)
			{
				matrixA(topIndex) = lambda / grid.dx2;
			}
		}
	}
	cout << "End Generate Linear System : matrix A" << endl;
}

template<class TT>
inline void SurfaceReconst<TT>::GenerateLinearSystem(const Field2D<double>& u, const Field2D<double>& f, const Field2D<Vector2D<double>>& d, const Field2D<Vector2D<double>>& b, VectorND<double>& vectorB)
{
	int index;
	int innerIRes = grid.iRes - 2;

	Field2D<Vector2D<double>> temp(grid);
	temp.dataArray = b.dataArray - d.dataArray;
	Field2D<double> div = Field2D<double>::Divegence(temp);

#pragma omp parallel for
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			index = (i - 1) + (j - 1)*innerIRes;

			vectorB(index) = -(mu*f(i, j) + lambda*div(i, j));

			if (i == grid.iStart + 1)
			{
				vectorB(index) += -lambda* u(i - 1, j) / grid.dx2;
			}
			if (i == grid.iEnd - 1)
			{
				vectorB(index) += -lambda* u(i + 1, j) / grid.dx2;
			}
			if (j == grid.jStart + 1)
			{
				vectorB(index) += -lambda* u(i, j - 1) / grid.dy2;
			}
			if (j == grid.jEnd - 1)
			{
				vectorB(index) += -lambda* u(i, j + 1) / grid.dy2;
			}
		}
	}
}

template<class TT>
inline void SurfaceReconst<TT>::OptimalU(const Field2D<double>& f, const Field2D<Vector2D<double>>& d, const Field2D<Vector2D<double>>& b, const CSR<double>& csrA, Field2D<double>& u)
{
	cout << "Start Optimal U" << endl;


	VectorND<double> vectorB = VectorND<double>((u.iRes - 2)*(u.jRes - 2));
	GenerateLinearSystem(u, f, d, b, vectorB);
	ofstream solutionFile1;
	//solutionFile1.open("D:\\Data/vectorB.txt", ios::binary);
	//for (int i = 0; i <= vectorB.iEnd; i++)
	//{
	//		solutionFile1 << vectorB(i) << endl;
	//}
	//solutionFile1.close();

	VectorND<double> uStar = CG(csrA, vectorB);

	int idx;
	int innerIRes = u.iRes - 2;

#pragma omp parallel for private(idx)
	for (int i = u.iStart + 1; i <= u.iEnd - 1; i++)
	{
		for (int j = u.jStart + 1; j <= u.jEnd - 1; j++)
		{
			idx = (i - 1) + (j - 1)*innerIRes;
			u(i, j) = uStar(idx);
		}
	}


	solutionFile1.open("D:\\Data/innerU.txt", ios::binary);
	for (int i = grid.iStart + 1; i <= grid.iEnd - 1; i++)
	{
		for (int j = grid.jStart + 1; j <= grid.jEnd - 1; j++)
		{
			solutionFile1 << i << " " << j << " " << grid(i, j) << " " << u(i, j) << endl;
		}
	}
	solutionFile1.close();

	cout << "End Optimal U" << endl;
	cout << endl;
}

template<class TT>
inline void SurfaceReconst<TT>::OptimalD(const Field2D<Vector2D<double>>& gradientU, const Field2D<Vector2D<double>>& b, Field2D<Vector2D<double>>& d)
{
#pragma omp parallel for
	for (int i = d.iStart; i <= d.iEnd; i++)
	{
		for (int j = d.jStart; j <= d.jEnd; j++)
		{
			d(i, j) = BregmanMethod<double>::Shrink(gradientU(i, j) + b(i, j), lambda);
		}
	}
}

template<class TT>
inline void SurfaceReconst<TT>::OptimalB(const Field2D<Vector2D<double>>& gradientU, const Field2D<Vector2D<double>>& d, Field2D<Vector2D<double>>& b)
{
#pragma omp parallel for
	for (int i = b.iStart; i <= b.iEnd; i++)
	{
		for (int j = b.jStart; j <= b.jEnd; j++)
		{
			b(i, j) = b(i, j) + gradientU(i, j) - d(i, j);
		}
	}
}
