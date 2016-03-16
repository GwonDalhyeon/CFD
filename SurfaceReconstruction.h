#pragma once
#include "AdvectionMethod2D.h"
#include "PoissonSolver.h"


double AdvectionMethod2D<double>::alpha;

template <class TT>
class SurfaceReconst
{
public:
	Grid2D grid;
	LevelSet2D levelSet;
	Field2D<double> distance;
	Field2D<double> velocity;

	Array2D<Vector2D<double>> givenPoint;
	int givenPointNum;
	double dt;
	int maxIteration;

	double LpNorm;
	double distanceThreshold;
	double curvatureThreshold;

	SurfaceReconst();
	~SurfaceReconst();

	SurfaceReconst(const Grid2D& ipGrid);
	SurfaceReconst(const Grid2D& ipGrid, const Array2D<Vector2D<double>>& ipData);

	void distance2Data();
	void distance2Data(const int& i, const int& j);
	void exactDistance();

	// surface reconstruction problem solver
	void initialCondition(int example);

	void computeVelocity();
	double computeIntegralTerm();
	void levelSetPropagatingTVDRK3();

	void surfaceReconstructionSolver(int example);

	bool stoppingCriterion();


	// write results
	void outputResult();
	void outputResult(const int& iter);
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
inline SurfaceReconst<TT>::SurfaceReconst(const Grid2D & ipGrid)
{
	grid = Grid2D(ipGrid);
	levelSet = LevelSet2D(ipGrid);
	distance = Field2D<double>(ipGrid);
	velocity = Field2D<double>(ipGrid);
	dt = grid.dx*grid.dx / 2.0;
	LpNorm = 2;
}

template<class TT>
inline SurfaceReconst<TT>::SurfaceReconst(const Grid2D & ipGrid, const Array2D<Vector2D<double>>& ipData)
{
	grid = Grid2D(ipGrid);
	levelSet = LevelSet2D(ipGrid);
	distance = Field2D<double>(ipGrid);
	velocity = Field2D<double>(ipGrid);
	givenPoint = ipData;
	givenPointNum = ipData.jEnd;
	dt = grid.dx*grid.dx / 2.0;
	LpNorm = 2;
}

template<class TT>
inline void SurfaceReconst<TT>::distance2Data()
{
	TT xMin, yMin;
	double h = min(grid.dx, grid.dy);

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
inline void SurfaceReconst<TT>::distance2Data(const int & i, const int & j)
{
	double tempDist;
	distance(i, j) = 100;
	for (int k = 0; k < givenPointNum; k++)
	{
		tempDist = (grid(i, j) - givenPoint(k)).magnitude();
		if (distance(i, j) > tempDist)
		{
			distance(i, j) = tempDist;
		}
	}
}

template<class TT>
inline void SurfaceReconst<TT>::exactDistance()
{
#pragma omp parallel for
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			distance2Data(i, j);
		}
	}
}

template<class TT>
inline void SurfaceReconst<TT>::initialCondition(int example)
{
	if (example == 1)
	{
		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		distance.dataArray = 100;
		velocity = Field2D<double>(grid);
		givenPoint = Array2D<Vector2D<double>>(1, 18);
		dt = grid.dx*grid.dy / 2.0;
		maxIteration = 500;
		LpNorm = 2;
		distanceThreshold = 10 * grid.dx;
		curvatureThreshold = (grid.xMax - grid.xMin) * 10;

		for (int i = 0; i < 5; i++)
		{
			givenPoint(i) = Vector2D<double>(0.3 + 0.1*double(i), 0.3) + grid.dx / 2;
			givenPoint(13 + i) = Vector2D<double>(0.3 + 0.1*double(i), 0.7) + grid.dx / 2;
		}
		for (int i = 0; i < 4; i++)
		{
			givenPoint(5 + i) = Vector2D<double>(0.3, 0.3 + 0.1*double(i)) + grid.dx / 2;
			givenPoint(9 + i) = Vector2D<double>(0.7, 0.3 + 0.1*double(i)) + grid.dx / 2;
		}
		distance2Data();

		// level set initialzation
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (abs(grid(i, j)(0) - 0.5)>abs(grid(i, j)(1) - 0.5))
				{
					levelSet(i, j) = abs(grid(i, j)(0) - 0.5) - 0.25;
				}
				else
				{
					levelSet(i, j) = abs(grid(i, j)(1) - 0.5) - 0.25;
				}
			}
		}




	}
	else if (example == 2)
	{

		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		distance.dataArray = 100;
		velocity = Field2D<double>(grid);
		givenPointNum = 100;
		givenPoint = Array2D<Vector2D<double>>(1, givenPointNum);
		dt = 5 * grid.dx*grid.dy / 2.0;
		maxIteration = 500;
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
		exactDistance();

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
	else if (example == 3)
	{
		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		distance.dataArray = 100;
		velocity = Field2D<double>(grid);
		givenPointNum = 100;
		givenPoint = Array2D<Vector2D<double>>(1, givenPointNum);
		dt = 5 * grid.dx*grid.dy / 2.0;
		maxIteration = 500;
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
		exactDistance();

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
		grid = Grid2D(0, 1, 101, 0, 1, 101);
		levelSet = LevelSet2D(grid);
		distance = Field2D<double>(grid);
		distance.dataArray = 100;
		velocity = Field2D<double>(grid);
		givenPointNum = 100;
		givenPoint = Array2D<Vector2D<double>>(1, givenPointNum);
		dt = 5 * grid.dx*grid.dy / 2.0;
		maxIteration = 5000;
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
		exactDistance();

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


}

template<class TT>
inline void SurfaceReconst<TT>::computeVelocity()
{
	double integralTerm = computeIntegralTerm();
	double lastTerm;
	Vector2D <double> gradPhi;

	double tempDist;
	double tempCurvature;

#pragma omp parallel for private(lastTerm, tempDist,tempCurvature, gradPhi)
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
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
			velocity(i, j) = gradPhi.magnitude()*integralTerm*pow(min(distance(i, j), distanceThreshold), LpNorm - 1)*lastTerm;

			//if (i < 1 && j < 1)
			//{
			//	cout << i << " " << j << endl;
			//	cout << "inte : " << integralTerm << endl;
			//	cout << "dist : " << distance(i, j) << endl;
			//	cout << "disG : " << distance.gradient(i, j) << endl;
			//	cout << "phiG : " << gradPhi << endl;
			//	cout << "phiG : " << gradPhi.magnitude() << endl;
			//	cout << "dotP : " << dotProduct(distance.gradient(i, j), gradPhi / (gradPhi.magnitude() + DBL_EPSILON)) << endl;
			//	cout << "curv : " << levelSet.meanCurvature(i, j) << endl;
			//	cout << "last : " << lastTerm << endl;
			//	cout << "velo : " << velocity(i, j) << endl;
			//	cout << endl;
			//}
		}
	}
}

template<class TT>
inline double SurfaceReconst<TT>::computeIntegralTerm()
{
	double sum = 0.0;
	levelSet.computeMeanCurvature();
#pragma omp parallel for reduction(+:sum)
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
		
			sum += pow(distance(i, j), LpNorm)*AdvectionMethod2D<double>::deltaFt(levelSet(i, j))*levelSet.gradient(i, j).magnitude()*grid.dx*grid.dy;
			//cout << i << " " << j << endl;
			//cout << "distance : " << pow(distance(i, j), LpNorm) << endl;
			//cout << "deltaft  : " << AdvectionMethod2D<double>::deltaFt(levelSet(i, j)) << endl;
			//cout << "gradient : " << levelSet.gradient(i, j) << endl;
			//cout << "magnitude: " << levelSet.gradient(i, j).magnitude() << endl;
			//cout << sum << endl;
			//cout << endl;
			//cout << i << " " << j << " " << sum << endl;
		}
	}
	return pow(sum, 1.0 / LpNorm - 1.0);
}

template<class TT>
inline void SurfaceReconst<TT>::levelSetPropagatingTVDRK3()
{
	LevelSet2D originLevelSet = levelSet;
	LevelSet2D tempLevelSet(originLevelSet.grid);

	Field2D<TT> k1(levelSet.grid);
	Field2D<TT> k2(levelSet.grid);
	Field2D<TT> k3(levelSet.grid);

	Field2D<TT> wenoXMinus(levelSet.grid);
	Field2D<TT> wenoXPlus(levelSet.grid);
	Field2D<TT> wenoYMinus(levelSet.grid);
	Field2D<TT> wenoYPlus(levelSet.grid);

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
inline void SurfaceReconst<TT>::surfaceReconstructionSolver(int example)
{
	initialCondition(example);

	AdvectionMethod2D<double>::alpha = min(grid.dx, grid.dy); // delta function parameter.

	//for (int j = 0; j < 100; j++)
	//{
	//	cout << "Reinitialization : " << j+1 << endl;
	//	AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
	//}

	outputResult(0);
	outputResult();
	
	int i = 0;
	while (!stoppingCriterion() && i < maxIteration)
	{
		i++;
		cout << "Iteration : " << i << endl;

		computeVelocity();

		levelSetPropagatingTVDRK3();
		//AdvectionMethod2D<double> ::levelSetPropagatingTVDRK3(levelSet, dt);
		//AdvectionMethod2D<double>::levelSetPropagatingEuler(levelSet, velocity, dt);


		for (int j = 0; j < 20; j++)
		{
			cout << "Reinitialization : " << i << "-" << j + 1 << endl;
			AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
		}

		if (i % 10 == 0)
		{
			outputResult(i);
		}
		cout << endl;
	}

	outputResult(i);
}

template<class TT>
inline bool SurfaceReconst<TT>::stoppingCriterion()
{
	bool criterion = true;

	Vector2D<int> test(1, 1);
	levelSet(test);

	for (int i = 0; i < givenPointNum; i++)
	{
		if (abs(levelSet.interpolation(givenPoint(i)))>grid.dx + grid.dy)
		{
			criterion = false;
			exit;
		}
	}

	return criterion;
}

template<class TT>
inline void SurfaceReconst<TT>::outputResult()
{
	//clock_t before;
	//double  result;
	//before = clock();

	ofstream solutionFile1;
	solutionFile1.open("D:\\Data/phi.txt", ios::binary);

	ofstream solutionFile2;
	solutionFile2.open("D:\\Data/distance.txt", ios::binary);
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			solutionFile1 << i << " " << j << " " << grid(i, j) << " " << levelSet(i, j) << endl;
			solutionFile2 << i << " " << j << " " << grid(i, j) << " " << distance(i, j) << endl;
		}
	}
	solutionFile1.close();
	solutionFile2.close();


	ofstream solutionFile3;
	solutionFile3.open("D:\\Data/pointData.txt", ios::binary);
	for (int i = givenPoint.iStart; i <= givenPoint.iEnd; i++)
	{
		for (int j = givenPoint.jStart; j <= givenPoint.jEnd; j++)
		{
			solutionFile3 << givenPoint(i, j) << endl;
		}
	}
	solutionFile3.close();


	//result = (double)(clock() - before) / CLOCKS_PER_SEC;
	//cout << "binary : " << result << "\n";
	////printf("걸린시간은 %5.2f 입니다.\n", result);
}

template<class TT>
inline void SurfaceReconst<TT>::outputResult(const int & iter)
{
	ofstream solutionFile1;
	solutionFile1.open("D:\\Data/phi" + to_string(iter) + ".txt", ios::binary);
	ofstream solutionFile2;
	solutionFile2.open("D:\\Data/velocity" + to_string(iter) + ".txt", ios::binary);

	//clock_t before;
	//double  result;

	//before = clock();
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			solutionFile1 << i << " " << j << " " << grid(i, j) << " " << levelSet(i, j) << endl;
			solutionFile2 << i << " " << j << " " << grid(i, j) << " " << velocity(i, j) << endl;
		}
	}
	solutionFile1.close();
	solutionFile2.close();

	//result = (double)(clock() - before) / CLOCKS_PER_SEC;
	//printf("걸린시간은 %5.2f 입니다.\n", result);
}
