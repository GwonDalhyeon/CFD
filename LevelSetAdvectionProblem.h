#pragma once
#include "CommonDef.h"
#include "Grid2D.h"
#include "Field2D.h"
#include "LevelSet2D.h"
#include "AdvectionMethod2D.h"

class LevelSetAdvection
{
public:
	Grid2D grid;
	LevelSet2D levelSet;

	double cflCondition;
	double dt;

	bool isPropagation, isReinitialization;

	// Level set propagation velocity.
	bool isVelocity;
	bool needReinitial;
	Field2D<double> velocityX;
	Field2D<double> velocityY;

	int reinitialIter;
	int propaMaxIteration;
	int propaWriteIter;


	// Zero level set point.
	int givenPointNum;
	VectorND<Vector2D<double>> givenPoint;


	// Reinitialization Test Variable
	LevelSet2D exactLevelSet;

	int reinitialMaxIteration;
	int reinitialWriteIter;


	// Surface reconstruction Variable
	Field2D<double> distance;
	Field2D<double> reconstructionVelocity;

	int reconstMaxIteration;
	int reconstWriteIter;

	double LpNorm;
	double distanceThreshold;
	double curvatureThreshold;




	LevelSetAdvection();
	~LevelSetAdvection();




	void InitialCondition(const int & example, const bool & propa, const bool & reinitial);
	void AdvectionSolver(const int & example, const bool & propa, const bool & reinitial, const double & cfl);




	// Propagation.
	void PropaInitialCondition(const int& example);

	// Reinitialization
	void ReinitializationInitialCondition(const int& example);


	// Surface reconstruction : Split Bregman Method
	// Notation : Geodesic Application of the Split Bregman Method ... Osher.
	double lambda;
	double mu;
	void SurfaceReconstructionSplitBregman(const int & example, const bool & propa, const bool & reinitial, const bool & surfReconst);
	void GenerateLinearSystem(Array2D<double>& matrixA);
	void GenerateLinearSystem(const Field2D<double>& u, const Field2D<double>& f, const Field2D<Vector2D<double>>& d, const Field2D<Vector2D<double>>& b, VectorND<double>& vectorB);
	void OptimalU(const Field2D<double>& f, const Field2D<Vector2D<double>>& d, const Field2D<Vector2D<double>>& b, const CSR<double>& csrA, Field2D<double>& u);
	void OptimalD(const Field2D<Vector2D<double>>& gradientU, const Field2D<Vector2D<double>>& b, Field2D<Vector2D<double>>& d);
	void OptimalB(const Field2D<Vector2D<double>>& gradientU, const Field2D<Vector2D<double>>& d, Field2D<Vector2D<double>>& b);


	// Adaptive time step functions.
	double AdaptiveTimeStep();
	double AdaptiveTimeStep(const Field2D<double>& velocity1);
	double AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2);

	void OutputResult(const int& iter);

	void OutputResult(const int& iter, const Field2D<double>& ipField, const string& fileName);

	void OutputResult(const int& iter, const VectorND<Vector2D<double>>& ipVector, const string& fileName);



private:

};

//#endif // !LevelSetAdvectionProblem_H




LevelSetAdvection::LevelSetAdvection()
{
}

LevelSetAdvection::~LevelSetAdvection()
{
}

inline void LevelSetAdvection::InitialCondition(const int& example, const bool & propa, const bool & reinitial)
{
	isPropagation = propa;
	isReinitialization = reinitial;

	if (isPropagation)
	{
		PropaInitialCondition(example);
	}
	else
	{
		ReinitializationInitialCondition(example);
	}
}


inline void LevelSetAdvection::AdvectionSolver(const int & example, const bool & propa, const bool & reinitial, const double & cfl)
{
	cflCondition = cfl;
	InitialCondition(example, propa, reinitial);

	grid.Variable();
	levelSet.phi.Variable("phi0");

	string str;
	const char*cmd;

	// Reinitialization
	if (isReinitialization)
	{
		OutputResult(0);
		for (int i = 1; i <= reinitialMaxIteration; i++)
		{

			cout << "Reinitialization : " << i << endl;
			dt = AdaptiveTimeStep();
			AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
			if (i%reinitialWriteIter == 0)
			{
				OutputResult(i);
			}
		}
	}


	// Propagation.
	if (isPropagation)
	{
		OutputResult(0);

		for (int i = 1; i <= propaMaxIteration; i++)
		{
			cout << "Level set advection : " << i << endl;

			if (isVelocity)
			{
				dt = AdaptiveTimeStep(velocityX, velocityY);
				AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, velocityX, velocityY, dt);
			}
			else
			{
				dt = AdaptiveTimeStep();
				AdvectionMethod2D<double>::levelSetPropagatingTVDRK3(levelSet, dt);
			}


			if (needReinitial)
			{
				for (int j = 0; j < reinitialIter; j++)
				{
					cout << "Reinitialization : " << i << "-" << j + 1 << endl;
					dt = AdaptiveTimeStep();
					AdvectionMethod2D<double>::levelSetReinitializationTVDRK3(levelSet, dt);
				}
			}

			if (i%propaWriteIter == 0)
			{
				OutputResult(i);
			}
		}
	}
}


inline void LevelSetAdvection::PropaInitialCondition(const int & example)
{
	if (example == 1)
	{
		isVelocity = true;
		needReinitial = false;

		cout << "Level set advection Test - Rotating a circle" << endl;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);

		levelSet = LevelSet2D(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt((grid(i, j).x - 0.5)*(grid(i, j).x - 0.5) + grid(i, j).y*grid(i, j).y) - 0.25;
			}
		}

		velocityX = Field2D<double>(grid);
		velocityY = Field2D<double>(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				velocityX(i, j) = -grid(i, j).y;
				velocityY(i, j) = grid(i, j).x;
			}
		}

		isVelocity = true;

		//dt = grid.dx*grid.dy;
		propaMaxIteration = 10000;
		propaWriteIter = 100;
	}
	else if (example == 2)
	{
		cout << "Level set advection Test - Spiral" << endl;

		isVelocity = false;
		needReinitial = false;

		grid = Grid2D(-1, 1, 201, -1, 1, 201);


		levelSet = LevelSet2D(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = sqrt((grid(i, j).x - 0.5)*(grid(i, j).x - 0.5) + grid(i, j).y*grid(i, j).y) - 0.25;
			}
		}


		//dt = grid.dx*grid.dy;
		propaMaxIteration = 1000;
		propaWriteIter = 10;
	}
	else if (example == 3)
	{
		cout << "Level set advection Test - Seven-point Star" << endl;

		isVelocity = true;
		needReinitial = false;

		grid = Grid2D(-0.25, 0.25, 101, -0.25, 0.25, 101);

		givenPointNum = 100;
		givenPoint = VectorND<Vector2D<double>>(givenPointNum);

		double s;
#pragma omp parallel for private (s)
		for (int i = 0; i < givenPointNum; i++)
		{
			s = double(i) / double(givenPointNum);
			givenPoint(i) = (0.1 + 0.065*sin(7 * 2 * PI*s))*Vector2D<double>(cos(2 * PI*s), sin(2 * PI*s));
		}

		levelSet = LevelSet2D(grid);
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.iEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{

				//levelSet(i, j) = distance2Data(i, j);
			}
		}

		isVelocity = false;

		//dt = grid.dx*grid.dy;
		propaMaxIteration = 1000;
		propaWriteIter = 10;
	}
}

inline void LevelSetAdvection::ReinitializationInitialCondition(const int & example)
{
	grid = Grid2D(-2, 2, 101, -2, 2, 101);

	//dt = grid.dx*grid.dy;
	reinitialMaxIteration = 1000;
	reinitialWriteIter = 10;

	exactLevelSet = LevelSet2D(grid);
	levelSet = LevelSet2D(grid);

	isVelocity = false;
	needReinitial = false;

	double a = 0.7;
	double r = 1.0;

	if (example == 1) //// Exmple1. A circle with center at the origen and radius 1.
	{
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.jEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				levelSet(i, j) = (grid(i, j).magnitude() - 1.0)*((grid(i, j) - 1.0).magnitude() + 0.1);
				exactLevelSet(i, j) = grid(i, j).magnitude() - 1.0;
			}
		}
	}
	else if (example == 2) //// Exmple2. Two circles of radius r are placed at (+-a,0)  and a sqruar on the plane. Let 0<a<r, sh that the two circles intersect each other.
	{

		double temp1, temp2, temp3;
#pragma omp parallel for private(temp1,temp2,temp3)
		for (int i = grid.iStart; i <= grid.jEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{

				temp1 = (a - grid(i, j).x) / sqrt((a - grid(i, j).x)*(a - grid(i, j).x) + grid(i, j).y * grid(i, j).y);
				temp2 = (a + grid(i, j).x) / sqrt((a + grid(i, j).x)*(a + grid(i, j).x) + grid(i, j).y * grid(i, j).y);
				if (temp1 >= a / r && temp2 >= a / r)
				{
					temp3 = min(grid(i, j).x * grid(i, j).x + (grid(i, j).y + sqrt(r*r - a*a))*(grid(i, j).y + sqrt(r*r - a*a)), grid(i, j).x * grid(i, j).x + (grid(i, j).y - sqrt(r*r - a*a))*(grid(i, j).y - sqrt(r*r - a*a)));
					temp3 = sqrt(temp3);
					levelSet(i, j) = temp3 * ((grid(i, j) - 1).magnitude2() + 0.1);
					exactLevelSet(i, j) = temp3;
				}
				else
				{
					temp3 = min(sqrt((grid(i, j).x + a)*(grid(i, j).x + a) + grid(i, j).y * grid(i, j).y), sqrt((grid(i, j).x - a)*(grid(i, j).x - a) + grid(i, j).y * grid(i, j).y)) - r;
					levelSet(i, j) = temp3 * ((grid(i, j) - 1).magnitude2() + 0.1);
					exactLevelSet(i, j) = temp3;
				}
			}
		}

	}
	else if (example == 3) ////Exmple3.Two circles of radius r are placed at(+-a, 0) on the plane.Let 0<a<r, sh that the two circles intersect each other.
	{
		//double temp1, temp2, temp3;
		//for (int i = grid.iStart; i <= grid.jEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		temp1 = (a - grid(i, j).x) / sqrt((a - grid(i, j).x)*(a - grid(i, j).x) + grid(i, j).y * grid(i, j).y);
		//		temp2 = (a + grid(i, j).x) / sqrt((a + grid(i, j).x)*(a + grid(i, j).x) + grid(i, j).y * grid(i, j).y);
		//		if (temp1 >= a / r && temp2 >= a / r)
		//		{
		//			temp3 = min(grid(i, j).x * grid(i, j).x + (grid(i, j).y + sqrt(r*r - a*a))*(grid(i, j).y + sqrt(r*r - a*a)), grid(i, j).x * grid(i, j).x + (grid(i, j).y - sqrt(r*r - a*a))*(grid(i, j).y - sqrt(r*r - a*a)));
		//			levelSet(i, j) = sqrt(temp3) * ((grid(i, j) - 1).magnitude2() + 0.1);
		//			exactLevelSet(i, j) = sqrt(temp3);
		//		}
		//		else
		//		{
		//			temp3 = min(sqrt((grid(i, j).x + a)*(grid(i, j).x + a) + grid(i, j).y * grid(i, j).y) - r, sqrt((grid(i, j).x - a)*(grid(i, j).x - a) + grid(i, j).y * grid(i, j).y) - r);
		//			levelSet(i, j) = temp3 * ((grid(i, j) - 1).magnitude2() + 0.1);
		//			exactLevelSet(i, j) = temp3;
		//		}
		//	}
		//}
	}
	else if (example == 4)
	{
#pragma omp parallel for
		for (int i = grid.iStart; i <= grid.jEnd; i++)
		{
			for (int j = grid.jStart; j <= grid.jEnd; j++)
			{
				if (grid(i, j).magnitude() < 1)
				{
					levelSet(i, j) = -1.0;

				}
				else
				{
					levelSet(i, j) = 0.5;

				}
				exactLevelSet(i, j) = grid(i, j).magnitude() - 1.0;
			}
		}
	}
	else if (example == 5)
	{
		//for (int i = grid.iStart; i <= grid.jEnd; i++)
		//{
		//	for (int j = grid.jStart; j <= grid.jEnd; j++)
		//	{
		//		if (abs(grid(i,j).x) < 0.5 || abs(grid(i, j).y)<0.5)
		//		{
		//			levelSet(i,j) = -1.0;
		//			exactLevelSet(i,j) = ;
		//		}
		//		else
		//		{
		//			levelSet(i, j) = 1.0;
		//		}
		//	}
		//}
	}
}



inline void LevelSetAdvection::SurfaceReconstructionSplitBregman(const int & example, const bool & propa, const bool & reinitial, const bool & surfReconst)
{

	InitialCondition(example, propa, reinitial);

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
			distance(i, j) = -distance(i,j);
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

inline void LevelSetAdvection::GenerateLinearSystem(Array2D<double>& matrixA)
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

inline void LevelSetAdvection::GenerateLinearSystem(const Field2D<double>& u, const Field2D<double>& f, const Field2D<Vector2D<double>>& d, const Field2D<Vector2D<double>>& b, VectorND<double>& vectorB)
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

inline void LevelSetAdvection::OptimalU(const Field2D<double>& f, const Field2D<Vector2D<double>>& d, const Field2D<Vector2D<double>>& b, const CSR<double>& csrA, Field2D<double>& u)
{
	cout << "Start Optimal U" << endl;


	VectorND<double> vectorB = VectorND<double>((u.iRes-2)*(u.jRes - 2));
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
	for (int i = u.iStart+1; i <= u.iEnd-1; i++)
	{
		for (int j = u.jStart+1; j <= u.jEnd-1; j++)
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

inline void LevelSetAdvection::OptimalD(const Field2D<Vector2D<double>>& gradientU, const Field2D<Vector2D<double>>& b, Field2D<Vector2D<double>>& d)
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

inline void LevelSetAdvection::OptimalB(const Field2D<Vector2D<double>>& gradientU, const Field2D<Vector2D<double>>& d, Field2D<Vector2D<double>>& b)
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


inline double LevelSetAdvection::AdaptiveTimeStep()
{
	return cflCondition*max(grid.dx, grid.dy);
}

inline double LevelSetAdvection::AdaptiveTimeStep(const Field2D<double>& velocity1)
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

inline double LevelSetAdvection::AdaptiveTimeStep(const Field2D<double>& velocity1, const Field2D<double>& velocity2)
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

inline void LevelSetAdvection::OutputResult(const int & iter)
{
	cout << "Write results" << endl;

	string fileName = "phi" + to_string(iter);
	levelSet.phi.WriteFile(fileName);

	if (false)
	{
		fileName = "velocity" + to_string(iter);
		reconstructionVelocity.WriteFile(fileName);
	}


	if (iter == 0)
	{
		if (isPropagation)
		{
			fileName = "velocityX";
			velocityX.WriteFile(fileName);
			fileName = "velocityY";
			velocityY.WriteFile(fileName);
		}

		if (false)//reconstruction
		{
			fileName = "distance";
			distance.WriteFile(fileName);

			fileName = "pointData";
			givenPoint.WriteFile(fileName);
		}


	}
}

inline void LevelSetAdvection::OutputResult(const int & iter, const Field2D<double>& ipField, const string & fileName)
{
	cout << "Write results" << endl;

	ofstream solutionFile1;
	solutionFile1.open("D:\\Data/"+ fileName + to_string(iter) + ".txt", ios::binary);
	for (int i = grid.iStart; i <= grid.iEnd; i++)
	{
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			solutionFile1 << i << " " << j << " " << grid(i, j) << " " << ipField(i, j) << endl;
		}
	}
	solutionFile1.close();
}

inline void LevelSetAdvection::OutputResult(const int & iter, const VectorND<Vector2D<double>>& ipVector, const string & fileName)
{
	cout << "Write results" << endl;

	ofstream solutionFile4;
	solutionFile4.open("D:\\Data/" + fileName + ".txt", ios::binary);
	for (int i = 0; i <= givenPoint.iEnd; i++)
	{
		solutionFile4 << givenPoint(i) << endl;
	}
	solutionFile4.close();
}
