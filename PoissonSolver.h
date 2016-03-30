#pragma once



//#ifndef PoissonSolver_H
//#define PoissonSolver_H
#include "CommonDef.h"
#include "VectorND.h"
#include "Grid2D.h"
#include "Array2D.h"
#include "CSR.h"
#include "LevelSet2D.h"
#include "LinearSolver.h"


class PoissonSolver
{
public:
	Grid2D grid;
	Grid2D innerGrid;
	Array2D<double> poissonMatrix;
	VectorND<double> poissonVector;

	Field2D<double> solution;
	VectorND<double> innerSolution;

	CSR<double> poissonCSR;
	LevelSet2D levelSet;

	PoissonSolver();
	~PoissonSolver();

	int index(int i, int j);
	int indexInner(int i, int j);
	//int indexVec(int i, int j);
	int indexMat(int i, int j);

	// A Boundary Condition Capturing Method
	PoissonSolver(const Grid2D& ipGrid, const LevelSet2D& ipLevelSet, const Field2D<double>& ipBeta, const Field2D<double>& ipF, const Field2D<double>& ipjCondition1, const Field2D<double>& ipjCondition2);

	void GenerateJumpCondi(int example, Field2D<double>& beta, Field2D<double>& f, Field2D<double>& jCondition1, Field2D<double>& jCondition2);
	void GeneratePoissonMatrixJumpCondi(const Field2D<double>& ipBeta, const Field2D<double>& ipF, const Field2D<double>& ipjCondition1, const Field2D<double>& ipjCondition2);
	void GeneratePoissonVectorJumpCondi(const Field2D<double>& ipBeta, const Field2D<double>& ipF, const Field2D<double>& ipjCondition1, const Field2D<double>& ipjCondition2);
	void SolvePoissonJumpCondi(int example, const Grid2D& ipGrid);
	void OutputResult();


private:

};

//#endif // !PoissonSolver_H



PoissonSolver::PoissonSolver()
{
}


PoissonSolver::~PoissonSolver()
{
}

PoissonSolver::PoissonSolver(const Grid2D& ipGrid, const LevelSet2D& ipLevelSet, const Field2D<double>& ipBeta, const Field2D<double>& ipF, const Field2D<double>& ipjCondition1, const Field2D<double>& ipjCondition2)
{
	grid = ipGrid;
	solution = Field2D<double>(grid);
	levelSet = ipLevelSet;

	Grid2D innerGrid = Grid2D(ipGrid.xMin + ipGrid.dx, ipGrid.xMax - ipGrid.dx, 1, ipGrid.iRes - 2, ipGrid.yMin + ipGrid.dy, ipGrid.yMax - ipGrid.dy, 1, ipGrid.jRes - 2);
	VectorND<double> innerSolution(innerGrid.iRes*innerGrid.jRes);
	Field2D<double> beta = ipBeta;;
	Field2D<double> f = ipF;
	Field2D<double> jCondition1 = ipjCondition1;
	Field2D<double> jCondition2 = ipjCondition2;
	//double leftBdry, rightBdry;

	poissonMatrix = Array2D<double>(1, innerGrid.iRes, 1, innerGrid.jRes);
	poissonVector = VectorND<double>(poissonMatrix.ijRes);
}



inline int PoissonSolver::index(int i, int j)
{
	return i - grid.iStart + (j - innerGrid.jStart)*grid.iRes;
}

inline int PoissonSolver::indexInner(int i, int j)
{
	return i - innerGrid.iStart + (j - innerGrid.jStart)*innerGrid.iRes;
}

inline int PoissonSolver::indexMat(int i, int j)
{
	return 0;
	//return (i - innerGrid.iStart)*innerGrid.iRes*innerGrid.jRes + (i - innerGrid.iStart) + (j - innerGrid.jStart)*innerGrid.iRes*(innerGrid.iRes*innerGrid.jRes + 1);

	//int matIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + (i - 1) + (j - 1)*innerGrid.iRes*(innerGrid.iRes*innerGrid.jRes + 1);
	//int matLeftIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + (i - 1 - 1) + (j - 1)*innerGrid.iRes*(innerGrid.iRes*innerGrid.jRes + 1);
	//int matRightIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + i + (j - 1)*innerGrid.iRes*(innerGrid.iRes*innerGrid.jRes + 1);
	//int matBottomIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + (i - 1) + (j - 1)*innerGrid.iRes*innerGrid.iRes*innerGrid.jRes + (j - 1 - 1)*innerGrid.iRes;
	//int matTopIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + (i - 1) + (j - 1)*innerGrid.iRes*innerGrid.iRes*innerGrid.jRes + (j)*innerGrid.iRes;
}


inline void PoissonSolver::GenerateJumpCondi(int example, Field2D<double>& beta, Field2D<double>& f, Field2D<double>& jCondition1, Field2D<double>& jCondition2)
{
	if (example == 1)
	{
		///////////////////////////////////////////////////////////////////////////////
		////
		////     Poisson equation 1D. Example 1
		////
		///////////////////////////////////////////////////////////////////////////////

		//double X0, X1, Y0, Y1, Z0, Z1;
		//double deltaX, deltaY, deltaZ;
		//int numX, numY, numZ;
		//X0 = 0; X1 = 1; numX = 101;
		//Y0 = 0; Y1 = 1; numY = 101;
		//Z0 = 0; Z1 = 1; numZ = 101;

		//levelSet = LevelSet(grid);
		//for (int i = 0; i < grid.numX; i++)
		//{
		//	levelSet.phi[i] = abs(grid.x[i] - 0.45) - 0.15 - grid.deltaX / 2;
		//}

		//leftBdry = 0, rightBdry = 0;

		//for (int i = 0; i < grid.numX; i++)
		//{
		//	if (levelSet.phi[i] <= 0)
		//	{
		//		beta[i] = 2;
		//	}
		//	else
		//	{
		//		beta[i] = 1;
		//	}
		//}

		//for (int i = 0; i < grid.numMatX; i++)
		//{
		//	if (levelSet.phi[i + 1] <= 0)
		//	{
		//		f[i] = (8 * grid.x[i + 1] * grid.x[i + 1] - 4)*exp(-grid.x[i + 1] * grid.x[i + 1]);
		//		poissonVector[i] = f[i];
		//		//cout<<i<<" "<<f[i]<<endl;
		//	}
		//	//else
		//	//{
		//	//	f[i] = 0;
		//	//}
		//}

		////for (int i = 0; i < grid.numMatX; i++)
		////{
		////	poissonVector[i] = f[i];
		////}

		////for (int i = 0; i < grid.numX; i++)
		////{
		////	jCondition1[i] = 0;
		////	jCondition2[i] = 0;
		////}

		//jCondition1[29] = -exp(-0.09);
		//jCondition1[30] = -exp(-0.09);
		//jCondition2[29] = -1.2*exp(-0.09);
		//jCondition2[30] = -1.2*exp(-0.09);

		//jCondition1[60] = -exp(-0.36);
		//jCondition1[61] = -exp(-0.36);
		//jCondition2[60] = 2.4*exp(-0.36);
		//jCondition2[61] = 2.4*exp(-0.36);
	}

	if (example == 2)
	{
		///////////////////////////////////////////////////////////////////////////////
		////
		////     Poisson equation 2D. Example 2
		////
		///////////////////////////////////////////////////////////////////////////////

		//double X0, X1, Y0, Y1, Z0, Z1;
		//double deltaX, deltaY, deltaZ;
		//int numX, numY, numZ;
		//X0 = 0; X1 = 1; numX = 11;
		//Y0 = 0; Y1 = 1; numY = 11;
		//Z0 = 0; Z1 = 1; numZ = 101;

		levelSet = LevelSet2D(grid);
		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			for (int i = grid.iStart; i <= grid.iEnd; i++)
			{
				levelSet.phi(i, j) = sqrt((grid(i, j)(0) - 0.5)*(grid(i, j)(0) - 0.5) + (grid(i, j)(1) - 0.5)*(grid(i, j)(1) - 0.5)) - 0.25 - grid.dx / 2;
			}
		}

		for (int j = grid.jStart; j <= grid.jEnd; j++)
		{
			for (int i = grid.iStart; i <= grid.iEnd; i++)
			{
				if (levelSet(i, j) <= 0)
				{
					beta(i, j) = 2;
				}
				else
				{
					beta(i, j) = 1;
				}
			}
		}

		for (int j = innerGrid.jStart; j <= innerGrid.jEnd; j++)
		{
			for (int i = innerGrid.iStart; i <= innerGrid.iEnd; i++)
			{
				if (levelSet(i, j) <= 0)
				{
					f(i, j) = 8 * (grid(i, j)(0) * grid(i, j)(0) + grid(i, j)(1) * grid(i, j)(1) - 1)*exp(-grid(i, j)(0) * grid(i, j)(0) - grid(i, j)(1) * grid(i, j)(1));
				}
			}
		}

		for (int j = grid.jStart; j <= grid.jEnd - 1; j++)
		{
			for (int i = grid.iStart; i <= grid.iEnd - 1; i++)
			{
				if ((levelSet(i, j) <= 0 && levelSet(i + 1, j) > 0) || (levelSet(i, j) > 0 && levelSet(i + 1, j) <= 0))
				{
					jCondition1(i, j) = -exp(-grid(i, j)(0) * grid(i, j)(0) - grid(i, j)(1) * grid(i, j)(1));
					jCondition1(i + 1, j) = -exp(-grid(i + 1, j)(0) * grid(i + 1, j)(0) - grid(i, j)(1) * grid(i, j)(1));
					jCondition2(i, j) = 8 * (2 * grid(i, j)(0) * grid(i, j)(0) + 2 * grid(i, j)(1) * grid(i, j)(1) - grid(i, j)(0) - grid(i, j)(1))*exp(-grid(i, j)(0) * grid(i, j)(0) - grid(i, j)(1) * grid(i, j)(1));
					jCondition2(i + 1, j) = 8 * (2 * grid(i + 1, j)(0) * grid(i + 1, j)(0) + 2 * grid(i, j)(1) * grid(i, j)(1) - grid(i + 1, j)(0) - grid(i, j)(1))*exp(-grid(i + 1, j)(0) * grid(i + 1, j)(0) - grid(i, j)(1) * grid(i, j)(1));
				}
				if ((levelSet(i, j) <= 0 && levelSet(i, j + 1) > 0) || (levelSet(i, j) > 0 && levelSet(i, j + 1) <= 0))
				{
					jCondition1(i, j) = -exp(-grid(i, j)(0) * grid(i, j)(0) - grid(i, j)(1) * grid(i, j)(1));;
					jCondition1(i, j + 1) = -exp(-grid(i, j)(0) * grid(i, j)(0) - grid(i, j + 1)(1) * grid(i, j + 1)(1));;
					jCondition2(i, j) = 8 * (2 * grid(i, j)(0) * grid(i, j)(0) + 2 * grid(i, j)(1) * grid(i, j)(1) - grid(i, j)(0) - grid(i, j)(1))*exp(-grid(i, j)(0) * grid(i, j)(0) - grid(i, j)(1) * grid(i, j)(1));
					jCondition2(i, j + 1) = 8 * (2 * grid(i, j)(0) * grid(i, j)(0) + 2 * grid(i, j + 1)(1) * grid(i, j + 1)(1) - grid(i, j)(0) - grid(i, j + 1)(1))*exp(-grid(i, j)(0) * grid(i, j)(0) - grid(i, j + 1)(1) * grid(i, j + 1)(1));
				}
			}
		}

	}
}

inline void PoissonSolver::GeneratePoissonMatrixJumpCondi(const Field2D<double>& beta, const Field2D<double>& f, const Field2D<double>& jCondition1, const Field2D<double>& jCondition2)
{
	double tempBeta = 0;

	int matIndex, matLeftIndex, matRightIndex, matTopIndex, matBottomIndex;

	for (int j = innerGrid.jStart; j <= innerGrid.jEnd; j++)
	{
		//poissonMatrix[i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*grid.numMatX*grid.numMatY + j*grid.numMatX ] = 0;
		//poissonMatrix[i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*(grid.numMatX*grid.numMatY + 1) ] = 0;
		for (int i = innerGrid.iStart; i <= innerGrid.iEnd; i++)
		{
			matIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + (i - 1) + (j - 1)*innerGrid.iRes*(innerGrid.iRes*innerGrid.jRes + 1);
			matLeftIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + (i - 1 - 1) + (j - 1)*innerGrid.iRes*(innerGrid.iRes*innerGrid.jRes + 1);
			matRightIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + i + (j - 1)*innerGrid.iRes*(innerGrid.iRes*innerGrid.jRes + 1);
			matBottomIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + (i - 1) + (j - 1)*innerGrid.iRes*innerGrid.iRes*innerGrid.jRes + (j - 1 - 1)*innerGrid.iRes;
			matTopIndex = (i - 1)*innerGrid.iRes*innerGrid.jRes + (i - 1) + (j - 1)*innerGrid.iRes*innerGrid.iRes*innerGrid.jRes + (j)*innerGrid.iRes;


			if (j > innerGrid.jStart && j < innerGrid.jEnd)
			{
				if ((levelSet(i, j) > 0 && levelSet(i, j + 1) <= 0) || (levelSet(i, j) <= 0 && levelSet(i, j + 1) > 0))
				{
					tempBeta = beta(i, j) * beta(i, j + 1) * (abs(levelSet(i, j)) + abs(levelSet(i, j + 1))) / (beta(i, j + 1) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i, j + 1)));
					poissonMatrix[matBottomIndex] = -1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j - 1)) / 2;
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j - 1)) / 2 + 1 / (grid.dy*grid.dy)*(tempBeta);
					poissonMatrix[matTopIndex] = -1 / (grid.dy*grid.dy)*tempBeta;
					//cout<<i<<" "<< tempBeta<<endl;
					//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
					//cout<<"";

				}
				else if ((levelSet(i, j - 1) > 0 && levelSet(i, j) <= 0) || (levelSet(i, j - 1) <= 0 && levelSet(i, j) > 0))
				{
					tempBeta = beta(i, j - 1) * beta(i, j) * (abs(levelSet(i, j - 1)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i, j - 1)) + beta(i, j - 1) * abs(levelSet(i, j)));
					poissonMatrix[matBottomIndex] = -1 / (grid.dy*grid.dy)*tempBeta;
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j + 1)) / 2 + 1 / (grid.dy*grid.dy)*(tempBeta);
					poissonMatrix[matTopIndex] = -1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j + 1)) / 2;
					//cout<<i<<" "<< tempBeta<<endl;
					//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
					//cout<<"";
				}
				else
				{
					poissonMatrix[matBottomIndex] = -1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j - 1)) / 2;
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j + 1)) / 2 + 1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j - 1)) / 2;
					poissonMatrix[matTopIndex] = -1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j + 1)) / 2;

					//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
					//cout<<"";
				}

				if (i > innerGrid.iStart && i < innerGrid.iEnd)
				{
					if ((levelSet(i, j) > 0 && levelSet(i + 1, j) <= 0) || (levelSet(i, j) <= 0 && levelSet(i + 1, j) > 0))
					{
						tempBeta = beta(i, j) * beta(i + 1, j) * (abs(levelSet(i, j)) + abs(levelSet(i + 1, j))) / (beta(i + 1, j) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i + 1, j)));
						poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(tempBeta);
						poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*tempBeta;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else if ((levelSet(i - 1, j) > 0 && levelSet(i, j) <= 0) || (levelSet(i - 1, j) <= 0 && levelSet(i, j) > 0))
					{
						tempBeta = beta(i - 1, j) * beta(i, j) * (abs(levelSet(i - 1, j)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i - 1, j)) + beta(i - 1, j) * abs(levelSet(i, j)));
						poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*tempBeta;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(tempBeta);
						poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
						poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
				}
				else if (i == innerGrid.iStart)
				{
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
					poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2;
					//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
					//cout<<"";
				}
				else
				{
					poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
					//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0<<endl;
					//cout<<"";
				}
			}
			else if (j == innerGrid.jStart)
			{
				poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j + 1)) / 2 + 1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j - 1)) / 2;
				poissonMatrix[matTopIndex] = -1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j + 1)) / 2;
				//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
				//cout<<"";
				if (i > innerGrid.iStart && i < innerGrid.iEnd)
				{
					if ((levelSet(i, j) > 0 && levelSet(i + 1, j) <= 0) || (levelSet(i, j) <= 0 && levelSet(i + 1, j) > 0))
					{
						tempBeta = beta(i, j) * beta(i + 1, j) * (abs(levelSet(i, j)) + abs(levelSet(i + 1, j))) / (beta(i + 1, j) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i + 1, j)));
						poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(tempBeta);
						poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*tempBeta;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else if ((levelSet(i - 1, j) > 0 && levelSet(i, j) <= 0) || (levelSet(i - 1, j) <= 0 && levelSet(i, j) > 0))
					{
						tempBeta = beta(i - 1, j) * beta(i, j) * (abs(levelSet(i - 1, j)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i - 1, j)) + beta(i - 1, j) * abs(levelSet(i, j)));
						poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*tempBeta;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(tempBeta);
						poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
						poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
				}
				else if (i == innerGrid.iStart)
				{
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
					poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2;
					//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
					//cout<<"";
				}
				else
				{
					poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
					//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0<<endl;
					//cout<<"";
				}
			}
			else
			{
				poissonMatrix[matBottomIndex] = -1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j - 1)) / 2;
				poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j + 1)) / 2 + 1 / (grid.dy*grid.dy)*(beta(i, j) + beta(i, j - 1)) / 2;
				//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0 <<endl;
				//cout<<"";
				if (i > innerGrid.iStart && i < innerGrid.iEnd)
				{
					if ((levelSet(i, j) > 0 && levelSet(i + 1, j) <= 0) || (levelSet(i, j) <= 0 && levelSet(i + 1, j) > 0))
					{
						tempBeta = beta(i, j) * beta(i + 1, j) * (abs(levelSet(i, j)) + abs(levelSet(i + 1, j))) / (beta(i + 1, j) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i + 1, j)));
						poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(tempBeta);
						poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*tempBeta;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else if ((levelSet(i - 1, j) > 0 && levelSet(i, j) <= 0) || (levelSet(i - 1, j) <= 0 && levelSet(i, j) > 0))
					{
						tempBeta = beta(i - 1, j) * beta(i, j) * (abs(levelSet(i - 1, j)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i - 1, j)) + beta(i - 1, j) * abs(levelSet(i, j)));
						poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*tempBeta;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(tempBeta);
						poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
						poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
				}
				else if (i == innerGrid.iStart)
				{
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
					poissonMatrix[matRightIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2;
					//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
					//cout<<"";
				}
				else
				{
					poissonMatrix[matLeftIndex] = -1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i + 1, j)) / 2 + 1 / (grid.dx*grid.dx)*(beta(i, j) + beta(i - 1, j)) / 2;
					//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0 <<endl;
					//cout<<"";
				}
			}
		}
	}

}

inline void PoissonSolver::GeneratePoissonVectorJumpCondi(const Field2D<double>& beta, const Field2D<double>& f, const Field2D<double>& jCondition1, const Field2D<double>& jCondition2)
{
	double aGamma = 0;
	double bGamma = 0;
	double theta = 0;
	double tempBeta = 0;
	//VectorND<double> fL(innerGrid.iRes*innerGrid.jRes);
	//VectorND<double> fR(innerGrid.iRes*innerGrid.jRes);
	//VectorND<double> fB(innerGrid.iRes*innerGrid.jRes);
	//VectorND<double> fT(innerGrid.iRes*innerGrid.jRes);
	double fL, fR, fB, fT;

	Vector2D<double> normalLeft;
	Vector2D<double> normalRight;
	Vector2D<double> normalCenter;
	Vector2D<double> normalBottom;
	Vector2D<double> normalTop;

	levelSet.computeUnitNormal();

	for (int j = innerGrid.jStart; j <= innerGrid.jEnd; j++)
	{
		for (int i = innerGrid.iStart; i <= innerGrid.iEnd; i++)
		{
			normalLeft = levelSet.unitNormal(i - 1, j);
			normalRight = levelSet.unitNormal(i + 1, j);
			normalCenter = levelSet.unitNormal(i, j);
			normalBottom = levelSet.unitNormal(i, j - 1);
			normalTop = levelSet.unitNormal(i, j + 1);

			if (levelSet(i - 1, j) > 0 && levelSet(i, j) <= 0)
			{
				theta = abs(levelSet(i - 1, j)) / (abs(levelSet(i - 1, j)) + abs(levelSet(i, j)));
				aGamma = (jCondition1(i - 1, j) * abs(levelSet(i, j)) + jCondition1(i, j) * abs(levelSet(i - 1, j))) / (abs(levelSet(i - 1, j)) + abs(levelSet(i, j)));
				bGamma = (jCondition2(i - 1, j) * normalLeft[0] * abs(levelSet(i, j)) + jCondition2(i, j) * normalCenter[0] * abs(levelSet(i - 1, j))) / (abs(levelSet(i - 1, j)) + abs(levelSet(i, j)));
				tempBeta = beta(i - 1, j) * beta(i, j) * (abs(levelSet(i - 1, j)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i - 1, j)) + beta(i - 1, j) * abs(levelSet(i, j)));
				fL = tempBeta*aGamma / (grid.dx*grid.dx) - tempBeta*bGamma*theta / (beta(i - 1, j) * grid.dx);
				//fL.values[matIndex] = tempBeta*aGamma / (grid.dx*grid.dx) - tempBeta*bGamma*theta / (beta(i - 1, j) * grid.dx);
				//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
				//cout<<i<<" " <<i+1 <<endl;
				//cout << i<< " "<< fL[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else if (levelSet(i - 1, j) <= 0 && levelSet(i, j) > 0)
			{
				theta = abs(levelSet(i - 1, j)) / (abs(levelSet(i - 1, j)) + abs(levelSet(i, j)));
				aGamma = (jCondition1(i - 1, j) * abs(levelSet(i, j)) + jCondition1(i, j) * abs(levelSet(i - 1, j))) / (abs(levelSet(i - 1, j)) + abs(levelSet(i, j)));
				bGamma = (jCondition2(i - 1, j) * normalLeft[0] * abs(levelSet(i, j)) + jCondition2(i, j) * normalCenter[0] * abs(levelSet(i - 1, j))) / (abs(levelSet(i - 1, j)) + abs(levelSet(i, j)));
				tempBeta = beta(i - 1, j) * beta(i, j) * (abs(levelSet(i - 1, j)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i - 1, j)) + beta(i - 1, j) * abs(levelSet(i, j)));
				fL = -tempBeta*aGamma / (grid.dx*grid.dx) + tempBeta*bGamma*theta / (beta(i - 1, j) * grid.dx);
				//fL.values[matIndex] = -tempBeta*aGamma / (grid.dx*grid.dx) + tempBeta*bGamma*theta / (beta(i - 1, j) * grid.dx);
				//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
				//cout<<i<<" " <<i+1 <<endl;
				//cout << i<< " "<< fL[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else
			{
				fL = 0;
			}

			if (levelSet(i, j) <= 0 && levelSet(i + 1, j) > 0)
			{
				theta = abs(levelSet(i + 1, j)) / (abs(levelSet(i, j)) + abs(levelSet(i + 1, j)));
				aGamma = (jCondition1(i, j) * abs(levelSet(i + 1, j)) + jCondition1(i + 1, j) * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i + 1, j)));
				bGamma = (jCondition2(i, j) * normalCenter[0] * abs(levelSet(i + 1, j)) + jCondition2(i + 1, j) * normalRight[0] * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i + 1, j)));
				tempBeta = beta(i, j) * beta(i + 1, j) * (abs(levelSet(i, j)) + abs(levelSet(i + 1, j))) / (beta(i + 1, j) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i + 1, j)));
				fR = tempBeta*aGamma / (grid.dx*grid.dx) + tempBeta*bGamma*theta / (beta(i + 1, j) * grid.dx);
				//fR.values[matIndex] = tempBeta*aGamma / (grid.dx*grid.dx) + tempBeta*bGamma*theta / (beta(i + 1, j) * grid.dx);
				//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
				//cout<<i+1<<" " <<i+1+1 <<endl;
				//cout << i<< " "<< fR[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;

			}
			else if (levelSet(i, j) > 0 && levelSet(i + 1, j) <= 0)
			{
				theta = abs(levelSet(i + 1, j)) / (abs(levelSet(i, j)) + abs(levelSet(i + 1, j)));
				aGamma = (jCondition1(i, j) * abs(levelSet(i + 1, j)) + jCondition1(i + 1, j) * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i + 1, j)));
				bGamma = (jCondition2(i, j) * normalCenter[0] * abs(levelSet(i + 1, j)) + jCondition2(i + 1, j) * normalRight[0] * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i + 1, j)));
				tempBeta = beta(i, j) * beta(i + 1, j) * (abs(levelSet(i, j)) + abs(levelSet(i + 1, j))) / (beta(i + 1, j) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i + 1, j)));
				fR = -tempBeta*aGamma / (grid.dx*grid.dx) - tempBeta*bGamma*theta / (beta(i + 1, j) * grid.dx);
				//fR.values[matIndex] = -tempBeta*aGamma / (grid.dx*grid.dx) - tempBeta*bGamma*theta / (beta(i + 1, j) * grid.dx);
				//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
				//cout<<i+1<<" " <<i+1+1 <<endl;
				//cout << i<< " "<< fR[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else
			{
				fR = 0;
			}

			if (levelSet(i, j - 1) > 0 && levelSet(i, j) <= 0)
			{
				theta = abs(levelSet(i, j - 1)) / (abs(levelSet(i, j - 1)) + abs(levelSet(i, j)));
				aGamma = (jCondition1(i, j - 1) * abs(levelSet(i, j)) + jCondition1(i, j) * abs(levelSet(i, j - 1))) / (abs(levelSet(i, j - 1)) + abs(levelSet(i, j)));
				bGamma = (jCondition2(i, j - 1) * normalBottom[1] * abs(levelSet(i, j)) + jCondition2(i, j) * normalCenter[1] * abs(levelSet(i, j - 1))) / (abs(levelSet(i, j - 1)) + abs(levelSet(i, j)));
				tempBeta = beta(i, j - 1) * beta(i, j) * (abs(levelSet(i, j - 1)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i, j - 1)) + beta(i, j - 1) * abs(levelSet(i, j)));
				fB = tempBeta*aGamma / (grid.dy*grid.dy) - tempBeta*bGamma*theta / (beta(i, j - 1) * grid.dy);
				//fB.values[matIndex] = tempBeta*aGamma / (grid.dy*grid.dy) - tempBeta*bGamma*theta / (beta(i, j - 1) * grid.dy);
				//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
				//cout<<i<<" " <<i+1 <<endl;
				//cout << i<< " "<< fL[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else if (levelSet(i, j - 1) <= 0 && levelSet(i, j) > 0)
			{
				theta = abs(levelSet(i, j - 1)) / (abs(levelSet(i, j - 1)) + abs(levelSet(i, j)));
				aGamma = (jCondition1(i, j - 1) * abs(levelSet(i, j)) + jCondition1(i, j) * abs(levelSet(i, j - 1))) / (abs(levelSet(i, j - 1)) + abs(levelSet(i, j)));
				bGamma = (jCondition2(i, j - 1) * normalBottom[1] * abs(levelSet(i, j)) + jCondition2(i, j) * normalCenter[1] * abs(levelSet(i, j - 1))) / (abs(levelSet(i, j - 1)) + abs(levelSet(i, j)));
				tempBeta = beta(i, j - 1) * beta(i, j) * (abs(levelSet(i, j - 1)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i, j - 1)) + beta(i, j - 1) * abs(levelSet(i, j)));
				fB = -tempBeta*aGamma / (grid.dy*grid.dy) + tempBeta*bGamma*theta / (beta(i, j - 1) * grid.dy);
				//fB.values[matIndex] = -tempBeta*aGamma / (grid.dy*grid.dy) + tempBeta*bGamma*theta / (beta(i, j - 1) * grid.dy);
				//cout<<levelSet[i]<<" "<<levelSet[i+1]<<endl;
				//cout<<i+1<<" " <<j+1 <<endl;
				//cout <<fB[matIndex]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else
			{
				fB = 0;
			}

			if (levelSet(i, j) <= 0 && levelSet(i, j + 1) > 0)
			{
				theta = abs(levelSet(i, j + 1)) / (abs(levelSet(i, j)) + abs(levelSet(i, j + 1)));
				aGamma = (jCondition1(i, j) * abs(levelSet(i, j + 1)) + jCondition1(i, j + 1) * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i, j + 1)));
				bGamma = (jCondition2(i, j) * normalCenter[1] * abs(levelSet(i, j + 1)) + jCondition2(i, j + 1) * normalTop[1] * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i, j + 1)));
				tempBeta = beta(i, j) * beta(i, j + 1) * (abs(levelSet(i, j)) + abs(levelSet(i, j + 1))) / (beta(i, j + 1) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i, j + 1)));
				fT = tempBeta*aGamma / (grid.dy*grid.dy) + tempBeta*bGamma*theta / (beta(i, j + 1) * grid.dy);
				//fT.values[matIndex] = tempBeta*aGamma / (grid.dy*grid.dy) + tempBeta*bGamma*theta / (beta(i, j + 1) * grid.dy);
				//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
				//cout<<i+1<<" " <<j+1 <<endl;
				//cout << fT[matIndex]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;

			}
			else if (levelSet(i, j) > 0 && levelSet(i, j + 1) <= 0)
			{
				theta = abs(levelSet(i, j + 1)) / (abs(levelSet(i, j)) + abs(levelSet(i, j + 1)));
				aGamma = (jCondition1(i, j) * abs(levelSet(i, j + 1)) + jCondition1(i, j + 1) * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i, j + 1)));
				bGamma = (jCondition2(i, j) * normalCenter[1] * abs(levelSet(i, j + 1)) + jCondition2(i, j + 1) * normalTop[1] * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i, j + 1)));
				tempBeta = beta(i, j) * beta(i, j + 1) * (abs(levelSet(i, j)) + abs(levelSet(i, j + 1))) / (beta(i, j + 1) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i, j + 1)));
				fT = -tempBeta*aGamma / (grid.dy*grid.dy) - tempBeta*bGamma*theta / (beta(i, j + 1) * grid.dy);
				//fT.values[matIndex] = -tempBeta*aGamma / (grid.dy*grid.dy) - tempBeta*bGamma*theta / (beta(i, j + 1) * grid.dy);
				//cout<<levelSet[i+1]<<" "<<levelSet[i+1+1]<<endl;
				//cout<<i+1<<" " <<j+1 <<endl;
				//cout <<fT[matIndex]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else
			{
				fT = 0;
			}

			poissonVector.values[indexInner(i, j)] = -f(i, j) - fR - fL - fB - fT;
		}

	}

	//ofstream ffl, ffr, ffb, fft, bbb, fff;
	//ffl.open("E:\Data/fl.txt");
	//ffr.open("E:\Data/fr.txt");
	//ffb.open("E:\Data/fb.txt");
	//fft.open("E:\Data/ft.txt");
	//bbb.open("E:\Data/b.txt");
	//fff.open("E:\Data/f.txt");

	//for (int j = 0; j < grid.numMatY; j++)
	//{
	//	for (int i = 0; i < grid.numMatX; i++)
	//	{
	//		ffl<<fL[i+j*grid.numMatX]<< " ";
	//		ffr<<fB[i+j*grid.numMatX]<< " ";
	//		ffb<<fB[i+j*grid.numMatX]<< " ";
	//		fft<<fT[i+j*grid.numMatX]<< " ";
	//		bbb<<-b[i+j*grid.numMatX]<< " ";
	//		fff<<f[i+j*grid.numMatX]<< " ";
	//	}
	//		ffl<< endl;
	//		ffr<< endl;
	//		ffb<< endl;
	//		fft<< endl;
	//		bbb<< endl;
	//		fff<< endl;
	//}
	//ffl.close();
	//ffr.close();
	//ffb.close();
	//fft.close();
	//bbb.close();
	//fff.close();	
}

inline void PoissonSolver::SolvePoissonJumpCondi(int example, const Grid2D& ipGrid)
{
	grid = ipGrid;
	solution = Field2D<double>(grid);

	innerGrid = Grid2D(ipGrid.xMin + ipGrid.dx, ipGrid.xMax - ipGrid.dx, 1, ipGrid.iRes - 2, ipGrid.yMin + ipGrid.dy, ipGrid.yMax - ipGrid.dy, 1, ipGrid.jRes - 2);
	VectorND<double> innerSolution(innerGrid.iRes*innerGrid.jRes);
	Field2D<double> beta(grid);
	Field2D<double> f(innerGrid);
	Field2D<double> jCondition1(grid);
	Field2D<double> jCondition2(grid);
	//double leftBdry, rightBdry;

	poissonMatrix = Array2D<double>(1, innerGrid.iRes*innerGrid.jRes, 1, innerGrid.iRes*innerGrid.jRes);
	poissonVector = VectorND<double>(poissonMatrix.iRes);

	GenerateJumpCondi(example, beta, f, jCondition1, jCondition2);

	GeneratePoissonMatrixJumpCondi(beta, f, jCondition1, jCondition2);

	GeneratePoissonVectorJumpCondi(beta, f, jCondition1, jCondition2);


	//ofstream qwer;
	//qwer.open("D:\Data/b1.txt", ios::binary);
	//for (int i = 0; i < grid.numMatX*grid.numMatY; i++)
	//{
	//	qwer << poissonVector[i] << endl;
	//}
	//qwer.close();


	//ofstream asdf;
	//asdf.open("D:\Data/A1.txt", ios::binary);
	//for (int i = innerGrid.iStart; i <= innerGrid.iEnd; i++)
	//{
	//	for (int j = innerGrid.jStart; j < innerGrid.jEnd; j++)
	//	{
	//		asdf << poissonMatrix[i*grid.numMatX*grid.numMatY + j] << " ";
	//	}
	//	asdf << endl;
	//}
	//asdf.close();


	poissonCSR = CSR<double>(poissonMatrix);
	innerSolution = CG<double>(poissonCSR, poissonVector);

	for (int i = innerGrid.iStart; i <= innerGrid.iEnd; i++)
	{
		for (int j = innerGrid.jStart; j < innerGrid.jEnd; j++)
		{
			solution(i, j) = innerSolution[indexInner(i, j)];
		}
	}


	OutputResult();

}

inline void PoissonSolver::OutputResult()
{
	//clock_t before;
	//double  result;
	//before = clock();

	ofstream solutionFile;
	solutionFile.open("D:\\Data/poisson.txt", ios::binary);

	for (int i = 0; i < grid.iRes; i++)
	{
		for (int j = 0; j < grid.jRes; j++)
		{
			solutionFile << i << " " << j << " " << grid(i, j) << " " << solution(i, j) << endl;
		}
	}

	solutionFile.close();

	//result = (double)(clock() - before) / CLOCKS_PER_SEC;
	//cout << "binary : " << result << "\n";
	////printf("걸린시간은 %5.2f 입니다.\n", result);
}
