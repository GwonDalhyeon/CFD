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

	FD solution;
	VectorND<double> innerSolution;

	CSR<double> poissonCSR;
	LS levelSet;

	PoissonSolver();
	~PoissonSolver();

	int index(int i, int j);
	int indexInner(int i, int j);
	//int indexVec(int i, int j);
	int indexMat(int i, int j);

	// A Boundary Condition Capturing Method
	PoissonSolver(const Grid2D& ipGrid, const LS& ipLevelSet, const FD& ipBeta, const FD& ipF, const FD& ipjCondition1, const FD& ipjCondition2);

	void GenerateJumpCondi(int example, FD& beta, FD& f, FD& jCondition1, FD& jCondition2);
	void GeneratePoissonMatrixJumpCondi(const FD& ipBeta, const FD& ipF, const FD& ipjCondition1, const FD& ipjCondition2);
	void GeneratePoissonVectorJumpCondi(const FD& ipBeta, const FD& ipF, const FD& ipjCondition1, const FD& ipjCondition2);
	void SolvePoissonJumpCondi(int example);


private:

};

//#endif // !PoissonSolver_H



PoissonSolver::PoissonSolver()
{
}


PoissonSolver::~PoissonSolver()
{
}

PoissonSolver::PoissonSolver(const Grid2D& ipGrid, const LS& ipLevelSet, const FD& ipBeta, const FD& ipF, const FD& ipjCondition1, const FD& ipjCondition2)
{
	grid = ipGrid;
	solution = FD(grid);
	levelSet = ipLevelSet;

	Grid2D innerGrid = Grid2D(ipGrid.xMin + ipGrid.dx, ipGrid.xMax - ipGrid.dx, 1, ipGrid.iRes - 2, ipGrid.yMin + ipGrid.dy, ipGrid.yMax - ipGrid.dy, 1, ipGrid.jRes - 2);
	VectorND<double> innerSolution(innerGrid.iRes*innerGrid.jRes);
	FD beta = ipBeta;;
	FD f = ipF;
	FD jCondition1 = ipjCondition1;
	FD jCondition2 = ipjCondition2;
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


inline void PoissonSolver::GenerateJumpCondi(int example, FD& beta, FD& f, FD& jCondition1, FD& jCondition2)
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

		levelSet = LS(grid);
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

inline void PoissonSolver::GeneratePoissonMatrixJumpCondi(const FD& beta, const FD& f, const FD& jCondition1, const FD& jCondition2)
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

inline void PoissonSolver::GeneratePoissonVectorJumpCondi(const FD& beta, const FD& f, const FD& jCondition1, const FD& jCondition2)
{
	double aGamma = 0;
	double bGamma = 0;
	double theta = 0;
	double tempBeta = 0;

	double fL, fR, fB, fT;

	VT normalLeft;
	VT normalRight;
	VT normalCenter;
	VT normalBottom;
	VT normalTop;

	levelSet.ComputeUnitNormal();

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
			}
			else if (levelSet(i - 1, j) <= 0 && levelSet(i, j) > 0)
			{
				theta = abs(levelSet(i - 1, j)) / (abs(levelSet(i - 1, j)) + abs(levelSet(i, j)));
				aGamma = (jCondition1(i - 1, j) * abs(levelSet(i, j)) + jCondition1(i, j) * abs(levelSet(i - 1, j))) / (abs(levelSet(i - 1, j)) + abs(levelSet(i, j)));
				bGamma = (jCondition2(i - 1, j) * normalLeft[0] * abs(levelSet(i, j)) + jCondition2(i, j) * normalCenter[0] * abs(levelSet(i - 1, j))) / (abs(levelSet(i - 1, j)) + abs(levelSet(i, j)));
				tempBeta = beta(i - 1, j) * beta(i, j) * (abs(levelSet(i - 1, j)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i - 1, j)) + beta(i - 1, j) * abs(levelSet(i, j)));
				fL = -tempBeta*aGamma / (grid.dx*grid.dx) + tempBeta*bGamma*theta / (beta(i - 1, j) * grid.dx);
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
			}
			else if (levelSet(i, j) > 0 && levelSet(i + 1, j) <= 0)
			{
				theta = abs(levelSet(i + 1, j)) / (abs(levelSet(i, j)) + abs(levelSet(i + 1, j)));
				aGamma = (jCondition1(i, j) * abs(levelSet(i + 1, j)) + jCondition1(i + 1, j) * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i + 1, j)));
				bGamma = (jCondition2(i, j) * normalCenter[0] * abs(levelSet(i + 1, j)) + jCondition2(i + 1, j) * normalRight[0] * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i + 1, j)));
				tempBeta = beta(i, j) * beta(i + 1, j) * (abs(levelSet(i, j)) + abs(levelSet(i + 1, j))) / (beta(i + 1, j) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i + 1, j)));
				fR = -tempBeta*aGamma / (grid.dx*grid.dx) - tempBeta*bGamma*theta / (beta(i + 1, j) * grid.dx);
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
			}
			else if (levelSet(i, j - 1) <= 0 && levelSet(i, j) > 0)
			{
				theta = abs(levelSet(i, j - 1)) / (abs(levelSet(i, j - 1)) + abs(levelSet(i, j)));
				aGamma = (jCondition1(i, j - 1) * abs(levelSet(i, j)) + jCondition1(i, j) * abs(levelSet(i, j - 1))) / (abs(levelSet(i, j - 1)) + abs(levelSet(i, j)));
				bGamma = (jCondition2(i, j - 1) * normalBottom[1] * abs(levelSet(i, j)) + jCondition2(i, j) * normalCenter[1] * abs(levelSet(i, j - 1))) / (abs(levelSet(i, j - 1)) + abs(levelSet(i, j)));
				tempBeta = beta(i, j - 1) * beta(i, j) * (abs(levelSet(i, j - 1)) + abs(levelSet(i, j))) / (beta(i, j) * abs(levelSet(i, j - 1)) + beta(i, j - 1) * abs(levelSet(i, j)));
				fB = -tempBeta*aGamma / (grid.dy*grid.dy) + tempBeta*bGamma*theta / (beta(i, j - 1) * grid.dy);
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
			}
			else if (levelSet(i, j) > 0 && levelSet(i, j + 1) <= 0)
			{
				theta = abs(levelSet(i, j + 1)) / (abs(levelSet(i, j)) + abs(levelSet(i, j + 1)));
				aGamma = (jCondition1(i, j) * abs(levelSet(i, j + 1)) + jCondition1(i, j + 1) * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i, j + 1)));
				bGamma = (jCondition2(i, j) * normalCenter[1] * abs(levelSet(i, j + 1)) + jCondition2(i, j + 1) * normalTop[1] * abs(levelSet(i, j))) / (abs(levelSet(i, j)) + abs(levelSet(i, j + 1)));
				tempBeta = beta(i, j) * beta(i, j + 1) * (abs(levelSet(i, j)) + abs(levelSet(i, j + 1))) / (beta(i, j + 1) * abs(levelSet(i, j)) + beta(i, j) * abs(levelSet(i, j + 1)));
				fT = -tempBeta*aGamma / (grid.dy*grid.dy) - tempBeta*bGamma*theta / (beta(i, j + 1) * grid.dy);
			}
			else
			{
				fT = 0;
			}

			poissonVector.values[indexInner(i, j)] = -f(i, j) - fR - fL - fB - fT;
		}

	}
}

inline void PoissonSolver::SolvePoissonJumpCondi(int example)
{
	grid = Grid2D(0, 1, 101, 0, 1, 101);
	grid.Variable();
	solution = FD(grid);

	innerGrid = Grid2D(grid.xMin + grid.dx, grid.xMax - grid.dx, 1, grid.iRes - 2, grid.yMin + grid.dy, grid.yMax - grid.dy, 1, grid.jRes - 2);
	VectorND<double> innerSolution(innerGrid.iRes*innerGrid.jRes);
	FD beta(grid);
	FD f(innerGrid);
	FD jCondition1(grid);
	FD jCondition2(grid);
	//double leftBdry, rightBdry;

	poissonMatrix = Array2D<double>(1, innerGrid.iRes*innerGrid.jRes, 1, innerGrid.iRes*innerGrid.jRes);
	poissonVector = VectorND<double>(poissonMatrix.iRes);

	GenerateJumpCondi(example, beta, f, jCondition1, jCondition2);

	GeneratePoissonMatrixJumpCondi(beta, f, jCondition1, jCondition2);

	GeneratePoissonVectorJumpCondi(beta, f, jCondition1, jCondition2);

	//poissonVector.Variable("poissonVector");
	//poissonMatrix.Variable("poissonMatrix");


	poissonCSR = CSR<double>(poissonMatrix);
	CGSolver::Solver(poissonCSR, poissonVector, innerSolution);

	#pragma omp parallel for
	for (int i = innerGrid.iStart; i <= innerGrid.iEnd; i++)
	{
		for (int j = innerGrid.jStart; j < innerGrid.jEnd; j++)
		{
			solution(i, j) = innerSolution[indexInner(i, j)];
		}
	}
	
	solution.Variable("poisson");
	MATLAB.Command("figure('units','normalized','outerposition',[0 0 1/2 1])");
	MATLAB.Command("x = reshape(X,1,size(X,1)*size(X,2))");
	MATLAB.Command("y = reshape(Y,1,size(X,1)*size(X,2))");
	MATLAB.Command("P = reshape(poisson,1,size(X,1)*size(X,2))");
	MATLAB.Command("plot3(x,y,P,'bo');");
	MATLAB.Command("grid on");
}
