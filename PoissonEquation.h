#pragma once

#include <iostream>
#include <time.h>

#include "CSR.h"
#include "GridInfomation.h"
#include "LevelSet.h"

#include "LinearSolver.h"
//#ifndef PoissonEquation_H
//#define PoissonEquation_H



class PoissonEquationSolver
{
public:
	double* poissonMatrix;
	double* poissonVector;
	double* beta;
	double* f;
	double* jCondition1;
	double* jCondition2;
	double* solution;
	double* tempSol;
	double leftBdry, rightBdry;

	csr sparsePoissonMatrix;
	GridInfo grid;
	LevelSet levelSet;

	PoissonEquationSolver();
	~PoissonEquationSolver();

	PoissonEquationSolver(GridInfo inputGrid);

	int index(int i, int j);
	int index(int i, int j, int k);
	int indexVec(int i, int j);
	int indexVec(int i, int j, int k);
	//int indexMat(int i, int j);
	//int indexMat(int i, int j, int k);

	// A Boundary Condition Capturing Method
	void generateJumpCondi(int example);
	void generatePoissonMatrixJumpCondi();
	void generatePoissonVectorJumpCondi();
	void solvePoissonEquationJumpCondi(int example);
	void outputResult();
private:

};

PoissonEquationSolver::PoissonEquationSolver()
{
	if (grid.dimension == 1)
	{
		poissonMatrix = new double[grid.numMatX*grid.numMatX];
		poissonVector = new double[grid.numMatX];
		beta = new double[grid.numX];
		f = new double[grid.numMatX];
		jCondition1 = new double[grid.numX];
		jCondition2 = new double[grid.numX];
		solution = new double[grid.numX];

		for (int i = 0; i < grid.numMatX*grid.numMatX; i++)
		{
			poissonMatrix[i] = 0;
		}
		for (int i = 0; i < grid.numMatX; i++)
		{
			poissonVector[i] = 0;
			f[i] = 0;
		}
		for (int i = 0; i < grid.numX; i++)
		{
			beta[i] = 0;
			jCondition1[i] = 0;
			jCondition2[i] = 0;
			solution[i] = 0;
		}
	}
	if (grid.dimension == 2)
	{
		poissonMatrix = new double[grid.numMatX*grid.numMatX*grid.numMatY*grid.numMatY];
		poissonVector = new double[grid.numMatX*grid.numMatY];
		beta = new double[grid.numX*grid.numY];
		f = new double[grid.numMatX*grid.numMatY];
		jCondition1 = new double[grid.numX*grid.numY];
		jCondition2 = new double[grid.numX*grid.numY];
		solution = new double[grid.numX*grid.numY];
		for (int i = 0; i < grid.numMatX*grid.numMatX*grid.numMatY*grid.numMatY; i++)
		{
			poissonMatrix[i] = 0;
		}
		for (int i = 0; i < grid.numMatX*grid.numMatY; i++)
		{
			poissonVector[i] = 0;
			f[i] = 0;
		}
		for (int i = 0; i < grid.numX*grid.numY; i++)
		{
			beta[i] = 0;
			jCondition1[i] = 0;
			jCondition2[i] = 0;
			solution[i] = 0;
		}
	}
}

PoissonEquationSolver::~PoissonEquationSolver()
{
	delete[] poissonMatrix, poissonVector, beta, f, jCondition1, jCondition2, solution, tempSol;
}



PoissonEquationSolver::PoissonEquationSolver(GridInfo inputGrid)
{
	grid = inputGrid;
	if (grid.dimension == 1)
	{
		poissonMatrix = new double[grid.numMatX*grid.numMatX];
		poissonVector = new double[grid.numMatX];
		beta = new double[grid.numX];
		f = new double[grid.numMatX];
		jCondition1 = new double[grid.numX];
		jCondition2 = new double[grid.numX];
		solution = new double[grid.numX];

		for (int i = 0; i < grid.numMatX*grid.numMatX; i++)
		{
			poissonMatrix[i] = 0;
		}
		for (int i = 0; i < grid.numMatX; i++)
		{
			poissonVector[i] = 0;
			f[i] = 0;
		}
		for (int i = 0; i < grid.numX; i++)
		{
			beta[i] = 0;
			jCondition1[i] = 0;
			jCondition2[i] = 0;
			solution[i] = 0;
		}
	}
	if (grid.dimension == 2)
	{
		poissonMatrix = new double[grid.numMatX*grid.numMatX*grid.numMatY*grid.numMatY];
		poissonVector = new double[grid.numMatX*grid.numMatY];
		beta = new double[grid.numX*grid.numY];
		f = new double[grid.numMatX*grid.numMatY];
		jCondition1 = new double[grid.numX*grid.numY];
		jCondition2 = new double[grid.numX*grid.numY];
		solution = new double[grid.numX*grid.numY];
		for (int i = 0; i < grid.numMatX*grid.numMatX*grid.numMatY*grid.numMatY; i++)
		{
			poissonMatrix[i] = 0;
		}
		for (int i = 0; i < grid.numMatX*grid.numMatY; i++)
		{
			poissonVector[i] = 0;
			f[i] = 0;
		}
		for (int i = 0; i < grid.numX*grid.numY; i++)
		{
			beta[i] = 0;
			jCondition1[i] = 0;
			jCondition2[i] = 0;
			solution[i] = 0;
		}
	}
}

inline int PoissonEquationSolver::index(int i, int j)
{
	return i + j*grid.numX;
}

inline int PoissonEquationSolver::index(int i, int j, int k)
{
	return i + j*grid.numX + k*grid.numX*grid.numY;
}

inline int PoissonEquationSolver::indexVec(int i, int j)
{
	return i + j*grid.numMatX;
}

inline int PoissonEquationSolver::indexVec(int i, int j, int k)
{
	return i + j*grid.numMatX + k*grid.numMatX*grid.numMatY;
}


inline void PoissonEquationSolver::generateJumpCondi(int example)
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

		levelSet = LevelSet(grid);
		for (int i = 0; i < grid.numX; i++)
		{
			levelSet.phi[i] = abs(grid.x[i] - 0.45) - 0.15 - grid.deltaX / 2;
		}

		leftBdry = 0, rightBdry = 0;

		for (int i = 0; i < grid.numX; i++)
		{
			if (levelSet.phi[i] <= 0)
			{
				beta[i] = 2;
			}
			else
			{
				beta[i] = 1;
			}
		}

		for (int i = 0; i < grid.numMatX; i++)
		{
			if (levelSet.phi[i + 1] <= 0)
			{
				f[i] = (8 * grid.x[i + 1] * grid.x[i + 1] - 4)*exp(-grid.x[i + 1] * grid.x[i + 1]);
				poissonVector[i] = f[i];
				//cout<<i<<" "<<f[i]<<endl;
			}
			//else
			//{
			//	f[i] = 0;
			//}
		}

		//for (int i = 0; i < grid.numMatX; i++)
		//{
		//	poissonVector[i] = f[i];
		//}

		//for (int i = 0; i < grid.numX; i++)
		//{
		//	jCondition1[i] = 0;
		//	jCondition2[i] = 0;
		//}

		jCondition1[29] = -exp(-0.09);
		jCondition1[30] = -exp(-0.09);
		jCondition2[29] = -1.2*exp(-0.09);
		jCondition2[30] = -1.2*exp(-0.09);

		jCondition1[60] = -exp(-0.36);
		jCondition1[61] = -exp(-0.36);
		jCondition2[60] = 2.4*exp(-0.36);
		jCondition2[61] = 2.4*exp(-0.36);
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

		levelSet = LevelSet(grid);
		for (int j = 0; j < grid.numY; j++)
		{
			for (int i = 0; i < grid.numX; i++)
			{
				levelSet.phi[index(i, j)] = sqrt((grid.x[i] - 0.5)*(grid.x[i] - 0.5) + (grid.y[j] - 0.5)*(grid.y[j] - 0.5)) - 0.25 - grid.deltaX / 2;
			}
		}

		for (int j = 0; j < grid.numY; j++)
		{
			for (int i = 0; i < grid.numX; i++)
			{
				if (levelSet.phi[index(i, j)] <= 0)
				{
					beta[index(i, j)] = 2;
				}
				else
				{
					beta[index(i, j)] = 1;
				}
			}
		}

		for (int j = 0; j < grid.numMatY; j++)
		{
			for (int i = 0; i < grid.numMatX; i++)
			{
				if (levelSet.phi[index(i + 1, j + 1)] <= 0)
				{
					f[indexVec(i, j)] = 8 * (grid.x[i + 1] * grid.x[i + 1] + grid.y[j + 1] * grid.y[j + 1] - 1)*exp(-grid.x[i + 1] * grid.x[i + 1] - grid.y[j + 1] * grid.y[j + 1]);
				}
				//else
				//{
				//	f[indexVec(i,j)] = 0;
				//}
			}
		}


		int centIndex, rightIndex, topIndex;
		for (int j = 0; j < grid.numY - 1; j++)
		{
			for (int i = 0; i < grid.numX - 1; i++)
			{
				centIndex = index(i, j);
				rightIndex = index(i + 1, j);
				topIndex = index(i, j + 1);
				if ((levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0) || (levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0))
				{
					jCondition1[centIndex] = -exp(-grid.x[i] * grid.x[i] - grid.y[j] * grid.y[j]);
					jCondition1[rightIndex] = -exp(-grid.x[i + 1] * grid.x[i + 1] - grid.y[j] * grid.y[j]);
					jCondition2[centIndex] = 8 * (2 * grid.x[i] * grid.x[i] + 2 * grid.y[j] * grid.y[j] - grid.x[i] - grid.y[j])*exp(-grid.x[i] * grid.x[i] - grid.y[j] * grid.y[j]);
					jCondition2[rightIndex] = 8 * (2 * grid.x[i + 1] * grid.x[i + 1] + 2 * grid.y[j] * grid.y[j] - grid.x[i + 1] - grid.y[j])*exp(-grid.x[i + 1] * grid.x[i + 1] - grid.y[j] * grid.y[j]);
				}
				if ((levelSet.phi[centIndex] <= 0 && levelSet.phi[topIndex]>0) || (levelSet.phi[centIndex]>0 && levelSet.phi[topIndex] <= 0))
				{
					jCondition1[centIndex] = -exp(-grid.x[i] * grid.x[i] - grid.y[j] * grid.y[j]);;
					jCondition1[topIndex] = -exp(-grid.x[i] * grid.x[i] - grid.y[j + 1] * grid.y[j + 1]);;
					jCondition2[centIndex] = 8 * (2 * grid.x[i] * grid.x[i] + 2 * grid.y[j] * grid.y[j] - grid.x[i] - grid.y[j])*exp(-grid.x[i] * grid.x[i] - grid.y[j] * grid.y[j]);
					jCondition2[topIndex] = 8 * (2 * grid.x[i] * grid.x[i] + 2 * grid.y[j + 1] * grid.y[j + 1] - grid.x[i] - grid.y[j + 1])*exp(-grid.x[i] * grid.x[i] - grid.y[j + 1] * grid.y[j + 1]);
				}
			}
		}

	}
}

inline void PoissonEquationSolver::generatePoissonMatrixJumpCondi()
{
	double tempBeta = 0;
	if (grid.dimension == 1)
	{
		for (int i = 0; i < grid.numMatX; i++)
		{
			if (i>0 && i<grid.numMatX - 1)
			{
				if ((levelSet.phi[i + 1]>0 && levelSet.phi[i + 1 + 1] <= 0) || (levelSet.phi[i + 1] <= 0 && levelSet.phi[i + 1 + 1]>0))
				{
					tempBeta = beta[i + 1] * beta[i + 1 + 1] * (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1])) / (beta[i + 1 + 1] * abs(levelSet.phi[i + 1]) + beta[i + 1] * abs(levelSet.phi[i + 1 + 1]));
					poissonMatrix[i*grid.numMatX + i - 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2;
					poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
					poissonMatrix[i*grid.numMatX + i + 1] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
					//cout<<i<<" "<< tempBeta<<endl;
				}
				else if ((levelSet.phi[i]>0 && levelSet.phi[i + 1] <= 0) || (levelSet.phi[i] <= 0 && levelSet.phi[i + 1]>0))
				{
					tempBeta = beta[i] * beta[i + 1] * (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1])) / (beta[i + 1] * abs(levelSet.phi[i]) + beta[i] * abs(levelSet.phi[i + 1]));
					poissonMatrix[i*grid.numMatX + i - 1] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
					poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
					poissonMatrix[i*grid.numMatX + i + 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
					//cout<<i<<" "<< tempBeta<<endl;
				}
				else
				{
					poissonMatrix[i*grid.numMatX + i - 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2;
					poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
					poissonMatrix[i*grid.numMatX + i + 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
				}

			}
			else if (i == 0)
			{
				poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*beta[i] + 1 / (grid.deltaX*grid.deltaX)*(beta[i] + beta[i + 1]) / 2;
				poissonMatrix[i*grid.numMatX + i + 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i + 1 + 1]) / 2;
			}
			else
			{
				poissonMatrix[i*grid.numMatX + i - 1] = -1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2;
				poissonMatrix[i*grid.numMatX + i] = 1 / (grid.deltaX*grid.deltaX)*(beta[i + 1] + beta[i]) / 2 + 1 / (grid.deltaX*grid.deltaX)*beta[i + 1];
			}
		}
	}
	else if (grid.dimension == 2)
	{
		int matIndex, matLeftIndex, matRightIndex, matTopIndex, matBottomIndex;
		int centIndex, leftIndex, rightIndex, topIndex, bottomIndex;

		for (int j = 0; j < grid.numMatY; j++)
		{
			//poissonMatrix[i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*grid.numMatX*grid.numMatY + j*grid.numMatX ] = 0;
			//poissonMatrix[i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*(grid.numMatX*grid.numMatY + 1) ] = 0;
			for (int i = 0; i < grid.numMatX; i++)
			{
				matIndex = i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*(grid.numMatX*grid.numMatY + 1);
				matLeftIndex = i*grid.numMatX*grid.numMatY + i - 1 + j*grid.numMatX*(grid.numMatX*grid.numMatY + 1);
				matRightIndex = i*grid.numMatX*grid.numMatY + i + 1 + j*grid.numMatX*(grid.numMatX*grid.numMatY + 1);
				matBottomIndex = i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*grid.numMatX*grid.numMatY + (j - 1)*grid.numMatX;
				matTopIndex = i*grid.numMatX*grid.numMatY + i + j*grid.numMatX*grid.numMatX*grid.numMatY + (j + 1)*grid.numMatX;

				centIndex = i + 1 + (j + 1)*grid.numX;
				leftIndex = i + (j + 1)*grid.numX;
				rightIndex = i + 1 + 1 + (j + 1)*grid.numX;
				bottomIndex = i + 1 + j*grid.numX;
				topIndex = i + 1 + (j + 1 + 1)*grid.numX;

				if (j>0 && j<grid.numMatY - 1)
				{
					if ((levelSet.phi[centIndex]>0 && levelSet.phi[topIndex] <= 0) || (levelSet.phi[centIndex] <= 0 && levelSet.phi[topIndex]>0))
					{
						tempBeta = beta[centIndex] * beta[topIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex])) / (beta[topIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[topIndex]));
						poissonMatrix[matBottomIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(tempBeta);
						poissonMatrix[matTopIndex] = -1 / (grid.deltaY*grid.deltaY)*tempBeta;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
						//cout<<"";

					}
					else if ((levelSet.phi[bottomIndex]>0 && levelSet.phi[centIndex] <= 0) || (levelSet.phi[bottomIndex] <= 0 && levelSet.phi[centIndex]>0))
					{
						tempBeta = beta[bottomIndex] * beta[centIndex] * (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[bottomIndex]) + beta[bottomIndex] * abs(levelSet.phi[centIndex]));
						poissonMatrix[matBottomIndex] = -1 / (grid.deltaY*grid.deltaY)*tempBeta;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(tempBeta);
						poissonMatrix[matTopIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2;
						//cout<<i<<" "<< tempBeta<<endl;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matBottomIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
						poissonMatrix[matTopIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2;

						//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
						//cout<<"";
					}

					if (i>0 && i<grid.numMatX - 1)
					{
						if ((levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0) || (levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0))
						{
							tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else if ((levelSet.phi[leftIndex]>0 && levelSet.phi[centIndex] <= 0) || (levelSet.phi[leftIndex] <= 0 && levelSet.phi[centIndex]>0))
						{
							tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else
						{
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
					}
					else if (i == 0)
					{
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
						//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0<<endl;
						//cout<<"";
					}
				}
				else if (j == 0)
				{
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
					poissonMatrix[matTopIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2;
					//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matTopIndex] <<endl;
					//cout<<"";
					if (i>0 && i<grid.numMatX - 1)
					{
						if ((levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0) || (levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0))
						{
							tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else if ((levelSet.phi[leftIndex]>0 && levelSet.phi[centIndex] <= 0) || (levelSet.phi[leftIndex] <= 0 && levelSet.phi[centIndex]>0))
						{
							tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else
						{
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
					}
					else if (i == 0)
					{
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
						//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0<<endl;
						//cout<<"";
					}
				}
				else
				{
					poissonMatrix[matBottomIndex] = -1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
					poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[topIndex]) / 2 + 1 / (grid.deltaY*grid.deltaY)*(beta[centIndex] + beta[bottomIndex]) / 2;
					//cout<<i<<" " << j<<" "<<poissonMatrix[matBottomIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0 <<endl;
					//cout<<"";
					if (i>0 && i<grid.numMatX - 1)
					{
						if ((levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0) || (levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0))
						{
							tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else if ((levelSet.phi[leftIndex]>0 && levelSet.phi[centIndex] <= 0) || (levelSet.phi[leftIndex] <= 0 && levelSet.phi[centIndex]>0))
						{
							tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*tempBeta;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(tempBeta);
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" "<< tempBeta<<endl;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
						else
						{
							poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
							poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
							//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
							//cout<<"";
						}
					}
					else if (i == 0)
					{
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matRightIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2;
						//cout<<i<<" " << j<<" "<<0<<" "<< poissonMatrix[matIndex]<<" " <<poissonMatrix[matRightIndex] <<endl;
						//cout<<"";
					}
					else
					{
						poissonMatrix[matLeftIndex] = -1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						poissonMatrix[matIndex] = poissonMatrix[matIndex] + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[rightIndex]) / 2 + 1 / (grid.deltaX*grid.deltaX)*(beta[centIndex] + beta[leftIndex]) / 2;
						//cout<<i<<" " << j<<" "<<poissonMatrix[matLeftIndex]<<" "<< poissonMatrix[matIndex]<<" " <<0 <<endl;
						//cout<<"";
					}
				}
			}
		}
	}
}

inline void PoissonEquationSolver::generatePoissonVectorJumpCondi()
{
	double aGamma = 0;
	double bGamma = 0;
	double theta = 0;
	double tempBeta = 0;

	if (grid.dimension == 1)
	{
		double* fL = new double[grid.numMatX];
		double* fR = new double[grid.numMatX];

		double normalLeft, normalCenter, normalRight;

		for (int i = 0; i < grid.numMatX; i++)
		{
			normalLeft = levelSet.unitNormal(i);
			normalCenter = levelSet.unitNormal(i + 1);
			normalRight = levelSet.unitNormal(i + 1 + 1);
			if (levelSet.phi[i]>0 && levelSet.phi[i + 1] <= 0)
			{
				theta = abs(levelSet.phi[i]) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				aGamma = (jCondition1[i] * abs(levelSet.phi[i + 1]) + jCondition1[i + 1] * abs(levelSet.phi[i])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				bGamma = (jCondition2[i] * normalLeft*abs(levelSet.phi[i + 1]) + jCondition2[i + 1] * normalCenter*abs(levelSet.phi[i])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				tempBeta = beta[i] * beta[i + 1] * (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1])) / (beta[i + 1] * abs(levelSet.phi[i]) + beta[i] * abs(levelSet.phi[i + 1]));
				fL[i] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[i] * grid.deltaX);
				//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
				//cout<<i<<" " <<i+1 <<endl;
				//cout << i<< " "<< fL[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;
			}
			else if (levelSet.phi[i] <= 0 && levelSet.phi[i + 1]>0)
			{
				theta = abs(levelSet.phi[i]) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				aGamma = (jCondition1[i] * abs(levelSet.phi[i + 1]) + jCondition1[i + 1] * abs(levelSet.phi[i])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				bGamma = (jCondition2[i] * normalLeft*abs(levelSet.phi[i + 1]) + jCondition2[i + 1] * normalCenter*abs(levelSet.phi[i])) / (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1]));
				tempBeta = beta[i] * beta[i + 1] * (abs(levelSet.phi[i]) + abs(levelSet.phi[i + 1])) / (beta[i + 1] * abs(levelSet.phi[i]) + beta[i] * abs(levelSet.phi[i + 1]));
				fL[i] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[i] * grid.deltaX);
				//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
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
				fL[i] = 0;
			}

			if (levelSet.phi[i + 1] <= 0 && levelSet.phi[i + 1 + 1]>0)
			{
				theta = abs(levelSet.phi[i + 1 + 1]) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				aGamma = (jCondition1[i + 1] * abs(levelSet.phi[i + 1 + 1]) + jCondition1[i + 1 + 1] * abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				bGamma = (jCondition2[i + 1] * normalCenter*abs(levelSet.phi[i + 1 + 1]) + jCondition2[i + 1 + 1] * normalRight*abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				tempBeta = beta[i + 1] * beta[i + 1 + 1] * (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1])) / (beta[i + 1 + 1] * abs(levelSet.phi[i + 1]) + beta[i + 1] * abs(levelSet.phi[i + 1 + 1]));
				fR[i] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[i + 1 + 1] * grid.deltaX);
				//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
				//cout<<i+1<<" " <<i+1+1 <<endl;
				//cout << i<< " "<< fR[i]<<endl;
				//cout<<"theta"<<theta<<endl;
				//cout<<"aGamma : "<<aGamma<<endl;
				//cout<<"bGamma : "<<bGamma<<endl;
				//cout<<"tempBeta : "<<tempBeta<<endl;
				//cout<<endl;

			}
			else if (levelSet.phi[i + 1]>0 && levelSet.phi[i + 1 + 1] <= 0)
			{
				theta = abs(levelSet.phi[i + 1 + 1]) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				aGamma = (jCondition1[i + 1] * abs(levelSet.phi[i + 1 + 1]) + jCondition1[i + 1 + 1] * abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				bGamma = (jCondition2[i + 1] * normalCenter*abs(levelSet.phi[i + 1 + 1]) + jCondition2[i + 1 + 1] * normalRight*abs(levelSet.phi[i + 1])) / (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1]));
				tempBeta = beta[i + 1] * beta[i + 1 + 1] * (abs(levelSet.phi[i + 1]) + abs(levelSet.phi[i + 1 + 1])) / (beta[i + 1 + 1] * abs(levelSet.phi[i + 1]) + beta[i + 1] * abs(levelSet.phi[i + 1 + 1]));
				fR[i] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[i + 1 + 1] * grid.deltaX);
				//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
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
				fR[i] = 0;
			}

			poissonVector[i] = -f[i] - fR[i] - fL[i];

		}

		delete[] fR, fL;
	}
	else if (grid.dimension == 2)
	{
		double* fL = new double[grid.numMatX*grid.numMatY];
		double* fR = new double[grid.numMatX*grid.numMatY];
		double* fB = new double[grid.numMatX*grid.numMatY];
		double* fT = new double[grid.numMatX*grid.numMatY];

		double* normalLeft = new double[2];
		double* normalRight = new double[2];
		double* normalCenter = new double[2];
		double* normalBottom = new double[2];
		double* normalTop = new double[2];

		int matIndex;
		int centIndex, leftIndex, rightIndex, topIndex, bottomIndex;

		for (int j = 0; j < grid.numMatY; j++)
		{
			for (int i = 0; i < grid.numMatX; i++)
			{
				matIndex = indexVec(i, j);// i + j*grid.numMatX;

				centIndex = index(i + 1, j + 1);		// i + 1 + (j + 1)*grid.numX;
				leftIndex = index(i, j + 1);			// i + (j + 1)*grid.numX;
				rightIndex = index(i + 1 + 1, j + 1);	// i + 1 + 1 + (j + 1)*grid.numX;
				bottomIndex = index(i + 1, j);			// i + 1 + j*grid.numX;
				topIndex = index(i + 1, j + 1 + 1);	// i + 1 + (j + 1 + 1)*grid.numX;

				levelSet.unitNormal(i - 1, j, normalLeft);
				levelSet.unitNormal(i + 1, j, normalRight);
				levelSet.unitNormal(i, j, normalCenter);
				levelSet.unitNormal(i, j - 1, normalBottom);
				levelSet.unitNormal(i, j + 1, normalTop);

				if (levelSet.phi[leftIndex]>0 && levelSet.phi[centIndex] <= 0)
				{
					theta = abs(levelSet.phi[leftIndex]) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					aGamma = (jCondition1[leftIndex] * abs(levelSet.phi[centIndex]) + jCondition1[centIndex] * abs(levelSet.phi[leftIndex])) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					bGamma = (jCondition2[leftIndex] * normalLeft[0] * abs(levelSet.phi[centIndex]) + jCondition2[centIndex] * normalCenter[0] * abs(levelSet.phi[leftIndex])) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
					fL[matIndex] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[leftIndex] * grid.deltaX);
					//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
					//cout<<i<<" " <<i+1 <<endl;
					//cout << i<< " "<< fL[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else if (levelSet.phi[leftIndex] <= 0 && levelSet.phi[centIndex]>0)
				{
					theta = abs(levelSet.phi[leftIndex]) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					aGamma = (jCondition1[leftIndex] * abs(levelSet.phi[centIndex]) + jCondition1[centIndex] * abs(levelSet.phi[leftIndex])) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					bGamma = (jCondition2[leftIndex] * normalLeft[0] * abs(levelSet.phi[centIndex]) + jCondition2[centIndex] * normalCenter[0] * abs(levelSet.phi[leftIndex])) / (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex]));
					tempBeta = beta[leftIndex] * beta[centIndex] * (abs(levelSet.phi[leftIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[leftIndex]) + beta[leftIndex] * abs(levelSet.phi[centIndex]));
					fL[matIndex] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[leftIndex] * grid.deltaX);
					//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
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
					fL[matIndex] = 0;
				}

				if (levelSet.phi[centIndex] <= 0 && levelSet.phi[rightIndex]>0)
				{
					theta = abs(levelSet.phi[rightIndex]) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet.phi[rightIndex]) + jCondition1[rightIndex] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[0] * abs(levelSet.phi[rightIndex]) + jCondition2[rightIndex] * normalRight[0] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
					fR[matIndex] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[rightIndex] * grid.deltaX);
					//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
					//cout<<i+1<<" " <<i+1+1 <<endl;
					//cout << i<< " "<< fR[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;

				}
				else if (levelSet.phi[centIndex]>0 && levelSet.phi[rightIndex] <= 0)
				{
					theta = abs(levelSet.phi[rightIndex]) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet.phi[rightIndex]) + jCondition1[rightIndex] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[0] * abs(levelSet.phi[rightIndex]) + jCondition2[rightIndex] * normalRight[0] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex]));
					tempBeta = beta[centIndex] * beta[rightIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[rightIndex])) / (beta[rightIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[rightIndex]));
					fR[matIndex] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[rightIndex] * grid.deltaX);
					//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
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
					fR[matIndex] = 0;
				}

				if (levelSet.phi[bottomIndex]>0 && levelSet.phi[centIndex] <= 0)
				{
					theta = abs(levelSet.phi[bottomIndex]) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					aGamma = (jCondition1[bottomIndex] * abs(levelSet.phi[centIndex]) + jCondition1[centIndex] * abs(levelSet.phi[bottomIndex])) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					bGamma = (jCondition2[bottomIndex] * normalBottom[1] * abs(levelSet.phi[centIndex]) + jCondition2[centIndex] * normalCenter[1] * abs(levelSet.phi[bottomIndex])) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					tempBeta = beta[bottomIndex] * beta[centIndex] * (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[bottomIndex]) + beta[bottomIndex] * abs(levelSet.phi[centIndex]));
					fB[matIndex] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[bottomIndex] * grid.deltaX);
					//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
					//cout<<i<<" " <<i+1 <<endl;
					//cout << i<< " "<< fL[i]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;
				}
				else if (levelSet.phi[bottomIndex] <= 0 && levelSet.phi[centIndex]>0)
				{
					theta = abs(levelSet.phi[bottomIndex]) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					aGamma = (jCondition1[bottomIndex] * abs(levelSet.phi[centIndex]) + jCondition1[centIndex] * abs(levelSet.phi[bottomIndex])) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					bGamma = (jCondition2[bottomIndex] * normalBottom[1] * abs(levelSet.phi[centIndex]) + jCondition2[centIndex] * normalCenter[1] * abs(levelSet.phi[bottomIndex])) / (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex]));
					tempBeta = beta[bottomIndex] * beta[centIndex] * (abs(levelSet.phi[bottomIndex]) + abs(levelSet.phi[centIndex])) / (beta[centIndex] * abs(levelSet.phi[bottomIndex]) + beta[bottomIndex] * abs(levelSet.phi[centIndex]));
					fB[matIndex] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[bottomIndex] * grid.deltaX);
					//cout<<levelSet.phi[i]<<" "<<levelSet.phi[i+1]<<endl;
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
					fB[matIndex] = 0;
				}

				if (levelSet.phi[centIndex] <= 0 && levelSet.phi[topIndex]>0)
				{
					theta = abs(levelSet.phi[topIndex]) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet.phi[topIndex]) + jCondition1[topIndex] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[1] * abs(levelSet.phi[topIndex]) + jCondition2[topIndex] * normalTop[1] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					tempBeta = beta[centIndex] * beta[topIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex])) / (beta[topIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[topIndex]));
					fT[matIndex] = tempBeta*aGamma / (grid.deltaX*grid.deltaX) + tempBeta*bGamma*theta / (beta[topIndex] * grid.deltaX);
					//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
					//cout<<i+1<<" " <<j+1 <<endl;
					//cout << fT[matIndex]<<endl;
					//cout<<"theta"<<theta<<endl;
					//cout<<"aGamma : "<<aGamma<<endl;
					//cout<<"bGamma : "<<bGamma<<endl;
					//cout<<"tempBeta : "<<tempBeta<<endl;
					//cout<<endl;

				}
				else if (levelSet.phi[centIndex]>0 && levelSet.phi[topIndex] <= 0)
				{
					theta = abs(levelSet.phi[topIndex]) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					aGamma = (jCondition1[centIndex] * abs(levelSet.phi[topIndex]) + jCondition1[topIndex] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					bGamma = (jCondition2[centIndex] * normalCenter[1] * abs(levelSet.phi[topIndex]) + jCondition2[topIndex] * normalTop[1] * abs(levelSet.phi[centIndex])) / (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex]));
					tempBeta = beta[centIndex] * beta[topIndex] * (abs(levelSet.phi[centIndex]) + abs(levelSet.phi[topIndex])) / (beta[topIndex] * abs(levelSet.phi[centIndex]) + beta[centIndex] * abs(levelSet.phi[topIndex]));
					fT[matIndex] = -tempBeta*aGamma / (grid.deltaX*grid.deltaX) - tempBeta*bGamma*theta / (beta[topIndex] * grid.deltaX);
					//cout<<levelSet.phi[i+1]<<" "<<levelSet.phi[i+1+1]<<endl;
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
					fT[matIndex] = 0;
				}

				poissonVector[matIndex] = -f[matIndex] - fR[matIndex] - fL[matIndex] - fB[matIndex] - fT[matIndex];
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

		delete[] fR, fL, fT, fB, normalBottom, normalLeft, normalRight, normalTop, normalCenter;
	}
}

inline void PoissonEquationSolver::solvePoissonEquationJumpCondi(int example)
{
	generateJumpCondi(example);

	generatePoissonMatrixJumpCondi();

	generatePoissonVectorJumpCondi();


	//ofstream qwer;
	//qwer.open("D:\Data/b1.txt", ios::binary);
	//for (int i = 0; i < grid.numMatX*grid.numMatY; i++)
	//{
	//	qwer << poissonVector[i] << endl;
	//}
	//qwer.close();


	//ofstream asdf;
	//asdf.open("D:\Data/A1.txt", ios::binary);
	//for (int i = 0; i < grid.numMatX*grid.numMatY; i++)
	//{
	//	for (int j = 0; j < grid.numMatX*grid.numMatY; j++)
	//	{
	//		asdf << poissonMatrix[i*grid.numMatX*grid.numMatY + j] << " ";
	//	}
	//	asdf << endl;
	//}
	//asdf.close();




	if (grid.dimension == 1)
	{
		tempSol = CG(grid.numMatX, poissonMatrix, poissonVector);

		solution[0] = leftBdry; solution[grid.numX - 1] = rightBdry;
		for (int i = 0; i < grid.numMatX; i++)
		{
			solution[i + 1] = tempSol[i];
			//cout<<i<<" "<<sol[i]<<endl;
		}
	}
	else if (grid.dimension == 2)
	{
		sparsePoissonMatrix = csr(grid.numMatX*grid.numMatY, grid.numMatX*grid.numMatY, poissonMatrix);
		tempSol = CG(sparsePoissonMatrix, poissonVector);

		for (int i = 0; i < grid.numMatX; i++)
		{
			for (int j = 0; j < grid.numMatY; j++)
			{
				solution[i + 1 + (j + 1)*grid.numX] = tempSol[i + j*grid.numMatX];
			}
		}
	}


	outputResult();

}

inline void PoissonEquationSolver::outputResult()
{
	//clock_t before;
	//double  result;
	//before = clock();

	ofstream solutionFile;
	solutionFile.open("D:\\Data/poisson.txt", ios::binary);
	if (grid.dimension == 1)
	{
		for (int i = 0; i < grid.numX; i++)
		{
			solutionFile << i << " " << grid.x[i] << " " << solution[i] << endl;
		}
	}
	if (grid.dimension == 2)
	{
		for (int i = 0; i < grid.numX; i++)
		{
			for (int j = 0; j < grid.numY; j++)
			{
				solutionFile << i << " " << j << " " << grid.x[i] << " " << grid.y[j] << " " << solution[i + j*grid.numX] << endl;
			}
		}
	}

	solutionFile.close();

	//result = (double)(clock() - before) / CLOCKS_PER_SEC;
	//cout << "binary : " << result << "\n";
	////printf("걸린시간은 %5.2f 입니다.\n", result);
}



//#endif // !PoissonEquation

