#pragma once

#include "CommonDef.h"
#include "Vector2D.h"
#include "VectorND.h"
#include "Grid2D.h"
#include "Array2D.h"


// Declare MATLAB variable.
// TT<Vector2D<tt>> type structure to MATLAB
template<class TT>
void VecND2DVariable(const char * varName, const VectorND<Vector2D<TT>>& ipVec)
{
	int rowNum = ipVec.iLength;
	//double* pointxy = new double[10000];
	double* pointxy = new double[rowNum * 2];
	mxArray* dataArray = mxCreateDoubleMatrix(rowNum, 2, mxREAL);

#pragma omp parallel for
	for (int i = ipVec.iStart; i <= ipVec.iEnd; i++)
	{
		pointxy[i - ipVec.iStart] = double(ipVec(i).x);;
		pointxy[i - ipVec.iStart + rowNum] = double(ipVec(i).y);
	}

	memcpy((void*)mxGetPr(dataArray), (void*)pointxy, sizeof(double) * rowNum * 2);
	engPutVariable(MATLAB.ME, varName, dataArray);
}

template<class TT>
void ArrayVec2DVariable(const char * varName, const Array2D<Vector2D<TT>>& ipArray)
{
	Array2D<double> tempArray1(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);
	Array2D<double> tempArray2(ipArray.iStart, ipArray.iRes, ipArray.jStart, ipArray.jRes);
	
#pragma omp parallel for
	for (int i = ipArray.iStart; i <= ipArray.iEnd; i++)
	{
		for (int j = ipArray.jStart; j <= ipArray.jEnd; j++)
		{
			tempArray1(i, j) = ipArray(i, j).x;
			tempArray2(i, j) = ipArray(i, j).y;
		}
	}

	
	string strX = string(varName) + "X";
	const char* varNameX = strX.c_str();
	tempArray1.Variable(varNameX);

	string strY = string(varName) + "Y";
	const char* varNameY = strY.c_str();
	tempArray2.Variable(varNameY);
}


//template<class TT>
//void FieldVec2DVariable(const char * varName, const Field2D<Vector2D<TT>>& ipField)
//{
//	ArrayVec2DVariable(varName, ipField.dataArray);
////	Array2D<double> tempArray1(ipField.iStart, ipField.iRes, ipField.jStart, ipField.jRes);
////	Array2D<double> tempArray2(ipField.iStart, ipField.iRes, ipField.jStart, ipField.jRes);
////
////#pragma omp parallel for
////	for (int i = ipField.iStart; i <= ipField.iEnd; i++)
////	{
////		for (int j = ipField.jStart; j <= ipField.jEnd; j++)
////		{
////			tempArray1(i, j) = ipField(i, j).x;
////			tempArray2(i, j) = ipField(i, j).y;
////		}
////	}
////
////
////	string strX = string(varName) + "X";
////	const char* varNameX = strX.c_str();
////	tempArray1.Variable(varNameX);
////
////	string strY = string(varName) + "Y";
////	const char* varNameY = strX.c_str();
////	tempArray2.Variable(varNameY);
//}