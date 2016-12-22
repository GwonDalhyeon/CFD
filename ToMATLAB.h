#pragma once

#include "CommonDef.h"

class ToMATLAB
{
public:
	Engine* ME;
	//static double* temp1;
	//static double temp;
	//static 
	ToMATLAB();
	~ToMATLAB();


	void Command(const char * command);
	void WriteImage(const char * filename, const int & fileNum, const char * fileFormat);

	template<class TT>
	void Variable(const char * varName, const TT & values);
	template<class TT>
	void Variable(const char * varName, const int & rowNum, const int & colNum, const TT * values);
private:

};

ToMATLAB::ToMATLAB()
{
	ME = engOpen("null");
}

ToMATLAB::~ToMATLAB()
{
	engClose(ME);
}

inline void ToMATLAB::Command(const char * command)
{
	engEvalString(ME, command);
}

inline void ToMATLAB::WriteImage(const char * filename, const int & fileNum, const char * fileFormat)
{
	string str;
	
	str = string("D:\\Data\\") + (string(filename) + to_string(fileNum) + string(".") + string(fileFormat));
	str = string("saveas(gcf, '") + str + string("', '") + string(fileFormat) + string("')");
	engEvalString(ME, str.c_str());
}



template<class TT>
inline void ToMATLAB::Variable(const char * varName, const TT & values)
{
	string str = string(varName) + "=double(" + to_string(values) + ");";
	engEvalString(ME, str.c_str());
}

template<class TT>
inline void ToMATLAB::Variable(const char * varName, const int & rowNum, const int & colNum, const TT * values)
{
	if (sizeof(TT) == 4)
	{
		mxArray* dataArray = mxCreateNumericMatrix(rowNum, colNum, mxINT32_CLASS, mxREAL);
		memcpy((int*)mxGetPr(dataArray), (int*)values, sizeof(int) * rowNum*colNum);
		engPutVariable(ME, varName, dataArray);
		string str = string(varName) + "=double(transpose(" + (varName)+"));";
		// Tranaspose array.
		engEvalString(ME, str.c_str());
		mxDestroyArray(dataArray);
		dataArray = 0;
	}
	else if (sizeof(TT) == 8)
	{
		mxArray* dataArray = mxCreateDoubleMatrix(rowNum, colNum, mxREAL);
		memcpy((void*)mxGetPr(dataArray), (void*)values, sizeof(double) * rowNum*colNum);
		engPutVariable(ME, varName, dataArray);
		string str = string(varName) + "=double(transpose(" + (varName)+"));";
		// Tranaspose array.
		engEvalString(ME, str.c_str());
		mxDestroyArray(dataArray);
		dataArray = 0;
	}
	else if (sizeof(TT) == 1)
	{
		mxArray* dataArray = mxCreateLogicalMatrix(rowNum, colNum);
		memcpy((TT*)mxGetPr(dataArray), (TT*)values, sizeof(TT) * rowNum*colNum);
		engPutVariable(ME, varName, dataArray);
		string str = string(varName) + "=double(transpose(" + (varName)+"));";
		// Tranaspose array.
		engEvalString(ME, str.c_str());
		mxDestroyArray(dataArray);
		dataArray = 0;
	}
	else
	{
		cout << "Input error" << endl;
		cout << "Size of TT = " << sizeof(TT) << endl;
		assert(false);
	}
}

ToMATLAB MATLAB;



///// Reference 1

//Grid2D grid(0, 1, 3, 0, 2, 5);
//grid.Variable();
//FDX(grid);
//FDY(grid);
//Array2D<int>row(grid);
//Array2D<int>col(grid);
//VectorND<double> vec(grid.iRes);
//
//vec.Variable("vec");
//X.Variable("xy");
//Y.Variable("yy");
//row.Variable("row");
//col.Variable("col");
//MATLAB.Command("figure('units','normalized','outerposition',[0 0 1 1])");
//MATLAB.Command("surf(X,Y,xy)");
//
//engEvalString(MATLAB.ME, "openvar('SinGraph')");
//engEvalString(MATLAB.ME, "workspace");
//engEvalString(MATLAB.ME, "edit");



///// Reference 2

//if (false)
//{
//	const int arraySize = 1000;
//	const double degToRad = 0.0174;
//
//	double* SinArray = new double[arraySize];
//	double* CosArray = new double[arraySize];
//	double* Degrees = new double[arraySize];
//
//	for (int i = 0; i < arraySize; i++)
//	{
//		Degrees[i] = double(i);
//		SinArray[i] = sin(double(i)*degToRad);
//		CosArray[i] = cos(double(i)*degToRad);
//	}
//
//	mxArray* SIN = mxCreateDoubleMatrix(arraySize, 1, mxREAL);
//	memcpy((void*)mxGetPr(SIN), (void*)SinArray, sizeof(double)*arraySize);
//	engPutVariable(MATLAB.ME, "SinGraph", SIN);
//
//	double* temp = new double[3 * 4];
//
//	for (int i = 0; i < 3 * 4; i++)
//	{
//		temp[i] = i;
//	}
//
//	delete[]SinArray, CosArray, Degrees, temp;
//}
//else if (false)
//{
//	const int arraySize = 1000;
//	const double degToRad = 0.0174;
//
//	mxArray* SIN = mxCreateDoubleMatrix(arraySize, 1, mxREAL);
//
//	double* SinArray = mxGetPr(SIN);
//
//	for (int i = 0; i < arraySize; i++)
//	{
//		SinArray[i] = sin(double(i)*degToRad);
//	}
//
//	memcpy((void*)mxGetPr(SIN), (void*)SinArray, sizeof(double)*arraySize);
//	engPutVariable(MATLAB.ME, "SinGraph", SIN);
//	engEvalString(MATLAB.ME, "figure('units','normalized','outerposition',[0 0 1 1])");
//	engEvalString(MATLAB.ME, "plot(SinGraph)");
//}
//else
//{
//	const int arraySize = 1000;
//	const double degToRad = 0.0174;
//
//	mxArray* SIN = mxCreateDoubleMatrix(1, 1, mxREAL);
//
//	double* SinArray = mxGetPr(SIN);
//
//	for (int i = 0; i < arraySize; i++)
//	{
//		SinArray[i] = sin(double(i)*degToRad);
//	}
//
//	memcpy((void*)mxGetPr(SIN), (void*)SinArray, sizeof(double)*arraySize);
//	engPutVariable(MATLAB.ME, "SinGraph", SIN);
//	engEvalString(MATLAB.ME, "figure('units','normalized','outerposition',[0 0 1 1])");
//	engEvalString(MATLAB.ME, "plot(SinGraph)");
//
//	delete[] SinArray;
//}