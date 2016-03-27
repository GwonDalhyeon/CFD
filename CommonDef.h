#pragma once

// Headers those are used a lot
#include <iostream>
#include <fstream>
#include <assert.h>
#include <limits>
#include <Windows.h>
#include <time.h>
#include <cmath>
#include <complex>
#include <algorithm>
#include <string>
#include <omp.h>
#include <stdlib.h>  

// MATLAB Engine
#include <engine.h>
#pragma comment(lib, "Libmx.lib")
#pragma comment(lib, "libmex.lib")
#pragma comment(lib, "libeng.lib")
#pragma comment(lib, "libmat.lib")
#pragma comment(lib, "mclmcr.lib")

//#include "MACROS.h"
//#include "VECTOR_2D.h"
//#include

using namespace std;

// Typedef
//typedef double			T;
//typedef VECTOR_2D<T>	VT;
//typedef VECTOR_2D<int>	VI;
//typedef complex<T>		Tcomp;
//typedef complex<int>    Icomp;


#define PI 3.141592