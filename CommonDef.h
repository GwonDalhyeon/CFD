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
#include <cstdlib>
#include <ctime>
#include <iomanip>

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



#define PI 3.141592
#define BC_FULL 0
#define BC_DIR -1
#define BC_OBJ -2
#define BC_NULL -3
#define BC_NEUM -4
#define BC_PER -5	
#define BC_OTHER -6
#define BC_IMPLICIT -7
#define BC_REFLECTION -8