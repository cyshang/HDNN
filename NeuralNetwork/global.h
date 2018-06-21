#ifndef GLOBAL_H_
#define GLOBAL_H_

#include <cmath>
#include <ctime>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <string>
#include "Parameter.h"

//#define DEBUG_MODE
#define OUTPUT_TO_SCREEN

#ifndef EIGEN_USE_MKL_ALL
//#define EIGEN_USE_MKL_ALL
#endif // !EIGEN_USE_MKL_ALL

#ifndef PI
#define PI 3.141592653589793
#endif

using std::endl;

extern std::ofstream debug;
extern Parameter parameter;

namespace RAND {

	double Uniform(double lower = 0, double upper = 1);
	double Normal(double mean = 0, double std_dev = 1);

};
#endif // !GLOBAL_H
