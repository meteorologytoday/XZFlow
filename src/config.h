
#define X_INCLUDE_EDC
#define X_INCLUDE_GNUPLOT
#define X_INCLUDE_IO

#include "xdmlab"
#include <iostream>
#include <fstream>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace X;
using namespace std;


typedef float T;

typedef DataIndex2<T,xpts,ypts> Data;
typedef DataIndex1<T,xpts> GridX;
typedef DataIndex1<T,ypts> GridY;

GridX x;
GridY y;
Data psi, vort, u, v;

xsize step = 1000;
float dt = 3;


int main() {

	init();
	
	for(xsize i = 0; i < step; ++i) {
		
	}



	return 0;
}

void init() {
	X::Util::dLinspace1(x,(T)0.0,Lx);
	X::Util::dLinspace1(y,(T)0.0,Ly);
	psi = 0.0;
	vort = 0.0;
}
