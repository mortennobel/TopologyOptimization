//
//  main.cpp
//  TopOpt
//
//  Created by Morten Nobel-Joergensen on 4/30/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <cassert>
#include <time.h>
#include "TopOpt.h"
#include "Matrix.h"
#include "MatrixBand.h"
#include "UnitTest.h"

#ifndef tfloat
#define tfloat double
#endif

int main(int argc, const char * argv[])
{
    clock_t init, final;
    
    init=clock();
    
#ifdef UNIT_TEST
    UnitTest test;
    test.run();
#else
    tfloat maxChange = 0.01;
    int nelx = 120;
    int nely = 40;
    tfloat volfrac = 0.5;
    tfloat rmin = 3.0;
#ifdef RAMP
    tfloat penal = 4.5;
    std::cout<<"Using RAMP "<<penal<<std::endl;
#else
    tfloat penal = 3.0;
    std::cout<<"Using SIMP "<<penal<<std::endl;
#endif
    TopOpt t(nelx, nely, volfrac, penal, rmin, maxChange);
#endif
    
    final=clock()-init;
    std::cout << "Run time: "<<(double)final / ((double)CLOCKS_PER_SEC);
    std::cout << std::flush;
    return 0;
}

