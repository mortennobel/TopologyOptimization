//
//  UnitTest.cpp
//  TopOpt
//
//  Created by Morten Nobel-Joergensen on 5/24/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>
#include <cassert>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>

#include "UnitTest.h"
#include "Matrix.h"
#include "MatrixBand.h"

/*
void assert(bool value){
    if (value == false){
        std::cerr << "Assertion failed "<<std::endl;
        
        void *array[10];
        size_t size;
        
        // get void*'s for all entries on the stack
        size = backtrace(array, 10);
        backtrace_symbols_fd(array, size, 2);
    }
}*/

void UnitTest::run(){
    matrixTest();
    matrixSolve();
    matrixSetDataColumnOrder();
    matrixSolveFreedofs();
    matrixBandSetDataTest();
    matrixBandTest();
    matrixBandTestSolve();
    matrixBandFreedofs();    
    matrixBandTestSolveFreeDofs();
    printf("Test completed succesfully\n");
    std::cout<<std::flush;
}

void UnitTest::matrixTest(){
    //    TopOpt t(45,30,0.5,3.0,1.5);
    Matrix A (2,3);
    A.set(0, 0, 2);
    A.set(0, 1, 3);
    A.set(0, 2, 4);
    A.set(1, 0, 1);
    
    Matrix B(3,2);
    B.set(0,1,1000);
    B.set(1,0,1);
    B.set(1,1,100);
    B.set(2,1,10);
    
    A.print();
    B.print();
    printf("\n");
    Matrix C = A.multiply(B);
    
    C.print();
    printf("\n");
    Matrix D(2,2);
    D.set(0,0,2);
    D.set(1,0,9);
    D.set(0,1,95);
    D.set(1,1,985);
    tfloat vector[] = {8,2};
    printf("vTMv %f\n",D.vTransposeMultMMultV(vector));
    assertEqualFloat(5732,D.vTransposeMultMMultV(vector));
    
    Matrix AA(3,3);
    AA.set(0, 0, 2);
    AA.set(1, 0,-1);
    AA.set(2, 0, 0);
    AA.set(0, 1,-1);
    AA.set(1, 1, 2);
    AA.set(2, 1,-1);
    AA.set(0, 2, 0);
    AA.set(1, 2,-1);
    AA.set(2, 2, 1);
    
    tfloat b[] = {0,0,1};
    tfloat x[] = {0,0,0};
    bool res = AA.solve(b, x);
    assert(res);
    
    assertEqualFloat(x[0], 1);
    assertEqualFloat(x[1], 2);
    assertEqualFloat(x[2], 3);
    tfloat bExpectedValue[] = {0,0,1};
    for (int i=0;i<3;i++){
        assertEqualFloat(b[i], bExpectedValue[i]);
    }
    
    x[0] = x[1] = x[2] = 0;
}

void UnitTest::matrixSolve(){
    Matrix K(5,5);
    tfloat kVal[] = {
        0.0618,    -0.0378,     0.0017,     0.0069,    -0.0223,
        -0.0378,     0.0618,     0.0223,   -0.0309,    -0.0017, 
        0.0017,     0.0223,     0.0618,    -0.0223,     0.0069, 
        0.0069,    -0.0309,    -0.0223,     0.0618,     0.0017, 
        -0.0223,    -0.0017,     0.0069,     0.0017,     0.0618
    };
    K.setDataColumnOrder(kVal);

    printf("matrixSolve K \n");
    K.print();
    tfloat F[] = {-1,0,0,0,0};
    tfloat res[5];
    tfloat extectedRes[] = {  -44.9778,
        -36.9778,
        13.7778,
        -8.0000,
        -18.5778};
    bool solved = K.solve(F, res);
    assert(solved);
    for (int i=0;i<5;i++){
        assertEqualFloat(res[i], extectedRes[i]);
    }
}

void UnitTest::matrixSolveFreedofs(){
    using namespace std;
    Matrix K(8,8);
    tfloat KDataColumnOrder[] = {
        0.0618,     0.0223,     0.0069,     0.0017,    -0.0378,    -0.0017,    -0.0309,    -0.0223, 
        0.0223,     0.0618,    -0.0017,    -0.0378,     0.0017,     0.0069,    -0.0223,    -0.0309, 
        0.0069,    -0.0017,     0.0618,    -0.0223,    -0.0309,     0.0223,    -0.0378,     0.0017, 
        0.0017,    -0.0378,    -0.0223,     0.0618,     0.0223 ,   -0.0309,    -0.0017,     0.0069, 
        -0.0378,     0.0017,    -0.0309,     0.0223,     0.0618,    -0.0223,     0.0069,    -0.0017, 
        -0.0017,     0.0069,     0.0223,    -0.0309,    -0.0223,     0.0618,     0.0017,    -0.0378, 
        -0.0309,    -0.0223,    -0.0378,    -0.0017,     0.0069,     0.0017,     0.0618,     0.0223, 
        -0.0223,    -0.0309,     0.0017,     0.0069,    -0.0017,    -0.0378,     0.0223,     0.0618 
    };
    K.setDataColumnOrder(KDataColumnOrder);
    tfloat F[] = {0,-1,0,0,
                0,0,0,0};
    tfloat x[] = {0,0,0,0,
        0,0,0,0};

    vector<int> fixeddofs;
    fixeddofs.push_back(0);
    fixeddofs.push_back(2);
    fixeddofs.push_back(7);
    
    K.solve(F, x, fixeddofs);
    printf("__ fixed dofs K __\n");
    K.print();
    assertEqualFloat(x[1], -45.0116);
    assertEqualFloat(x[3], -37.0146);
    assertEqualFloat(x[4],  13.7833);
    assertEqualFloat(x[5], -7.9971);
    assertEqualFloat(x[6], -18.5792);
}

void UnitTest::matrixSetDataColumnOrder(){
    tfloat arrayColumnOrder[] = 
    {
        0,1,2,
        3,4,5,
        6,7,8
    };
    Matrix m(3,3);
    m.setDataColumnOrder(arrayColumnOrder);
    assertEqualFloat(m.get(2, 0), 6);
}

void UnitTest::matrixBandTest(){
    /// Column-major order storage of band-matrix 
    /// Example m x n (3 x 3 - bandwidth 2)
    ///  1 4 0
    ///  4 5 8
    ///  0 8 9 
    /// is stored as banded  
    ///  1 5 9 
    ///  4 8 0
    /// is stored in array
    ///  [1 4 5 8 9 0]
    
    MatrixBand mb(3,2);
    tfloat mbVal[] = {
        1, 4, 0,
        4, 5, 8,
        0, 8, 9 
    };
    mb.setDataColumnOrder(mbVal);
    tfloat expectedValues[] = {1, 4, 5, 8, 9, 0};
    tfloat * dataPointer = mb.getDataPointer();
    for (int i=0;i<6;i++){
        assertEqualFloat(expectedValues[i], dataPointer[i]);
    }
    printf("Band matrix\n");
    mb.print();
    
    ///  2 -1  0
    /// -1  2 -1
    ///  0 -1  1 
    MatrixBand BB(3,2);
    tfloat BBval[] = {
        2,2,1,-1,-1,0
    };
    
    BB.setBandDataColumnOrder(BBval);

    
    printf("Band matrix 2 \n");
    BB.print();
    tfloat expectedValues2[] = {2, -1, 2, -1, 1, 0};
    tfloat * dataPointer2 = BB.getDataPointer();
    for (int i=0;i<6;i++){
        printf("%f == %f\n", expectedValues2[i], dataPointer2[i]);
        assertEqualFloat(expectedValues2[i], dataPointer2[i]);
    }
    tfloat b[] = {0,0,1};
    tfloat x[] = {0,0,0};
    bool res = BB.solve(b, x);
    assert(res);
    
    assertEqualFloat(x[0], 1);
    assertEqualFloat(x[1], 2);
    assertEqualFloat(x[2], 3);
    
    printf("Solution is: %f, %f, %f \n", x[0],x[1],x[2]);

}

void UnitTest::matrixBandSetDataTest(){
    MatrixBand K1(5,5);
    
    tfloat kVal[] = {
        0.0618,     0.0618,     0.0618,     0.0618,      0.0618,             
        -0.0378,     0.0223,    -0.0223,     0.0017,    0,
        0.0017,    -0.0309,     0.0069,     0,          0,
        0.0069,    -0.0017,     0,          0,          0, 
        -0.0223,    0,          0,          0,          0
    };
    K1.setBandDataColumnOrder(kVal);
    
    MatrixBand K2(5,5);
    

     tfloat kVal2[] = {
         0.0618,    -0.0378,     0.0017,     0.0069,    -0.0223,
         -0.0378,     0.0618,     0.0223,   -0.0309,    -0.0017, 
         0.0017,     0.0223,     0.0618,    -0.0223,     0.0069, 
         0.0069,    -0.0309,    -0.0223,     0.0618,     0.0017, 
         -0.0223,    -0.0017,     0.0069,     0.0017,     0.0618
     };
     K2.setDataColumnOrder(kVal2);
    
    for (int i = 0; i < 5 * 5; i++) {
        assertEqualFloat(K1.getDataPointer()[i], K2.getDataPointer()[i]);
    }
}

void UnitTest::matrixBandTestSolve(){
    MatrixBand K(5,5);
    
    tfloat kVal[] = {
        0.0618,     0.0618,     0.0618,     0.0618,      0.0618,             
        -0.0378,     0.0223,    -0.0223,     0.0017,    0,
        0.0017,    -0.0309,     0.0069,     0,          0,
        0.0069,    -0.0017,     0,          0,          0, 
        -0.0223,    0,          0,          0,          0
    };
    
    K.setBandDataColumnOrder(kVal);
    printf("Max band\n");;
    K.print();
    tfloat F[] = {-1,0,0,0,0};
    tfloat res[5];
    tfloat extectedRes[] = {  
        -44.9778,
        -36.9778,
        13.7778,
        -8.0000,
        -18.5778};
    bool solved = K.solve(F, res);
    assert(solved);
    for (int i=0;i<5;i++){
        assertEqualFloat(res[i], extectedRes[i]);
    }
}

void UnitTest::matrixBandFreedofs(){
    tfloat vals[] = {
        0.0618,     0.0446,     0.0137,     0.0034,    -0.0755,    -0.0034,    -0.0618,    -0.0446, 
        0.0446,     0.0618,    -0.0034,    -0.0755,     0.0034,     0.0137,    -0.0446,    -0.0618, 
        0.0137,    -0.0034,     0.0618,    -0.0446,    -0.0618,     0.0446,    -0.0755,     0.0034, 
        0.0034,    -0.0755,    -0.0446,     0.0618,     0.0446,    -0.0618,    -0.0034,     0.0137, 
        -0.0755,     0.0034,    -0.0618,     0.0446,     0.0618,    -0.0446,     0.0137,    -0.0034, 
        -0.0034,     0.0137,     0.0446,    -0.0618,    -0.0446,     0.0618,     0.0034,    -0.0755, 
        -0.0618,    -0.0446,    -0.0755,    -0.0034,     0.0137,     0.0034,     0.0618,     0.0446, 
        -0.0446,    -0.0618,     0.0034,     0.0137,    -0.0034,    -0.0755,     0.0446,     0.0618         
    };
    MatrixBand mb(8,8);
    mb.setDataColumnOrder(vals);
    std::vector<int> fixeddofs;
    fixeddofs.push_back(0);
    fixeddofs.push_back(2);
    fixeddofs.push_back(7);
    mb.fixedDofs(fixeddofs);
    
    printf("MB fixed dofs:\n");
    mb.print();
    
    
    Matrix m(8,8);
    m.setDataColumnOrder(vals);
    m.fixedDofs(fixeddofs);
    printf("m fixed dofs:\n");
    m.print();
    
    for (int i=0;i<8;i++){
        for (int j=0;j<8;j++){
            assertEqualFloat(m.get(i,j), mb.get(i,j));
        }
    }
}

void UnitTest::matrixBandTestSolveFreeDofs(){
    /* Matlab
     A = [0.0618,     0.0446,     0.0137,     0.0034,    -0.0755,    -0.0034,    -0.0618,    -0.0446
     0.0446,     0.0618,    -0.0034,    -0.0755,     0.0034,     0.0137,    -0.0446,    -0.0618, 
     0.0137,    -0.0034,     0.0618,    -0.0446,    -0.0618,     0.0446,    -0.0755,     0.0034, 
     0.0034,    -0.0755,    -0.0446,     0.0618,     0.0446,    -0.0618,    -0.0034,     0.0137, 
     -0.0755,     0.0034,    -0.0618,     0.0446,     0.0618,    -0.0446,     0.0137,    -0.0034, 
     -0.0034,     0.0137,     0.0446,    -0.0618,    -0.0446,     0.0618,     0.0034,    -0.0755, 
     -0.0618,    -0.0446,    -0.0755,    -0.0034,     0.0137,     0.0034,     0.0618,     0.0446, 
     -0.0446,    -0.0618,     0.0034,     0.0137,    -0.0034,    -0.0755,     0.0446,     0.0618 ]
     B = [0 -1 0 0 0 0 0 0]'
     freedofs = [2 4 5 6 7]
     A(freedofs,freedofs)\B(freedofs,:)
     */
    
    tfloat vals[] = {
    0.0618,     0.0446,     0.0137,     0.0034,    -0.0755,    -0.0034,    -0.0618,    -0.0446, 
    0.0446,     0.0618,    -0.0034,    -0.0755,     0.0034,     0.0137,    -0.0446,    -0.0618, 
    0.0137,    -0.0034,     0.0618,    -0.0446,    -0.0618,     0.0446,    -0.0755,     0.0034, 
    0.0034,    -0.0755,    -0.0446,     0.0618,     0.0446,    -0.0618,    -0.0034,     0.0137, 
    -0.0755,     0.0034,    -0.0618,     0.0446,     0.0618,    -0.0446,     0.0137,    -0.0034, 
    -0.0034,     0.0137,     0.0446,    -0.0618,    -0.0446,     0.0618,     0.0034,    -0.0755, 
    -0.0618,    -0.0446,    -0.0755,    -0.0034,     0.0137,     0.0034,     0.0618,     0.0446, 
    -0.0446,    -0.0618,     0.0034,     0.0137,    -0.0034,    -0.0755,     0.0446,     0.0618 
    };
    MatrixBand mb(8,8);
    mb.setDataColumnOrder(vals);
    
    
    tfloat F[] = {0,-1,0,0,
        0,0,0,0};
    tfloat x[] = {0,0,0,0,
        0,0,0,0};
    
    std::vector<int> fixeddofs;
    fixeddofs.push_back(0);
    fixeddofs.push_back(2);
    fixeddofs.push_back(7);
    
    mb.solve(F, x, fixeddofs);
    assertEqualFloat(x[1], 0);
    assertEqualFloat(x[3], 16.1812);
    assertEqualFloat(x[4], 0);
    assertEqualFloat(x[5], 16.1812);
    assertEqualFloat(x[6], 0);    
}

#define IS_NAN(f) (f!=f)
void UnitTest::assertEqualFloat(tfloat f1, tfloat f2){
    float delta = abs(f1-f2);
    assert(!IS_NAN(f1));
    assert(!IS_NAN(f2));
    assert( delta < 0.001);
}