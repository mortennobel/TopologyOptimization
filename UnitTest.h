//
//  UnitTest.h
//  TopOpt
//
//  Created by Morten Nobel-Joergensen on 5/24/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef TopOpt_UnitTest_h
#define TopOpt_UnitTest_h

#ifndef tfloat
#define tfloat double
#endif

class UnitTest {
public:
    void run();
    void matrixTest();
    void matrixSetDataColumnOrder();
    void matrixSolve();
    void matrixSolveFreedofs();
    void matrixBandTest();
    void matrixBandSetDataTest();
    void matrixBandTestSolve();
    void matrixBandFreedofs();    
    void matrixBandTestSolveFreeDofs();
private:
    void assertEqualFloat(tfloat f1, tfloat f2);
};

#endif
