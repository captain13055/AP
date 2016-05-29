//
//  fieldAnalysis.hpp
//  AP
//
//  Created by yuya on 2015/10/30.
//  Copyright © 2015年 hy. All rights reserved.
//

#ifndef fieldAnalysis_hpp
#define fieldAnalysis_hpp

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <math.h>
#include <vector>
#include <time.h>
#include "data.h"
#include "function.h"

using namespace std;


void BL(double dtheta, double r0, double r1);
void rPowerField(double r0, double r1); // center of one cell.
void rPowerField(double r0, double r1, const double betaD, const double betaF, const double betaGap, const double betaStraight);
void halfCellMagField(double r, double z);


#endif /* fieldAnalysis_hpp */
