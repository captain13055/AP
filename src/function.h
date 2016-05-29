//
//  function.h
//  AP
//
//  Created by yuya on 2015/08/28.
//  Copyright (c) 2015年 hy. All rights reserved.
//

#ifndef __AP__function__
#define __AP__function__

#include <stdio.h>
#include "data.h"
#include "equation.h"

extern double rClosedOrbit, rClosedOrbitAt;
extern double prClosedOrbitAt;
extern double initialCondition[6];
extern double closedOrbitLength;
extern double timeMeasure;
extern vector<double> rClosedOrbitInCell;

void initialize(double*);

void existCOchekingByAccel(double *, double, int, int, int);
void findClosedOrbit(double *, double, int rAtTheta);
void getClosedOrbitLength(double *, double);
void drawClosedOrbit(double *, double);
void drawClosedOrbitDistortion(double *initialCondition, double dtheta);

void betatron(double * initial_condition, double dtheta, double dr, double dz, int totalTurn, int cellPosition);
void tuneShiftCalculation(double *, double dtheta, double dr, double dz, int iE, int eE, int totalTurn);
void tuneShiftCalculationBeta(double *, double dtheta, double dr, double dz, int iE, int eE, int totalTurn);
void Tune(double *, double dtheta, double dr, double dz, int, double &htune, double &vtune);
double Hune(double *, double dtheta, double dr, int);
double Vune(double *, double dtheta, double dz, int);

#endif /* defined(__AP__function__) */
