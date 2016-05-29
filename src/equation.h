//
//  equation.h
//  AP
//
//  Created by yuya on 2015/08/02.
//  Copyright (c) 2015年 hy. All rights reserved.
//

#ifndef __AP__equation__
#define __AP__equation__

#include <stdio.h>
#include <vector>
#define e 1.60217657e-19
#define cl 2.99792458e8
#define SIMass(mass) (mass * e * 1.0e6 / (cl * cl))
#define SIMom(mom) (mom * e * 1.0e6 / cl)
#define NatMass(mass) (mass * cl * cl / (e * 1.0e6) )
#define NatMom(mom) (mom * cl / (e * 1.0e6) )

using namespace std;

const double Q = e;
const double mass = SIMass(938.2720813);
extern double kinetic;
extern double momentum;
extern int itheta;
extern double flux, fluy, fluz;
extern double porrection, kick;
extern vector<double> r_var;
extern vector<double> z_var;
extern vector<double> t_var;
extern vector<double> pr;
extern vector<double> pth;
extern vector<double> pz;
extern vector<double> betaFunction;

void runge_kutta(double (*)(double, double, double, double, double, double), double (*)(double, double, double, double, double, double), double (*)(double, double, double, double, double, double),
                                 double (*)(double, double, double, double, double, double), double (*)(double, double, double, double, double, double), double (*)(double, double, double, double, double, double), double);
double drdth(double, double, double, double, double, double);
double dzdth(double, double, double, double, double, double);
double dtdth(double, double, double, double, double, double);
double dprdth(double, double, double, double, double, double);
double dpthdth(double, double, double, double, double, double);
double dpzdth(double, double, double, double, double, double);
double Bx();
double By();
double Bz();
double gamma();
double lorentx(double px, double py, double pz);
double lorenty(double px, double py, double pz);
double lorentz(double px, double py, double pz);

#endif /* defined(__AP__equation__) */
