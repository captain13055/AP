//
//  equation.cpp
//  AP
//
//  Created by yuya on 2015/08/02.
//  Copyright (c) 2015年 hy. All rights reserved.
//

#include "equation.h"
#include "data.h"


double flux;
double fluy;
double fluz;

int itheta;
double kinetic;
double momentum;
double porrection;
double kick;

vector<double> r_var;
vector<double> z_var;
vector<double> t_var;
vector<double> pr;
vector<double> pth;
vector<double> pz;

void runge_kutta(double (*dxdt)(double, double, double, double, double, double), double (*dydt)(double, double, double, double, double, double), double (*dzdt)(double, double, double, double, double, double),
                                 double (*dpxdt)(double, double, double, double, double, double), double (*dpydt)(double, double, double, double, double, double), double (*dpzdt)(double, double, double, double, double, double), double h)
{
    double x, y, z, pxw, pyw, pzw;
    double kpx[4], kpy[4], kpz[4], kx[4], ky[4], kz[4], k[6];
    double p_check;
    double t;

    x = r_var.back();
    y = t_var.back();
    z = z_var.back();
    pxw = pr.back();
    pyw = pth.back();
    pzw = pz.back();
    
    t = sqrt(pow(pxw, 2) + pow(pyw, 2) + pow(pzw, 2));
    
    kpx[0] = h * (*dpxdt)(pxw, pyw, pzw, x, y, z);
    kpy[0] = h * (*dpydt)(pxw, pyw, pzw, x, y, z);
    kpz[0] = h * (*dpzdt)(pxw, pyw, pzw, x, y, z);
    kx[0] = h * (*dxdt)(pxw, pyw, pzw, x, y, z);
    ky[0] = h * (*dydt)(pxw, pyw, pzw, x, y, z);
    kz[0] = h * (*dzdt)(pxw, pyw, pzw, x, y, z);
    
    kpx[1] = h * (*dpxdt)(pxw + kpx[0]/2, pyw + kpy[0]/2, pzw + kpz[0]/2, x + kx[0]/2, y + ky[0]/2, z + kz[0]/2);
    kpy[1] = h * (*dpydt)(pxw + kpx[0]/2, pyw + kpy[0]/2, pzw + kpz[0]/2, x + kx[0]/2, y + ky[0]/2, z + kz[0]/2);
    kpz[1] = h * (*dpzdt)(pxw + kpx[0]/2, pyw + kpy[0]/2, pzw + kpz[0]/2, x + kx[0]/2, y + ky[0]/2, z + kz[0]/2);
    kx[1] = h * (*dxdt)(pxw + kpx[0]/2, pyw + kpy[0]/2, pzw + kpz[0]/2, x + kx[0]/2, y + ky[0]/2, z + kz[0]/2);
    ky[1] = h * (*dydt)(pxw + kpx[0]/2, pyw + kpy[0]/2, pzw + kpz[0]/2, x + kx[0]/2, y + ky[0]/2, z + kz[0]/2);
    kz[1] = h * (*dzdt)(pxw + kpx[0]/2, pyw + kpy[0]/2, pzw + kpz[0]/2, x + kx[0]/2, y + ky[0]/2, z + kz[0]/2);
    
    kpx[2] = h * (*dpxdt)(pxw + kpx[1]/2, pyw + kpy[1]/2, pzw + kpz[1]/2, x + kx[1]/2, y + ky[1]/2, z + kz[1]/2);
    kpy[2] = h * (*dpydt)(pxw + kpx[1]/2, pyw + kpy[1]/2, pzw + kpz[1]/2, x + kx[1]/2, y + ky[1]/2, z + kz[1]/2);
    kpz[2] = h * (*dpzdt)(pxw + kpx[1]/2, pyw + kpy[1]/2, pzw + kpz[1]/2, x + kx[1]/2, y + ky[1]/2, z + kz[1]/2);
    kx[2] = h * (*dxdt)(pxw + kpx[1]/2, pyw + kpy[1]/2, pzw + kpz[1]/2, x + kx[1]/2, y + ky[1]/2, z + kz[1]/2);
    ky[2] = h * (*dydt)(pxw + kpx[1]/2, pyw + kpy[1]/2, pzw + kpz[1]/2, x + kx[1]/2, y + ky[1]/2, z + kz[1]/2);
    kz[2] = h * (*dzdt)(pxw + kpx[1]/2, pyw + kpy[1]/2, pzw + kpz[1]/2, x + kx[1]/2, y + ky[1]/2, z + kz[1]/2);
    
    kpx[3] = h * (*dpxdt)(pxw + kpx[2], pyw + kpy[2], pzw + kpz[2], x + kx[2], y + ky[2], z + kz[2]);
    kpy[3] = h * (*dpydt)(pxw + kpx[2], pyw + kpy[2], pzw + kpz[2], x + kx[2], y + ky[2], z + kz[2]);
    kpz[3] = h * (*dpzdt)(pxw + kpx[2], pyw + kpy[2], pzw + kpz[2], x + kx[2], y + ky[2], z + kz[2]);
    kx[3] = h * (*dxdt)(pxw + kpx[2], pyw + kpy[2], pzw + kpz[2], x + kx[2], y + ky[2], z + kz[2]);
    ky[3] = h * (*dydt)(pxw + kpx[2], pyw + kpy[2], pzw + kpz[2], x + kx[2], y + ky[2], z + kz[2]);
    kz[3] = h * (*dzdt)(pxw + kpx[2], pyw + kpy[2], pzw + kpz[2], x + kx[2], y + ky[2], z + kz[2]);
    
    k[0] = (kpx[0] + 2*kpx[1] + 2*kpx[2] + kpx[3]) / 6;
    k[1] = (kpy[0] + 2*kpy[1] + 2*kpy[2] + kpy[3]) / 6;
    k[2] = (kpz[0] + 2*kpz[1] + 2*kpz[2] + kpz[3]) / 6;
    k[3] = (kx[0] + 2*kx[1] + 2*kx[2] + kx[3]) / 6;
    k[4] = (ky[0] + 2*ky[1] + 2*ky[2] + ky[3]) / 6;
    k[5] = (kz[0] + 2*kz[1] + 2*kz[2] + kz[3]) / 6;
    
    p_check = sqrt(pow(pxw + k[0], 2) + pow(pyw + k[1], 2) + pow(pzw + k[2], 2));
    
    pr.push_back(pxw + k[0]);
    pth.push_back(pyw + k[1]);
    pz.push_back(pzw + k[2]);
    r_var.push_back(x + k[3]);
    t_var.push_back(y + k[4]);
    z_var.push_back(z + k[5]);
    
    pr.back() = pr.back() * porrection / p_check;
    pth.back() = pth.back() * porrection / p_check;
    pz.back() = pz.back() * porrection / p_check;
    
    /*
    if (useToscaField) {
        if (r_var.back() > r_max || r_var.back() < r_min || z_var.back() > z_max || z_var.back() < z_min) {
            cout << "r: " << r_var.back() << endl
            << "theta: " << (double)itheta / dth_decimal << endl
            << "z: " << z_var.back() << endl;
            cout << "Divergence KUTTA" << endl;
            cout << "Magnet Field" << endl;
            cout << "r: " << r_min << " ~ " << r_max << endl;
            cout << "z: " << z_min << " ~ " << z_max << endl;
            exit(0);;
        }
	}*/
}


//differential equations
double Bx()
{
    double bx;
    if (useToscaField == true) {
        bx = flux;
    }
    else {
        bx = 0.0;
    }
    
    return bx;
}

double By()
{
    double by;
    if (useToscaField == true) {
        by = fluy;
    }
    else {
        by = 0.0;
    }
    
    return by;
}

double Bz()
{
    double bz, r0, b0;
    if (useToscaField == true) {
        bz = fluz;
    }
    else {
        r0 = 0.456557;
        b0 = -1.0;
        bz = b0 * pow(r_var.back() / r0, -0.5);
    }
    return bz;
}


double lorentx(double px, double py, double pz)
{
    return (py * Bz() - pz * By() );
}
double lorenty(double px, double py, double pz)
{
    return (pz * Bx() - px * Bz() );
}
double lorentz(double px, double py, double pz)
{
    return (px * By() - py * Bx() );
}
double gamma() {
    double g;
    
    g = (kinetic + NatMass(mass)) / NatMass(mass);
    return g;
}


double drdth(double pr, double pth, double pz, double r, double z, double t)
{
    return r * pr / pth;
}
double dzdth(double pr, double pth, double pz, double r, double z, double t)
{
    return r * pz / pth;
}
double dtdth(double pr, double pth, double pz, double r, double z, double t)
{
    return r * gamma() * mass / pth;
}
double dprdth(double pr, double pth, double pz, double r, double z, double t)
{
    return Q * r * lorentx(pr, pth, pz) / pth + pth;
}
double dpthdth(double pr, double pth, double pz, double r, double z, double t)
{
    return Q * r * lorenty(pr, pth, pz) / pth - pr;
}
double dpzdth(double pr, double pth, double pz, double r, double z, double t)
{
    return Q * r * lorentz(pr, pth, pz) / pth;
}




