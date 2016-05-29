//
//  function.cpp
//  AP
//
//  Created by yuya on 2015/08/28.
//  Copyright (c) 2015年 hy. All rights reserved.
//

#include "function.h"

double rClosedOrbit, rClosedOrbitAt;
double prClosedOrbitAt;
double initialCondition[6];
double closedOrbitLength;
double timeMeasure;
vector<double> rClosedOrbitInCell;


void initialize(double *init) {
    r_var.clear();
    z_var.clear();
    t_var.clear();
    pr.clear();
    pth.clear();
    pz.clear();
    r_var.push_back(init[0]);
    t_var.push_back(init[1]);
    z_var.push_back(init[2]);
    pr.push_back(init[3]);
    pth.push_back(init[4]);
    pz.push_back(init[5]);
    porrection = sqrt(pow(initialCondition[3], 2) + pow(initialCondition[4], 2) + pow(initialCondition[5], 2));
}

void existCOchekingByAccel(double *initialCondition, double dth, int iE, int eE, int dE) {
    ofstream ofs("existCO.csv");
    for (int ikinetic = iE; ikinetic <= eE; ikinetic += dE) {
        cout << "CO checking in " << ikinetic << endl;
        kinetic = ikinetic;
        initialCondition[4] = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic)); //pth
        porrection = sqrt(pow(initialCondition[3], 2) + pow(initialCondition[4], 2) + pow(initialCondition[5], 2));
        findClosedOrbit(initialCondition, dth, 0);
        ofs << kinetic << "," << rClosedOrbit << endl;
    }
}

void findClosedOrbit(double *initialCondition, double dth, int rAtTheta) {
    int repeatTime = 0;
    double tmpClosed;
    
    tmpClosed = 0.0;
    initialize(initialCondition);
    if (useToscaField) {
        if (r_var.at(0) > r_max || r_var.at(0) < r_min || z_var.at(0) > z_max || z_var.at(0) < z_min) {
            cout << "r: " << r_var.back() << endl
            << "theta: " << (double)itheta / dth_decimal << endl
            << "z: " << z_var.back() << endl;
            cout << "Divergence in findClosedOrbit" << endl;
            cout << "Magnet Field" << endl;
            cout << "r: " << r_min << " ~ " << r_max << endl;
            cout << "z: " << z_min << " ~ " << z_max << endl;
            rClosedOrbit = 0;
            return;
        }
    }
repeat:
    rClosedOrbitInCell.clear();
    for (int ith = 0 * dth_decimal; ith <= 360 / symmetryNumber * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        if (useToscaField) {
            getAveFields();
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dth * M_PI / 180);
        if (itheta == rAtTheta * dth_decimal && rAtTheta != 0) {
            rClosedOrbitAt = r_var.back();
            prClosedOrbitAt = pr.back();
        }
        rClosedOrbitInCell.push_back(r_var.back());
        //Divergence action
        if (useToscaField) {
            if (r_var.back() > r_max || r_var.back() < r_min || z_var.back() > z_max || z_var.back() < z_min) {
                cout << "r: " << r_var.back() << endl
                << "theta: " << (double)itheta / dth_decimal << endl
                << "z: " << z_var.back() << endl;
                cout << "Divergence in findClosedOrbit" << endl;
                cout << "Magnet Field" << endl;
                cout << "r: " << r_min << " ~ " << r_max << endl;
                cout << "z: " << z_min << " ~ " << z_max << endl;
                rClosedOrbit = 0;
                exit(0);
            }
        }
    }
    
    if (useToscaField) {
        if (fabs(r_var.at(0) - r_var.back()) > 1e-10 || NatMom(fabs(pth.at(0) - pth.back())) > 1e-6) {
            if (repeatTime > 100) {
                cout << "Error too large." << endl;
                cout << "pr at    0 deg: " << setprecision(10) << NatMom(pth.at(0)) << endl;
                cout << "pr at cell deg: " << setprecision(10) << NatMom(pth.back()) << endl;
                rClosedOrbit = 0;
                exit(0);;
            }
            tmpClosed = (r_var.at(0) + r_var.back()) / 2;
            initialize(initialCondition);
            r_var.at(0) = tmpClosed;
            repeatTime++;
            goto repeat;
        }
        else {
            tmpClosed = r_var.back();
            cout << "closed orbit found." << endl
            << "initial r: " << fixed << setprecision(10) << tmpClosed << endl;
        
            rClosedOrbit = tmpClosed;
        }
    }
    else if (!useToscaField) {
        // Weak focusing
        rClosedOrbit = 0.456557;
    }

}

void getClosedOrbitLength(double *initialCondition, double dth) {
    double length = 0;
    double dtheta = M_PI / (180 * dth_decimal);
    closedOrbitLength = 0;
    initialCondition[4] = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic)); //pth
    initialize(initialCondition);
    r_var.at(0) = rClosedOrbit;
    for (int ith = 0 * dth_decimal; ith <= 360 * dth_decimal; ith += 1) {
        itheta = ith;
        if (useToscaField) {
            getAveFields();
        }
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dth * M_PI / 180);
        timeMeasure += t_var.back();
        length += r_var.back() * dtheta;
    }
    
    closedOrbitLength = length;
}


void drawClosedOrbit(double *initialCondition, double dtheta) {
    ofstream closedOrbitXY("closedOrbitXY.csv");
    ofstream excursion("excursion.csv");
    initialize(initialCondition);
    r_var.at(0) = rClosedOrbit;
    z_var.at(0) = 0;

    for (int ith = 0; ith <= 360 * dth_decimal; ith += 1) {
        itheta = ith;
        if (useToscaField) {
            getAveFields();
        }
        if (ith % 10 == 0) {
            closedOrbitXY << r_var.back() * cos((double)ith / dth_decimal * M_PI / 180) * 100 << ","
            << r_var.back() * sin((double)ith / dth_decimal * M_PI / 180) * 100 << endl;
            
            if (ith <= 360 / symmetryNumber * dth_decimal) {
                excursion << (double)ith / dth_decimal << "," << r_var.back() * 100 << endl;
            }
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
    cout << "excursion " << 100 * (*max_element(r_var.begin(), r_var.end()) - *min_element(r_var.begin(), r_var.end())) << endl;
}

void drawClosedOrbitDistortion(double *initialCondition, double dtheta) {
    ofstream COD("COD.csv");
    ofstream orbit("COD_ORBIT.csv");
    double tmpClosed;
    bool cod = false;
    
    initialize(initialCondition);
    findClosedOrbit(initialCondition, dtheta, 0);
    initialize(initialCondition);
    r_var.at(0) = rClosedOrbit;
    z_var.at(0) = 0;
    cout << "Injection" << endl;
    cout << r_var.back() << endl << z_var.back() << endl << NatMom(pr.back()) << endl << NatMom(pth.back()) << endl << NatMom(pz.back()) << endl;
    if (useToscaField) {
        if (r_var.at(0) > r_max || r_var.at(0) < r_min || z_var.at(0) > z_max || z_var.at(0) < z_min) {
            cout << "r: " << r_var.back() << endl
            << "theta: " << (double)itheta / dth_decimal << endl
            << "z: " << z_var.back() << endl;
            cout << "Divergence in findClosedOrbit" << endl;
            cout << "Magnet Field" << endl;
            cout << "r: " << r_min << " ~ " << r_max << endl;
            cout << "z: " << z_min << " ~ " << z_max << endl;
            rClosedOrbit = 0;
            return;
        }
    }
repeat:
    for (int ith = 0 * dth_decimal; ith <= 360 * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        if (useToscaField) {
            getAveFields();
        }
        flux = 0;
        fluy = 0;
        if (ith == 180 * dth_decimal) {
            momentum  = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic));
            pr.back() += momentum * kick;
            pth.back() = sqrt(pow(momentum, 2) - pow(pr.back(), 2));
        }
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
        if (cod) {
            if (ith % 10 == 0) {
                COD << ith / (double)dth_decimal << "," << r_var.back() - rClosedOrbitInCell.at(ith % angleInCell) << endl;
                orbit << r_var.back() * cos((double)ith / dth_decimal * M_PI / 180) * 100 << "," << r_var.back() * sin((double)ith / dth_decimal * M_PI / 180) * 100 << endl;
            }
        }
        //Divergence action
        if (useToscaField) {
            if (r_var.back() > r_max || r_var.back() < r_min || z_var.back() > z_max || z_var.back() < z_min) {
                cout << "r: " << r_var.back() << endl
                << "theta: " << (double)itheta / dth_decimal << endl
                << "z: " << z_var.back() << endl;
                cout << "Divergence in findClosedOrbit" << endl;
                cout << "Magnet Field" << endl;
                cout << "r: " << r_min << " ~ " << r_max << endl;
                cout << "z: " << z_min << " ~ " << z_max << endl;
                rClosedOrbit = 0;
                exit(0);
            }
        }
    }
    
    if (useToscaField) {
        if (fabs(r_var.at(0) - r_var.back()) > 1e-10) {
            tmpClosed = (r_var.at(0) + r_var.back()) / 2;
            initialize(initialCondition);
            r_var.at(0) = tmpClosed;
            z_var.at(0) = 0;
            cout << tmpClosed << endl;
            goto repeat;
        }
        else {
            tmpClosed = r_var.back();
            initialize(initialCondition);
            r_var.at(0) = tmpClosed;
            z_var.at(0) = 0;
            cout << "closed orbit found." << endl
            << "initial r: " << fixed << setprecision(10) << tmpClosed << endl;
            if (!cod) {
                cod = true;
                rClosedOrbit = tmpClosed;
                cout << "closed orbit distorsion found." << endl;
                goto repeat;
            }
        }
    }
    cout << "360" << endl;
    cout << r_var.back() << endl << z_var.back() << endl << NatMom(pr.back()) << endl << NatMom(pth.back()) << endl << NatMom(pz.back()) << endl;
}

void betatron(double * initial_condition, double dtheta, double dr, double dz, int totalTurn, int cellPosition) {
    int tmp_theta, thPeriod;
    stringstream fileName;
    ofstream lattice;
    initialize(initial_condition);
    r_var.at(0) = rClosedOrbit + dr;
    z_var.at(0) = dz;
    thPeriod = cellPosition * dth_decimal;
    fileName << cellPosition << "_" << "lattice.csv";
    lattice.open(fileName.str());
    
    lattice << "theta [deg]" << "," << "r [m]" << "," << "z [m]" << ","
    << "pr [MeV/c]" << "," << "pth [Mev/c]" << "," << "pz [MeV/c]" << "," << "pr / pth" << endl;
    for (int ith = 0; ith <= 360 * totalTurn * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        tmp_theta = ith;

        if (useToscaField) {
            getAveFields();
        }
        if (ith % (360 * 20 * dth_decimal) == 0 && ith != 0) {
            cout << ith / (360 * dth_decimal) << endl;
        }
        
        if (tmp_theta == thPeriod) {
            lattice << tmp_theta / dth_decimal << ","
            << fixed << setprecision(20)
            << r_var.back() - rClosedOrbitAt << ","
            << z_var.back() << ","
            << NatMom(pr.back()) - NatMom(prClosedOrbitAt) << ","
            << NatMom(pth.back()) << ","
            << NatMom(pz.back()) << ","
            << NatMom(pr.back()) / NatMom(pth.back()) << endl;
            
            thPeriod += 360 / symmetryNumber * dth_decimal;
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
}

void tuneShiftCalculation(double *initialCondition, double dtheta, double dr, double dz, int iE, int eE, int totalTurn) {
    double tuneHorizontal, tuneVertical;
    ofstream tuneShift("tuneShift.csv");
    cout << "****   Tune Shift Calculation starts   ****" << endl;
    tuneShift << "E" << "," << "Qh" << "," << "Qv" << endl;
    for (int ikinetic = iE; ikinetic <= eE; ikinetic ++) {
        kinetic = ikinetic;
        initialCondition[4] = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic)); //pth
        porrection = sqrt(pow(initialCondition[3], 2) + pow(initialCondition[4], 2) + pow(initialCondition[5], 2));
        
        if (fabs(kick) < 1e-6) {
            findClosedOrbit(initialCondition, dtheta, 0);
        }
        
        tuneHorizontal = Hune(initialCondition, dtheta, dr, totalTurn);
        tuneVertical = Vune(initialCondition, dtheta, dz, totalTurn);
        tuneShift << ikinetic << "," << tuneHorizontal << "," << tuneVertical << endl;
    }
}


void tuneShiftCalculationBeta(double *initialCondition, double dtheta, double dr, double dz, int iE, int eE, int totalTurn) {
    double tuneHorizontal, tuneVertical;
    ofstream tuneShift("tuneShiftBeta.csv");
    
    tuneShift << "E" << "," << "Qh" << "," << "Qv" << "r" << endl;
    for (int ikinetic = iE; ikinetic <= eE; ikinetic++) {
        kinetic = ikinetic;
        initialCondition[4] = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic)); //pth
        porrection = sqrt(pow(initialCondition[3], 2) + pow(initialCondition[4], 2) + pow(initialCondition[5], 2));
        findClosedOrbit(initialCondition, dtheta, 0);
        Tune(initialCondition, dtheta, dr, dz, totalTurn, tuneHorizontal, tuneVertical);
		cout << endl;
		tuneShift << ikinetic << "," << tuneHorizontal << "," << tuneVertical << "," << rClosedOrbit << endl;
    }
}

void Tune(double * initial_condition, double dtheta, double dr, double dz, int totalTurn, double &htune, double &vtune) {
    int tmp_theta, thPeriod;
    double r1, r2, z1, z2, p1, p2, l1, l2, l3, x;
    double summuh, summuv;
    initialize(initial_condition);
    r_var.at(0) = rClosedOrbit + dr;
    z_var.at(0) = dz;
    summuh = 0;
    summuv = 0;
    htune = 0;
    vtune = 0;
    thPeriod = 360 * dth_decimal / symmetryNumber ;
    cout << "calculating Tune at " << (int)(kinetic + 1e-4) << " MeV ...." << endl;
    for (int ith = 0; ith <= 360 * totalTurn * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        tmp_theta = ith;
        
        if (useToscaField) {
            getAveFields();
        }
        
        if (ith % (360 * 20 * dth_decimal) == 0 && ith != 0) {
	  cout << ith / (360 * dth_decimal) << endl;
        }
        
        if ((tmp_theta % thPeriod) == 0 && (tmp_theta != 0)) {
	  momentum  = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic));
	  r1 = (rClosedOrbit - r_var.at(ith)) * 1000;
	  r2 = (rClosedOrbit - r_var.at(ith - thPeriod)) * 1000;
	  p1 = NatMom(pr.at(ith)) / NatMom(momentum) * 1000;
	  p2 = NatMom(pr.at(ith - thPeriod)) / NatMom(momentum) * 1000;
	  l1 = sqrt(pow(r1, 2) + pow(p1, 2));
	  l2 = sqrt(pow(r2, 2) + pow(p2, 2));
	  l3 = sqrt(pow(r1 - r2, 2) + pow(p1 - p2, 2));
	  x = ((l1*l1 + l2*l2 - l3*l3) / (2 * l1 * l2));
	  summuh += acos(x);
	  //    from here why different Tune and Hune's x
	  z1 = z_var.at(ith) * 1000;
	  z2 = z_var.at(ith - thPeriod) * 1000;
	  p1 = NatMom(pz.at(ith)) / NatMom(momentum) * 1000;
	  p2 = NatMom(pz.at(ith - thPeriod)) / NatMom(momentum) * 1000;
	  l1 = sqrt(pow(z1, 2) + pow(p1, 2));
	  l2 = sqrt(pow(z2, 2) + pow(p2, 2));
	  l3 = sqrt(pow(z1 - z2, 2) + pow(p1 - p2, 2));
	  x = ((l1*l1 + l2*l2 - l3*l3) / (2 * l1 * l2));
	  summuv += acos(x);
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
    
    htune = summuh / (totalTurn * 2 * M_PI);
    vtune = summuv / (totalTurn * 2 * M_PI);
}

double Hune(double * initial_condition, double dtheta, double dr, int totalTurn) {
    int tmp_theta, thPeriod;
    double r1, r2, p1, p2, l1, l2, l3, x;
    double summu, HTune;
    initialize(initial_condition);
    r_var.at(0) = rClosedOrbit + dr;
    cout << "Hune starts at r:" << endl;
    cout << rClosedOrbit << endl;
    summu = 0;
    if (fabs(kick) < 1e-6) {
        thPeriod = 360 / symmetryNumber * dth_decimal;
    }
    else
        thPeriod = 360 * dth_decimal;
    
    cout << "calculating Hune at " << (int)(kinetic + 1e-4) << " MeV ...." << endl;
    for (int ith = 0; ith <= 360 * totalTurn * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        tmp_theta = ith;
        
        if (useToscaField) getAveFields();
        
        if (ith % (360 * 20 * dth_decimal) == 0 && ith != 0) {
            cout << ith / (360 * dth_decimal) << endl;
        }
        if (ith % (180 * dth_decimal) == 0 && ith != 0 && ith % (360 * dth_decimal) != 0) {
            momentum  = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic));
            pr.back() += momentum * kick;
            pth.back() = sqrt(pow(momentum, 2) - pow(pr.back(), 2));
            pz.back() = 0;
        }
        
        if ((tmp_theta % thPeriod) == 0 && (tmp_theta != 0)) {
            r1 = (rClosedOrbit - r_var.at(ith)) * 1000;
            r2 = (rClosedOrbit - r_var.at(ith - thPeriod)) * 1000;
            p1 = NatMom(pr.at(ith)) / NatMom(momentum) * 1000;
            p2 = NatMom(pr.at(ith - thPeriod)) / NatMom(momentum) * 1000;
            l1 = sqrt(pow(r1, 2) + pow(p1, 2));
            l2 = sqrt(pow(r2, 2) + pow(p2, 2));
            l3 = sqrt(pow(r1 - r2, 2) + pow(p1 - p2, 2));
            x = ((l1*l1 + l2*l2 - l3*l3) / (2 * l1 * l2));
           
            if (fabs(x - 1) < 1e-12) {
                summu += acos(1);
            }
            else if (fabs(x + 1) < 1e-12) {
                summu += acos(-1);
            }
            else {
                summu += acos(x);
            }
            
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
    HTune = summu / (totalTurn * 2 * M_PI);
    
    return HTune;
}

double Vune(double * initial_condition, double dtheta, double dz, int totalTurn) {
    int tmp_theta, thPeriod;
    double z1, z2, p1, p2, l1, l2, l3, x;
    double summu, VTune;
    initialize(initial_condition);
    r_var.at(0) = rClosedOrbit;
    z_var.at(0) = dz;
    summu = 0;
    
    if (fabs(kick) < 1e-6) {
        thPeriod = 360 / symmetryNumber * dth_decimal;
    }
    else
        thPeriod = 360 * dth_decimal;
    
    cout << "calculating Vune at " << (int)(kinetic + 1e-4) << " MeV ...." << endl;
    for (int ith = 0; ith <= 360 * totalTurn * dth_decimal; ith += 1) {
        //getAveFieldsで現在の粒子の位置から、実効磁場を求める。
        //dtheta毎に、getAveFieldsをコールする必要がある。
        itheta = ith;
        tmp_theta = ith;
        
        if (useToscaField) {
            getAveFields();
        }
        if (ith % (360 * 20 * dth_decimal) == 0 && ith != 0) {
            cout << ith / (360 * dth_decimal) << endl;
        }
        /* COD kick not complete. 20160215
        if (ith % (180 * dth_decimal) == 0 && ith != 0 && ith % (360 * dth_decimal) != 0) {
            momentum  = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic));
            pr.back() += momentum * kick;
            pth.back() = sqrt(pow(momentum, 2) - pow(pr.back(), 2));
        }*/
        if ((tmp_theta % thPeriod) == 0 && (tmp_theta != 0)) {
            z1 = z_var.at(ith) * 1000;
            z2 = z_var.at(ith - thPeriod) * 1000;
            p1 = NatMom(pz.at(ith)) / NatMom(momentum) * 1000;
            p2 = NatMom(pz.at(ith - thPeriod)) / NatMom(momentum) * 1000;
            l1 = sqrt(pow(z1, 2) + pow(p1, 2));
            l2 = sqrt(pow(z2, 2) + pow(p2, 2));
            l3 = sqrt(pow(z1 - z2, 2) + pow(p1 - p2, 2));
            x = ((l1*l1 + l2*l2 - l3*l3) / (2 * l1 * l2));
            
            if (fabs(x - 1) < 1e-12) {
                summu += acos(1);
            }
            else if (fabs(x + 1) < 1e-12) {
                summu += acos(-1);
            }
            else {
                summu += acos(x);
            }
            
        }
        
        runge_kutta(drdth, dtdth, dzdth, dprdth, dpthdth, dpzdth, dtheta * M_PI / 180);
    }
    VTune = summu / (totalTurn * 2 * M_PI);
    
    return VTune;
}
