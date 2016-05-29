//
//  fieldAnalysis.cpp
//  AP
//
//  Created by yuya on 2015/10/30.
//  Copyright © 2015年 hy. All rights reserved.
//

#include "fieldAnalysis.hpp"

void BL(double dtheta, double r0, double r1) {
    ofstream BLF("BLF.csv");
    ofstream BLD("BLD.csv");
	ofstream Grad("Grad.csv");
    int bF, bD, bG, bS, bCell;
    double sum;
    double drAve, dBLAve;
    const int aveNum = 1;
    vector<double> local_k, local_kR, vec_BL_F, vec_BL_D;
    
    bCell = 360 * dth_decimal / symmetryNumber ;
    bF = (int)(angleFocus * dth_decimal + 1e-7);
    bD = (int)(angleDefocus * dth_decimal + 1e-7);
    bG = (int)(angleGap * dth_decimal + 1e-7);
    bS = (int)(angleStraight * dth_decimal + 1e-7);
    if (magnetTriplet == DFD) {
        bD = (int)(angleDefocus * 2 * dth_decimal + 1e-7);
        bF = (int)(angleFocus * dth_decimal / 2 + 1e-7);
        bG = (int)(angleGap * dth_decimal + 1e-7);
        bS = (int)(angleStraight * dth_decimal + 1e-7);
    }
    
    if (bCell != bD + bG + bF + bS + bF + bG + bD && magnetTriplet == FDF) {
        cout << "symmetry not matched" << endl;
        cout << "Cell angle: " << bCell << endl;
        cout << "Total input angle: " << bD + bG + bF + bS + bF + bG + bD << endl;
        exit(0);
    }
    
    else if (bCell != bD + bG + bF + bS + bF + bG + bD && magnetTriplet == DFD) {
        cout << "symmetry not matched" << endl;
        cout << "Cell angle: " << bCell << endl;
        cout << "Total input angle: " << bD + bG + bF + bS + bF + bG + bD << endl;
        exit(0);;
    }
    
    r_var.clear();
    z_var.clear();
    if (magnetTriplet == DFD) {
        // F section
        for (double r = r0; r <= r1; r += 0.01) {
            cout << "F r:" << r << endl;
            sum = 0;
            r_var.push_back(r);
            z_var.push_back(0);
            for (int ith = -bF - bG / 2, i = 0; ith <= bF + bG / 2; ith++, i++) {
                itheta = ith;
                getAveFields();
                sum += fluz * r_var.back() * rad(dtheta);
            }
            vec_BL_F.push_back(sum);
        }
        for (int i = 1; i <= r_var.size() - 2; i++) {
            local_kR.push_back(r_var.at(i));
            local_k.push_back(r_var.at(i) / vec_BL_F.at(i) * (vec_BL_F.at(i + 1) - vec_BL_F.at(i - 1)) /
							  (r_var.at(i + 1) - r_var.at(i - 1)));
        }
        
        for (int i = 0; i <= local_k.size() - 1; i++) {
            BLF << local_kR.at(i) << "," << local_k.at(i) << endl;
        }
        r_var.clear();
        z_var.clear();
        local_k.clear();
		local_kR.clear();
        // D section
        for (double r = r0; r <= r1; r += 0.01) {
            cout << "D r:" << r << endl;
            sum = 0;
            r_var.push_back(r);
            z_var.push_back(0);
            for (int ith = bF + bG / 2; ith <= bF + bG + bD + bS / 2; ith++) {
                itheta = ith;
                getAveFields();
                sum += fluz * r_var.back() * rad(dtheta);
            }
            vec_BL_D.push_back(sum);
        }
        for (int i = 1; i <= r_var.size() - 2; i++) {
            local_kR.push_back(r_var.at(i));
            local_k.push_back(r_var.at(i) / vec_BL_D.at(i) * (vec_BL_D.at(i + 1) - vec_BL_D.at(i - 1)) /
							  (r_var.at(i + 1) - r_var.at(i - 1)));
			Grad << r_var.at(i) << "," <<
				(vec_BL_D.at(i + 1) - vec_BL_D.at(i - 1)) / (r_var.at(i + 1) - r_var.at(i - 1)) << endl;
        }
        
        for (int i = 0; i <= local_k.size() - 1; i++) {
            BLD << local_kR.at(i) << "," << local_k.at(i) << endl;
        }
        r_var.clear();
        z_var.clear();
        vec_BL_D.clear();
    }

    else if (magnetTriplet == FDF) {
        // F section
        for (double r = r0; r <= r1; r += 0.005) {
            cout << "F r:" << r << endl;
            sum = 0;
            r_var.push_back(r);
            z_var.push_back(0);
            for (int ith = 0; ith <= bS / 2 + bF + bG / 2; ith++) {
                itheta = ith;
                getAveFields();
                sum += fluz * r_var.back() * rad(dtheta);
            }
            vec_BL_F.push_back(sum);
        }
        for (int i = aveNum; i <= r_var.size() - 2 - aveNum; i++) {
            drAve = 0;
            dBLAve = 0;
            /*
            for (int j = i - aveNum; j <= i + aveNum; j++) {
                drAve += pow(r_var.at(j + 1) - r_var.at(j), 2);
                dBLAve += pow(vec_BL_F.at(j + 1) - vec_BL_F.at(j), 2);
            }
            drAve = sqrt(drAve);
            dBLAve = sqrt(dBLAve);
            drAve = drAve / (2 * aveNum);
            dBLAve = dBLAve / (2 * aveNum);*/
            for (int j = i - aveNum; j <= i + aveNum; j++) {
                drAve += r_var.at(j + 1) - r_var.at(j);
                dBLAve += vec_BL_F.at(j + 1) - vec_BL_F.at(j);
            }
            drAve = drAve / (2 * aveNum);
            dBLAve = dBLAve / (2 * aveNum);
            
            local_k.push_back(fabs(r_var.at(i) / vec_BL_F.at(i) * dBLAve / drAve));
        }
        for (int i = 0; i <= local_k.size() - 1; i++) {
            BLF << r_var.at(i) << "," << local_k.at(i) << endl;
        }
        r_var.clear();
        z_var.clear();
        local_k.clear();
        // D section
        for (double r = r0; r <= r1; r += 0.005) {
            cout << "D r:" << r << endl;
            sum = 0;
            r_var.push_back(r);
            z_var.push_back(0);
            for (int ith = bS / 2 + bF + bG / 2; ith <= bS / 2 + bF + bG + bD + bD + bG / 2; ith++) {
                itheta = ith;
                getAveFields();
                sum += fluz * r_var.back() * rad(dtheta);
            }
            vec_BL_D.push_back(sum);
        }
        
        for (int i = aveNum; i <= r_var.size() - 2 - aveNum; i++) {
            drAve = 0;
            dBLAve = 0;
            /*
            for (int j = i - aveNum; j <= i + aveNum; j++) {
                drAve += pow(r_var.at(j + 1) - r_var.at(j), 2);
                dBLAve += pow(vec_BL_D.at(j + 1) - vec_BL_D.at(j), 2);
            }
            drAve = sqrt(drAve);
            dBLAve = sqrt(dBLAve);
            drAve = drAve / (2 * aveNum);
            dBLAve = dBLAve / (2 * aveNum);
            */
            for (int j = i - aveNum; j <= i + aveNum; j++) {
                drAve += r_var.at(j + 1) - r_var.at(j);
                dBLAve += vec_BL_D.at(j + 1) - vec_BL_D.at(j);
            }
            drAve = drAve / (2 * aveNum);
            dBLAve = dBLAve / (2 * aveNum);
            
            local_k.push_back(fabs(r_var.at(i) / vec_BL_D.at(i) * dBLAve / drAve));
        }
        
        for (int i = 0; i <= local_k.size() - 1; i++) {
            BLD << r_var.at(i) << "," << local_k.at(i) << endl;
        }
        r_var.clear();
        z_var.clear();
        vec_BL_D.clear();
    }
    
}

void rPowerField(double r0, double r1) {
    double drAve, dBAve;
    int aveNum = 10;
    vector<double> vec_Bz, local_k;
    ofstream centMag("rocalCenterMagnet.csv");
    ofstream LK("localk.csv");
    r_var.clear();
    z_var.clear();
    
    for (double r = r0; r < r1; r += 0.001) {
        itheta = 360 * dth_decimal / symmetryNumber / 2;
        r_var.push_back(r);
        z_var.push_back(0);
        getAveFields();
        vec_Bz.push_back(fluz);
        centMag << r_var.back() << "," << z_var.back() << "," << (double)itheta / dth_decimal << "," << flux << "," << fluy << "," << fabs(fluz) << endl;
    }

    for (int i = aveNum; i <= r_var.size() - 2 - aveNum; i++) {
        drAve = 0;
        dBAve = 0;
        for (int j = i - aveNum; j <= i + aveNum; j++) {
            drAve += pow(r_var.at(j + 1) - r_var.at(j), 2);
            dBAve += pow(vec_Bz.at(j + 1) - vec_Bz.at(j), 2);
        }
        drAve = sqrt(drAve);
        dBAve = sqrt(dBAve);
        drAve = drAve / (2 * aveNum);
        dBAve = dBAve / (2 * aveNum);
        
        local_k.push_back(fabs(r_var.at(i) / vec_Bz.at(i) * dBAve / drAve));
    }
    
    for (int i = 0; i <= local_k.size() - 1; i++) {
        LK << r_var.at(i) << "," << local_k.at(i) << endl;
    }

    r_var.clear();
    z_var.clear();
}
// decide r and r to the power k.
// theta is center of each magnet
void rPowerField(double r0, double r1, double betaD, double betaF, double betaGap, double betaStraight) {
    ofstream dMag("rocalDield.csv");
    ofstream fMag("rocalField.csv");
    ofstream LKF("focalk.csv");
    ofstream LKD("docalk.csv");
    vector<double> vec_F_Bz, vec_D_Bz;
    vector<double> local_k;
    double fCenterBeta, dCenterBeta;
    double drAve, dBAve;
    int aveNum = 10;
    
    r_var.clear();
    z_var.clear();
    if (magnetTriplet == FDF) {
        fCenterBeta = betaStraight + betaF / 2;
        dCenterBeta = fCenterBeta + betaF / 2 + betaGap + betaD;
        for (double r = r0; r < r1; r += 0.001) {
            itheta = fCenterBeta * dth_decimal;
            r_var.push_back(r);
            z_var.push_back(0);
            getAveFields();
            vec_F_Bz.push_back(fabs(fluz));
            fMag << r_var.back() << "," << z_var.back() << "," << (double)itheta / dth_decimal << "," << flux << "," << fluy << "," << fabs(fluz) << endl;
        }
        
        for (int i = aveNum; i <= r_var.size() - 2 - aveNum; i++) {
            drAve = 0;
            dBAve = 0;
            for (int j = i - aveNum; j <= i + aveNum; j++) {
                drAve += pow(r_var.at(j + 1) - r_var.at(j), 2);
                dBAve += pow(vec_F_Bz.at(j + 1) - vec_F_Bz.at(j), 2);
            }
            drAve = sqrt(drAve);
            dBAve = sqrt(dBAve);
            drAve = drAve / (2 * aveNum);
            dBAve = dBAve / (2 * aveNum);
            
            local_k.push_back(fabs(r_var.at(i) / vec_F_Bz.at(i) * dBAve / drAve));
        }
        
        for (int i = 0; i <= local_k.size() - 1; i++) {
            LKF << r_var.at(i) << "," << local_k.at(i) << endl;
        }
        
        r_var.clear();
        z_var.clear();
        local_k.clear();
        
        for (double r = r0; r < r1; r += 0.001) {
            itheta = dCenterBeta * dth_decimal;
            r_var.push_back(r);
            z_var.push_back(0);
            getAveFields();
            vec_D_Bz.push_back(fabs(fluz));
            dMag << r_var.back() << "," << z_var.back() << "," << (double)itheta / dth_decimal << "," << flux << "," << fluy << "," << fabs(fluz) << endl;
        }
        
        for (int i = aveNum; i <= r_var.size() - 2 - aveNum; i++) {
            drAve = 0;
            dBAve = 0;
            for (int j = i - aveNum; j <= i + aveNum; j++) {
                drAve += pow(r_var.at(j + 1) - r_var.at(j), 2);
                dBAve += pow(vec_D_Bz.at(j + 1) - vec_D_Bz.at(j), 2);
            }
            drAve = sqrt(drAve);
            dBAve = sqrt(dBAve);
            drAve = drAve / (2 * aveNum);
            dBAve = dBAve / (2 * aveNum);
            
            local_k.push_back(fabs(r_var.at(i) / vec_D_Bz.at(i) * dBAve / drAve));
        }

        
        for (int i = 0; i <= local_k.size() - 1; i++) {
            LKD << r_var.at(i) << "," << local_k.at(i) << endl;
        }
        
        r_var.clear();
        z_var.clear();
        local_k.clear();
    }
    else if (magnetTriplet == DFD) {
        dCenterBeta = betaStraight + betaD;
        fCenterBeta = dCenterBeta + betaD + betaGap + betaF / 2;
        
        cout << dCenterBeta * dth_decimal << endl;
        cout << fCenterBeta * dth_decimal << endl;
        for (double r = r0; r < r1; r += 0.001) {
            itheta = fCenterBeta * dth_decimal;
            r_var.push_back(r);
            z_var.push_back(0);
            getAveFields();
            vec_F_Bz.push_back(fabs(fluz));
            fMag << r_var.back() << "," << z_var.back() << "," << (double)itheta / dth_decimal << "," << flux << "," << fluy << "," << fabs(fluz) << endl;
        }
        
        for (int i = aveNum; i <= r_var.size() - 2 - aveNum; i++) {
            drAve = 0;
            dBAve = 0;
            for (int j = i - aveNum; j <= i + aveNum; j++) {
                drAve += pow(r_var.at(j + 1) - r_var.at(j), 2);
                dBAve += pow(vec_F_Bz.at(j + 1) - vec_F_Bz.at(j), 2);
            }
            drAve = sqrt(drAve);
            dBAve = sqrt(dBAve);
            drAve = drAve / (2 * aveNum);
            dBAve = dBAve / (2 * aveNum);
            
            local_k.push_back(fabs(r_var.at(i) / vec_F_Bz.at(i) * dBAve / drAve));
        }
        
        for (int i = 0; i <= local_k.size() - 1; i++) {
            LKF << r_var.at(i) << "," << local_k.at(i) << endl;
        }
        
        r_var.clear();
        z_var.clear();
        local_k.clear();
        
        for (double r = r0; r < r1; r += 0.001) {
            itheta = dCenterBeta * dth_decimal;
            r_var.push_back(r);
            z_var.push_back(0);
            getAveFields();
            vec_D_Bz.push_back(fabs(fluz));
            dMag << r_var.back() << "," << z_var.back() << "," << (double)itheta / dth_decimal << "," << flux << "," << fluy << "," << fabs(fluz) << endl;
        }
        
        for (int i = aveNum; i <= r_var.size() - 2 - aveNum; i++) {
            drAve = 0;
            dBAve = 0;
            for (int j = i - aveNum; j <= i + aveNum; j++) {
                drAve += pow(r_var.at(j + 1) - r_var.at(j), 2);
                dBAve += pow(vec_D_Bz.at(j + 1) - vec_D_Bz.at(j), 2);
            }
            drAve = sqrt(drAve);
            dBAve = sqrt(dBAve);
            drAve = drAve / (2 * aveNum);
            dBAve = dBAve / (2 * aveNum);
            
            local_k.push_back(fabs(r_var.at(i) / vec_D_Bz.at(i) * dBAve / drAve));
        }
        
        
        for (int i = 0; i <= local_k.size() - 1; i++) {
            LKD << r_var.at(i) << "," << local_k.at(i) << endl;
        }
        
        r_var.clear();
        z_var.clear();
        local_k.clear();
    }
    
}

// changed for 22.5
void halfCellMagField(double r, double z) {
    ofstream field ("halfCellField.csv");
    
    for (int ith = - 360 * dth_decimal / symmetryNumber / 2; ith <= 360 * dth_decimal / symmetryNumber / 2; ith += 1) {
        if (ith % 10 != 0) continue;
        itheta = ith;
        r_var.push_back(r);
        z_var.push_back(z);
        getAveFields();
        field << r_var.back() << "," << z_var.back() << "," << (double)itheta / dth_decimal << ","
        << flux << "," << fluy << "," << fluz << endl;
    }
    r_var.clear();
    z_var.clear();
}
