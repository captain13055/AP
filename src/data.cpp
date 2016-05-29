//
//  data.cpp
//  AP
//
//  Created by yuya on 2015/08/02.
//  Copyright (c) 2015年 hy. All rights reserved.
//

#include "data.h"

map fieldMap;
triplet magnetTriplet;
int dataIndex, symmetryNumber;
int rDec, thDec, zDec, iRDelta, iThDelta, iZDelta;
int dth_decimal;
int r_size, th_size, z_size;
int ir_max, ir_min, iz_max, iz_min;
int angleInCell;
double r_max, r_min, z_max, z_min;
bool useToscaField;

double angleDefocus;
double angleGap;
double angleFocus;
double angleStraight;

string filename;

// [0] : x , [1] : y, [2] : z [cm, gauss]
vector<double> pos[3];
vector<double> mag[3];
vector<double> dat[6];
vector<double> cylin[2];   //[0]: r, [1]: theta [cm]


void getData() {
    char ch, cho;
    double x_tmp, y_tmp;
    coordinateSystem cs;
    dataIndex = 0;
    ifstream ifs;
    string str, trash, fieldMapName;
    stringstream ss;
    if (fieldMap == Kyushu) {
        symmetryNumber = 12;
        cs = xyz;
        magnetTriplet = DFD;
        angleDefocus = 1.715;
        angleGap = 1.7;
        angleFocus = 10.24;
        angleStraight = 9.5;
        fieldMapName = filename;
    }
    else if (fieldMap == MERIT) {
        symmetryNumber = 8;
        cs = xzy;
        magnetTriplet = FDF;
        angleDefocus = 4.1025;
        angleGap = 1.984;
        angleFocus = 11.643;
        angleStraight = 4.7705 * 2;
        fieldMapName = filename;
    }
    else if (fieldMap == Kyoto) {
        symmetryNumber = 12;
        cs = xyz;
        magnetTriplet = DFD;
        angleDefocus = 1.715;
        angleGap = 1.7;
        angleFocus = 10.24;
        angleStraight = 9.5;
        fieldMapName = filename;
    }
    else {
        cout << "incorrect field map was chosen." << endl;
        exit(0);
    }
    ifs.open(fieldMapName);
    
    
    if (ifs.fail()) {
        cerr << "File does not exist.\n";
        exit(0);
    }
    while (!ifs.eof()) {
        cho = ifs.get();
        if (cho == '\t') {
            cho = ' ';
        }
        ss << cho;
    }

    for (int i = 1; i <= 8; i++) {
        getline(ss, trash, '\n');
    }
    
    while (!ss.eof()) {
        dataIndex++;
        for (int i = 0; i < 5; i++) {
            while ((ch = ss.get()) == ' ');
            getline(ss, str, ' ');
            str = ch + str;
            try {
                stod(str);
            } catch (invalid_argument) {
	        break;
            }
            dat[i].push_back(stod(str));
        }
        while ((ch = ss.get()) == ' ');
        getline(ss, str, '\n');
        str = ch + str;
        try {
            stod(str);
        } catch (invalid_argument) {
            break;
        }
        dat[5].push_back(stod(str));
        
    }
    dataIndex--;
    if (cs == xzy) {
        for (int i = 0; i <= dataIndex - 1; i++) {
            pos[0].push_back(dat[0][i]);
            pos[1].push_back(dat[2][i]);
            pos[2].push_back(dat[1][i]);
            mag[0].push_back(-dat[3][i]);
            mag[1].push_back(-dat[5][i]);
            mag[2].push_back(-dat[4][i]);
        }
    }
    else if (cs == yxz) {
        for (int i = 0; i <= dataIndex - 1; i++) {
            pos[0].push_back(dat[1][i]);
            pos[1].push_back(-dat[0][i]);
            pos[2].push_back(dat[2][i]);
            mag[0].push_back(dat[4][i]);
            mag[1].push_back(-dat[3][i]);
            mag[2].push_back(dat[5][i]);
        }
    }
    else if (cs == xyz) {
		if (fieldMap == Kyushu) {
			for (int i = 0; i <= dataIndex - 1; i++) {
				pos[0].push_back(dat[0][i]);
				pos[1].push_back(dat[1][i]);
				pos[2].push_back(dat[2][i]);
				mag[0].push_back(dat[3][i]);
				mag[1].push_back(dat[4][i]);
				mag[2].push_back(dat[5][i]);
			}
        }
		else if(fieldMap == Kyoto) {
			for (int i = 0; i <= dataIndex - 1; i++) {
				pos[0].push_back(dat[0][i]);
				pos[1].push_back(dat[1][i]);
				pos[2].push_back(dat[2][i]);
				mag[0].push_back(-dat[3][i]);
				mag[1].push_back(-dat[4][i]);
				mag[2].push_back(-dat[5][i]);
			}
		}
    }
    else if (cs == rtz) {
        for (int i = 0; i <= dataIndex - 1; i++) {
            pos[2].push_back(dat[2][i]);
            mag[0].push_back(dat[3][i]);
            mag[1].push_back(dat[4][i]);
            mag[2].push_back(dat[5][i]);
        }
    }
    else {
        cout << "Cheat coordinate system." << endl;
        exit(0);
    }

    if (cs == xzy || cs == yxz || cs == xyz) {
        for (int i = 0; i <= dataIndex - 1; i++) {
            cylin[0].push_back(sqrt(pow(pos[0][i], 2) + pow(pos[1][i], 2)));
            cylin[1].push_back(atan(pos[1][i] / pos[0][i]) * 180 / M_PI);
            x_tmp = mag[0].at(i);
            y_tmp = mag[1].at(i);

            mag[0].at(i) = x_tmp * cos(cylin[1][i] * M_PI / 180) + y_tmp * sin(cylin[1][i] * M_PI / 180);
            mag[1].at(i) = -x_tmp * sin(cylin[1][i] * M_PI / 180) + y_tmp * cos(cylin[1][i] * M_PI / 180);
        
            if (fabs(pos[2][i]) < 1e-9) {
                mag[0].at(i) = 0;
                mag[1].at(i) = 0;
            }
        }
    }
    else if (cs == rtz) {
        for (int i = 0; i <= dataIndex - 1; i++) {
            cylin[0].push_back(dat[0][i]);
            cylin[1].push_back(dat[1][i]);
        }
    }
    
    r_max = cylin[0].back() * 1e-2;
    r_min = cylin[0].at(0) * 1e-2;
    z_max = pos[2].back() * 1e-2;
    z_min = -z_max;
    
    ir_min = (int)(r_min * 100 + 1e-7);
    ir_max = (int)(r_max * 100 + 1e-7);
    iz_min = (int)(z_min * 100 + 1e-7);
    iz_max = (int)(z_max * 100 + 1e-7);
    
    
}

void getSize(double dth) {
    double rDelta, thDelta, zDelta;
    bool r_in, th_in, z_in;
    r_in = false;
    th_in = false;
    z_in = false;
    for (int i = 0; i <= dataIndex - 1; i++) {
        if (fabs(cylin[0][i] - cylin[0][i+1]) > 1e-5 && !r_in) {
            r_size = i;
            r_size++;
            r_in = true;
        }
        
        if (fabs(pos[2][i] - pos[2][i+1]) > 1e-5 && !z_in) {
            z_size = i;
            z_size++;
            z_in = true;
        }
        
        if (fabs(cylin[1][i] - cylin[1][i+1]) > 1e-5 && !th_in) {
            th_size = i;
            th_size++;
            th_in = true;
        }
        
        if (r_in && z_in && th_in) {
            break;
        }
    }
    
    rDec = 1;
    thDec = 1;
    zDec = 1;
    dth_decimal = 1;
    
    rDelta = fabs(cylin[0][r_size] - cylin[0][0]);
    thDelta = fabs(cylin[1][th_size] - cylin[1][0]);
    zDelta = fabs(pos[2][z_size] - pos[2][0]);
    while ((dth - (1 - 1e-9) < 1e-10)) {
        dth_decimal *= 10;
        dth *= 10;
    }
    while ((rDelta - (1 - 1e-9) < 1e-10)) {
        rDec *= 10;
        rDelta *= 10;
    }
    while ((thDelta - (1 - 1e-9) < 1e-10)) {
        thDec *= 10;
        thDelta *= 10;
    }
    while ((zDelta - (1 - 1e-9) < 1e-10)) {
        zDec *= 10;
        zDelta *= 10;
    }
    iRDelta = (int)(rDelta + 1e-8);
    iThDelta = (int)(thDelta + 1e-8);
    iZDelta = (int)(zDelta + 1e-8);

}


void getAveFields() {
    int ir_0 = 0, ir_1 = 0;
    int iz_0[2], iz_1[2];
    int ith_0[4], ith_1[4];
    int rInit, thInit, zInit, rCur, thCur, zCur;
    int symHalfTheta, symTheta;
    bool inv_x_field = false, inv_y_field = false;
    double m_vol[8], m_par[3][8], m_field[3][8], rev_rate;
    double r_tmp, th_tmp, z_tmp;
    int ith_tmp;
    
    rInit = (int)(cylin[0].at(0) * rDec + 1e-7);
    thInit = (int)(cylin[1].at(0) * thDec + 1e-7);
    zInit = (int)(pos[2].at(0) * zDec + 1e-7);
    symTheta = 360 * dth_decimal / symmetryNumber;
    symHalfTheta = symTheta / 2;
    
    if (useToscaField) {
        if (r_var.back() > r_max || r_var.back() < r_min || z_var.back() > z_max || z_var.back() < z_min) {
            cout << "r: " << r_var.back() * 100 << " [cm]" << endl
            << "theta: " << (double)itheta / dth_decimal << endl
            << "z: " << z_var.back() * 100 << " [cm]" << endl;
            cout << "Divergence in getAveField" << endl;
            cout << "Magnet Field" << endl;
            cout << "r: " << r_min * 100 << " ~ " << r_max * 100 << " [cm]" << endl;
            cout << "z: " << z_min * 100 << " ~ " << z_max * 100 << " [cm]" << endl;
            exit(0);
        }
    }
    
    
    r_tmp = r_var.back();
    ith_tmp = itheta;
    th_tmp = (double)ith_tmp / dth_decimal;
    z_tmp = z_var.back();
    r_tmp *= 100;  // m to cm
    z_tmp *= 100;  // m to cm
    
    /****  theta を　-15 ~ 15に収める  ****/
    if (ith_tmp > symHalfTheta) {
        if (ith_tmp % symHalfTheta == 0 && ith_tmp % symTheta != 0) {
            ith_tmp = symHalfTheta;
        }
        
        else if (ith_tmp % symTheta== 0) {
            ith_tmp = 0;
        }
        
        else {
            ith_tmp = ith_tmp % symTheta;
            if (ith_tmp <= symTheta && ith_tmp > symHalfTheta){
                ith_tmp = ith_tmp - symTheta;
            }
        }
    }
    //thetaの最終判断 (0 ~ 15に収まる)
    if (-symHalfTheta <= ith_tmp && ith_tmp <= 0) {
        ith_tmp *= -1;
        inv_y_field = true;
    }
    /****                            ****/
    //z 対照
    if (z_tmp < 0) {
        z_tmp *= -1;
        inv_x_field = true;
    }
    
    rCur = (int)(r_tmp * rDec + 1e-7);
    thCur = ith_tmp;
    zCur = (int)(z_tmp * zDec + 1e-7);
    ir_0 = (rCur - rInit) / iRDelta  * r_size;
    ir_1 = ir_0 + r_size;
    iz_0[0] = ir_0 + (zCur - zInit) / iZDelta * z_size;
    iz_0[1] = iz_0[0] + z_size;
    iz_1[0] = iz_0[0] + r_size;
    iz_1[1] = iz_0[1] + r_size;
    
    if (itheta == symHalfTheta) {
        ith_0[0] = iz_0[0] + (thCur * thDec / dth_decimal - thInit) / iThDelta * th_size;
        ith_1[0] = ith_0[0] - 1;
        ith_0[1] = iz_0[1] + (thCur * thDec / dth_decimal - thInit) / iThDelta * th_size;
        ith_1[1] = ith_0[1] - 1;
        ith_0[2] = iz_1[0] + (thCur * thDec / dth_decimal - thInit) / iThDelta * th_size;
        ith_1[2] = ith_0[2] - 1;
        ith_0[3] = iz_1[1] + (thCur * thDec / dth_decimal - thInit) / iThDelta * th_size;
        ith_1[3] = ith_0[3] - 1;
    }
    else {
        ith_0[0] = iz_0[0] + (thCur * thDec / dth_decimal - thInit) / iThDelta * th_size;
        ith_1[0] = ith_0[0] + 1;
        ith_0[1] = iz_0[1] + (thCur * thDec / dth_decimal - thInit) / iThDelta * th_size;
        ith_1[1] = ith_0[1] + 1;
        ith_0[2] = iz_1[0] + (thCur * thDec / dth_decimal - thInit) / iThDelta * th_size;
        ith_1[2] = ith_0[2] + 1;
        ith_0[3] = iz_1[1] + (thCur * thDec / dth_decimal - thInit) / iThDelta * th_size;
        ith_1[3] = ith_0[3] + 1;
    }
    
    
    for (int i = 0; i <= 3; i++) {
        m_par[0][i] = cylin[0][ith_0[i]];
        m_par[1][i] = cylin[1][ith_0[i]];
        m_par[2][i] = pos[2][ith_0[i]];
        m_field[0][i] = mag[0][ith_0[i]];
        m_field[1][i] = mag[1][ith_0[i]];
        m_field[2][i] = mag[2][ith_0[i]];
    }
    for (int i = 0; i <= 3; i++) {
        m_par[0][i+4] = cylin[0][ith_1[i]];
        m_par[1][i+4] = cylin[1][ith_1[i]];
        m_par[2][i+4] = pos[2][ith_0[i]];
        m_field[0][i+4] = mag[0][ith_1[i]];
        m_field[1][i+4] = mag[1][ith_1[i]];
        m_field[2][i+4] = mag[2][ith_1[i]];
    }
    
    rev_rate = 0;
    for (int i = 0; i <= 7; i++) {
        m_vol[7-i] = fabs(r_tmp - m_par[0][i]) * m_par[0][i] * fabs(th_tmp - m_par[1][i]) * M_PI / 180 * fabs(z_tmp - m_par[2][i]);
        rev_rate += m_vol[7-i];
    }

    flux = 0;
    fluy = 0;
    fluz = 0;

    for (int i = 0; i <= 7; i++) {
        flux += m_vol[i] / rev_rate * m_field[0][i];
        fluy += m_vol[i] / rev_rate * m_field[1][i];
        fluz += m_vol[i] / rev_rate * m_field[2][i];
    }
    
    
    if (inv_x_field == true) {
        flux *= -1;
        fluy *= -1;
    }
    
    if (inv_y_field == true) {
        fluy *= -1;
    }
    
    //gauss to Tesla
    flux *= 1e-4;
    fluy *= 1e-4;
    fluz *= 1e-4;
    
    if (fabs(flux) < 1e-10) {
        flux = 0;
    }
    if (fabs(fluy) < 1e-10) {
        fluy = 0;
    }
    if (fabs(fluz) < 1e-10) {
        fluz = 0;
    }
    
    
}
