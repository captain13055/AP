//
//  run.cpp
//  AP
//
//  Created by yuya on 2015/09/15.
//  Copyright (c) 2015年 hy. All rights reserved.
//

#include "run.h"

void run() {
  int totalTurn, cellPosition, iE, eE;
  double dr, dtheta, dz;
  double rFieldCheck, zFieldCheck;
  clock_t t1, t2;
  vector<string> vinp, vpar, vfunc, value, function, vcont, vfcont;
  ifstream ifs("ap.inp");
  string is, sdata;
  
  int skin, siE, seE, stn;
  double sir, sdr, sdz, sdth, sdthbl, sps, smfr, smfz, scod;
  string sfmt, sfmn;

  int ihalf, iECOAccel, iFCO, iPCO, iPCOD, iPPS, iCTune;
  
  sdata = "";

  while (getline(ifs, is)) {
    sdata += is + "\n";
  }
  sdata.erase(--sdata.end());
  split(sdata, "#", vinp);
  vinp.at(0).erase(--vinp.at(0).end());
  vinp.at(1).erase(vinp.at(1).begin());
  split(vinp.at(0), "\n", vpar);
  split(vinp.at(1), "\n", vfunc);
  
  for (string s: vpar) {
    split(s, ":", value);
    DeleteSpace(value.back());
    vcont.push_back(value.back());
  }

  for (string s: vfunc) {
    split(s, ":", function);
    DeleteSpace(function.back());
    vfcont.push_back(function.back());
  }

  skin = stoi(vcont.at(0));
  sfmt = vcont.at(1);
  sfmn = vcont.at(2);
  sir = stod(vcont.at(3));
  sdr = stod(vcont.at(4));
  sdz = stod(vcont.at(5));
  sdth = stod(vcont.at(6));
  sdthbl = stod(vcont.at(7));
  siE = stoi(vcont.at(8));
  seE = stoi(vcont.at(9));
  stn = stoi(vcont.at(10));
  sps = stod(vcont.at(11));
  smfr = stod(vcont.at(12));
  smfz = stod(vcont.at(13));
  scod = stod(vcont.at(14));

  ihalf = stoi(vfcont.at(0));
  iECOAccel = stoi(vfcont.at(1));
  iFCO = stoi(vfcont.at(2));
  iPCO = stoi(vfcont.at(3));
  iPCOD = stoi(vfcont.at(4));
  iPPS = stoi(vfcont.at(5));
  iCTune = stoi(vfcont.at(6));
    
  t1 = clock();
  /*******************************/
  /*******************************/
  /**** CALCULATION CONDITION ****/
  // fieldMap means which field do you want?
  // Look at getData() in data.cpp
  if (sfmt == "Kyushu") {
    fieldMap = Kyushu;
  }
  else if (sfmt == "MERIT") {
    fieldMap = MERIT;
  }
  else if (sfmt == "Kyoto") {
    fieldMap = Kyoto;
  }
  else {
    cout << "No such a type" << endl;
    cout << "Choose" << endl;
    cout << "\"Kyushu\"" << endl;
    cout << "\"MERIT\"" << endl;
    cout << "\"Kyoto\"" << endl;
  }
  filename = sfmn;
  // useToscaField should be always true.
  useToscaField = true;
  
  // essential condition
  // kinetic energy.
  kinetic = skin;// [MeV]
  momentum = sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic);
  // initialCondition 0 ~ 5 = {r [m], t [sec], z[m], pr [MeV/c], pth [MeV/c], pz [MeV/c]}
  if (useToscaField) {
    initialCondition[0] = sir;
  }
  else {
    initialCondition[0] = sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic) / ( 300 * 1.0);
  }
  initialCondition[1] = 0; //t
  initialCondition[2] = 0; //z
  initialCondition[3] = 0; //pr
  initialCondition[4] = SIMom(sqrt(2 * NatMass(mass) * kinetic + kinetic * kinetic)); //pth
  initialCondition[5] = 0; //pz
  //sub essential condition
  // dr [m], dz [m] mean deviation from closed orbit.
  dr = sdr;
  dz = sdz;
  // dtheta [deg] -> each step in runge-kutta, converted into rad in calculation.
  dtheta = sdth;
  // In calculation of tune shift, iE and eE are used.
  iE = siE;
  eE = seE;
  // totalTurn -> How many turns do you want?
  totalTurn = stn;
  //sub non-essential condition
  // cellPosition -> later.
  cellPosition = sps;
  // r, zFieldCheck -> youcand check mag field at these points.
  rFieldCheck = smfr;
  zFieldCheck = smfz;
  //COD KICK r'
  kick = scod; //mrad
  /****     END     CONDITION ****/
  /*******************************/
  /*******************************/
  
  /************************************************************/
  /************************************************************/
  /************************************************************/
  /*******               MAIN CALCULATION               *******/
  cout << "1, kinetic energy  :" << kinetic << endl
       << "2, field map type  :" << sfmt << endl
       << "3, field map name  :" << sfmn << endl
       << "4, initial r       :" << initialCondition[0] << endl
       << "5, dr              :" << dr << endl
       << "6, dz              :" << dz << endl
       << "7, dtheta          :" << dtheta << endl
       << "8, dtheta for BL   :" << sdthbl << endl
       << "9,  Accel init E   :" << iE << endl
       << "10, Accel end  E   :" << eE << endl
       << "11, totalTurn      :" << totalTurn << endl
       << "12, Look PS at     :" << cellPosition << endl
       << "13, Look MF at r   :" << rFieldCheck << endl
       << "14, Look MF at z   :" << zFieldCheck << endl
       << "15, COD Kick       :" << kick << endl;
  cout << "Start ? (Press enter to start, others to quit)" << endl;
  if (cin.get() != '\n') exit(0);

  cout << "******Start******" << endl;
  
  dtheta = sdthbl; // Only for BL
  cout << "Getting a data now ....." << endl;
  getData();
  cout << "Complete." << endl;
  cout << "Getiing a data size now ....." << endl;
  getSize(dtheta);
  cout << "Complete." << endl;
  initialize(initialCondition);
  cout << "Start Calculating." << endl;
  // Checking magnetic field.
  
  BL(dtheta, 4.3, 5.4); //500, 550 == REAL, MERIT
  dtheta = sdth;
  getSize(dtheta);
  
  //    rPowerField(r_min, r_max);
  //    rPowerField(r_min, r_max, angleDefocus, angleFocus, angleGap, angleStraight / 2);
  if (ihalf) {
    cout << "Start hallCellMagField" << endl;
    halfCellMagField(rFieldCheck, zFieldCheck);
    cout << "End" << endl << endl;
  }

  if (iECOAccel) {
    cout << "Start E vs CO radius by Accel" << endl;
    existCOchekingByAccel(initialCondition, dtheta, iE, eE, 10);
    cout << "End" << endl << endl;
  }
  if (iFCO) {
    cout << "Start Find Closed Orbit" << endl;
    findClosedOrbit(initialCondition, dtheta, cellPosition);
    cout << "End" << endl << endl;
  }
  if (iPCO) {
    cout << "Start Plot Closed Orbit" << endl;
    drawClosedOrbit(initialCondition, dtheta);
    cout << "End" << endl << endl;
  }
  if (iPCOD) {
    cout << "Start Plot Closed Orbit Distorsion" << endl;
    drawClosedOrbitDistortion(initialCondition, dtheta);
    cout << "End" << endl << endl;
  }
  if (iPPS) {
    cout << "Start Plot Phase Space" << endl;
    betatron(initialCondition, dtheta, dr, dz, totalTurn, cellPosition);
    cout << "End" << endl << endl;
  }
  // In not beta, the calculation is done by shifting dr, dz from closed orbit separately.
  // In beta, at the same time.
  if(iCTune) {
    cout << "Start Calculate Tunes iE to eE" << endl;
    //    tuneShiftCalculation(initialCondition, dtheta, dr, dz, iE, eE, totalTurn);
    tuneShiftCalculationBeta(initialCondition, dtheta, dr, dz, iE, eE, totalTurn);
    cout << "End" << endl << endl;
  }
  
  
  /*******             END MAIN CALCULATION             *******/
  /************************************************************/
  /************************************************************/
  /************************************************************/
  
  t2 = clock();
  cout << "finish!!" << endl;
  cout << "calculation time: " << (int)((double)(t2 - t1) / CLOCKS_PER_SEC) / 60 << " m "
       << (int)((double)(t2 - t1) / CLOCKS_PER_SEC) % 60 << "s" << endl;
}
