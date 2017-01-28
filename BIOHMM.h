#ifndef BIOHMM_H
#define BIOHMM_H

#include <cstdlib>
#include <cmath>
#include <iostream>
#include "General.h"
//#include "NN.h"
//#include "NNt.h"

#define MAX 2500


// BIOHMM ver. 1.0 (16/9/2005)
// Copyright (C) Gianluca Pollastri 2005




class BIOHMM {
 private:

  int NI;
  int NO;
  int NF;
  int NB;

//  NN* NetOut;
//  NNt* NetF;
//  NNt* NetB;

  Float* mua;
  Float* mub;
  Float* mu1;

  // junction tree cluster potentials
  Float***** OFBI;
  Float***** FFBI;
  Float***** BBFI;

  // junction tree separators potentials
  Float**** FBIa;
  Float**** FBIb;
  Float*** FB1;

  // CPTs
  Float**** P_OFBI;
  Float*** P_FFI;
  Float*** P_BBI;

  // CPT sufficient stats?
  Float**** P_OFBIss;
  Float*** P_FFIss;
  Float*** P_BBIss;

  Float* P_F;
  Float* P_B;

  Float* Y;

  double error;

  void alloc();
  void initialiseTables();
  void normaliseTables();
  void resetTables();
  void resetStats();
  void printStats();
  void printTables();
  void readTables(std::istream& is);
  void writeTables(std::ostream& os);

  void tree_alloc(int length);
  void tree_dealloc(int length);

  void attachPOFBI(int t);
  void attachPFFI(int t);
  void attachPBBI(int t);
  void attachPF();
  void attachPB(int T);
  void attachCPT(int length);

  void injectOut(int* y, int length);
  void injectIn(int* x, int length);

  void saveOutput(int length);

//  int max(Float*a,int m);

 public:

  BIOHMM(int NI, int NO, int NF, int NB);
  BIOHMM(std::istream& is);
  void read(std::istream& is);
  void write(std::ostream& os);
  void propagate(int length);
  void sufficientStats(int length);

//  void extimation(Float* seq, int* y, int length, int backp=0);
  void extimation(int* seq, int* y, int length);

//  void extimation(Float* seq, Float* y, int length, int backp=0);
//  void extimation(int* seq, Float* y, int length, int backp=0);

  void maximization(Float att=0.1, Float prior=0.05);

//  void Feed(Float* seq, int length);
  void Feed(int* seq, int length);
//  void predict(Float* seq, int length);
  void predict(int* seq, int length);


  Float* out() {return Y;}

  Float getError() {
	return error;
  };
  void resetError() {
	error=0.0;
  };


};


#endif // BIOHMM_H
