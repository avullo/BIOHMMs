/*

  A class to represent bidirectional IOHMMs

 */

#ifndef BIOHMM_H
#define BIOHMM_H

#include <cstdlib>
#include <cmath>
#include <iostream>

#define MAX 2500


// BIOHMM ver. 1.0 (16/9/2005)
// Copyright (C) Gianluca Pollastri 2005

class BIOHMM {
 public:

  BIOHMM(int, int, int, int);
  BIOHMM(std::istream& is);
  ~BIOHMM();
  void read(std::istream& is);
  void write(std::ostream& os);
  void propagate(int);
  void sufficientStats(int);

  void extimation(int*, int*, int);
  void maximization(Float = .1, Float = .05);
  void predict(int*, int);

  Float* out() { return Y; }
  Float getError() { return error; }
  void resetError() { error =.0; }

 private:
  /*
    Input, output and forward/backward states are multinomial random variables
    Store the number of values for each of them.
   */
  int NI; // number of values of I_t (random variable representing input at time t)
  int NO; // the same for O_t (output at time t)
  int NF; // number of realisations of F_t (forward state variable at time t)
  int NB; // the same for B_t (backward state variable at t)

  /*
    Model parameters

    The parameters of a BN specify the local conditional distribution of each
    variable given its parents. In this case, the distributions are:

    P(O_t|F_t, B_t, I_t)
    P(F_t|F_t-1, I_t)
    P(B_t|B_t+1, I_t)

    In the discrete case, these can be explicitly represented using conditional
    propability tables (CPTs). 
    The assumption here is that the model is stationary, i.e. the above CPTs do
    not vary over time, which can be seen as a form of parameter sharing which
    significantly reduces the degrees of freedom.
    
    NOTE: the unconditional distribution of root nodes (i.e. P(U_t)) do not need 
    to be modeled since we assume there are no missing data in the input sequence.
  */
  Float**** P_OFBI;
  Float*** P_FFI;
  Float*** P_BBI;

  // To be complete, the model needs to specify the distributions of the boundary
  // varialbles, i.e. P(F_0), P(B_T+1)
  Float* P_F;
  Float* P_B;

  // Junction tree cluster (i.e. clique) potentials
  // Cliques are (I_t, F_t, B_t, Y_t), (I_t, B_t, F_t, F_t-1), (I_t, F_t, B_t, B_T+1)
  Float***** OFBI;
  Float***** FFBI;
  Float***** BBFI;

  // junction tree separators potentials
  Float**** FBIa;
  Float**** FBIb;
  Float*** FB1;

  Float* mua;
  Float* mub;
  Float* mu1;

  /*
    Learning is formulated in the context of maximum likelihood (ML) and is
    solved by the expectation maximisation (EM) algorithm.

    In the case of belief networks, the E-step consists of computing the 
    expected sufficient statistics for the parameters. In the case of multinomial
    tables, the M-step consists in updating each entry of a table by the
    corresponding normalised expected counts.
    
    These tables store the expected counts.
   */
  Float**** P_OFBIss;
  Float*** P_FFIss;
  Float*** P_BBIss;

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
  
};


#endif // BIOHMM_H
