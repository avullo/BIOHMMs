#ifndef MODEL_H
#define MODEL_H

#include "Instance.h"
#include "Parameterisation.h"
#include <cstdlib>
#include <cmath>
#include <string>
#include <iostream>

#define MAX 2500

class Model {
 public:
  Model();
  Model(int, int, int, int, const std::string& = "MT");
  Model(std::istream& is);
  ~Model();
  
  void read(std::istream&);
  void write(std::ostream&);

  int inputDim() const { return NI; }
  int outputDim() const { return NO; }
  int forwardStateDim() const { return NF; }
  int backwardStateDim() const { return NB; }

  void randomize(int seed);
  void e_step(Instance*);
  void m_step(Float = 0.1, Float = .05);

  void predict(Instance*);
  int getClasses() const { return NO; }
  Float getError() const { return error; }
  void resetError() { error = .0; }

 private:
  /*
    Input, output and forward/backward states are multinomial random variables
    Store the number of values for each of them.
  */
  int NI; // number of values of I_t (random variable representing input at time t)
  int NO; // the same for O_t (output at time t)
  int NF; // number of realisations of F_t (forward state variable at time t)
  int NB; // the same for B_t (backward state variable at t)

  // the output after E-step (to check the error) and prediction
  Float* Y;
  double error;

  /*
    the underlying implementation of the model parameters, i.e.
    P(O_t|F_t, B_t, I_t), P(F_t|F_t-1, I_t), P(B_t|B_t+1, I_t)

   set fixed to multinomial tables until I develop a generic framework
  */
  std::string ptype;
  Parameterisation* model; 

  /*
    Junction tree cluster (i.e. clique) potentials
    Cliques are (I_t, F_t, B_t, O_t), (I_t, B_t, F_t, F_t-1), (I_t, F_t, B_t, B_T+1)
  */
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

  void allocStats();
  void deallocStats();
  void resetStats();
  void sufficientStats(int);
  
  // juction tree routines
  void tree_alloc(int);
  void tree_dealloc(int);

  void attachCPT(int);
  void attachPOFBI(int);
  void attachPFFI(int);
  void attachPBBI(int);
  void attachPF();
  void attachPB(int T);

  void injectIn(int*, int);
  void injectOut(int*, int);
  
  void propagate(int);
  void saveOutput(int);
  
};

#endif // MODEL_H
