#ifndef MODEL_H
#define MODEL_H

#include <cstdlib>
#include <cmath>
#include "Sequence.h"
#include "BIOHMM.h"


class Model {
  int NU;
  int NY;

  int NF;
  int NB;

  Float threshold;

  BIOHMM* Net;

  int** Conf;

  Float temp_error;
  int temp_aas;
  
  int* counted;
  double squared_error;
  int nerrors;
  int* nerrors_;

  Float epsilon;

  void alloc();

 public:

  Model(int NF, int NB);
  Model(std::istream& is);

  void read(std::istream& is);
  void write(std::ostream& os);
  
  void randomize(int seed);
  void extimation(Sequence* seq);
  void maximization(Float att = 0.1, Float prior = .05);

  void predict(Sequence* seq);
  Float* out() { return Net->out(); } 
  int** getConf() { return  Conf; }
  int getNErrors() { return nerrors; }
  int getNErrors_(int i) { return nerrors_[i]; }
  int getClasses() { return NY; }
  int* getCounted() { return counted; }
  void resetNErrors();
  Float get_tot_error() { return temp_error; }
  Float get_squared_error() { return Net->getError(); }
  void reset_squared_error() { Net->resetError(); }
  
};

#endif // MODEL_H
