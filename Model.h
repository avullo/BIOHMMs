#ifndef MODEL_H
#define MODEL_H

#include <cstdlib>
#include <cmath>
#include "Instance.h"
#include "CPTParameterisation.h"

class Model {
  int NU;
  int NY;
  int NF;
  int NB;

  // the underlying model implementation,
  // set fixed to CPTs at the moment
  CPTParameterisation* model; 

  Float threshold;
  Float epsilon;
  
 public:

  Model(int NU, int NY, int NF, int NB);
  Model(std::istream& is);
  ~Model();
  void read(std::istream& is);
  void write(std::ostream& os);
  
  void randomize(int seed);
  void extimation(Instance* seq);
  void maximization(Float att = 0.1, Float prior = .05);

  void predict(Instance* seq);
  Float* out() { return model->out(); } 
  int getClasses() { return NY; }
  Float get_squared_error() { return model->getError(); }
  void reset_squared_error() { model->resetError(); }
  
};

#endif // MODEL_H
