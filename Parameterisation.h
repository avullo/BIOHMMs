#ifndef PARAMETERISATION_H
#define PARAMETERISATION_H

#include "General.h"
#include <string>
#include <iostream>

/*

  Abtract class representing the underlying parameterisation of the
  conditional probabilities of a BIOHMMs, i.e.
  P(O_t|F_t, B_t, I_t), P(F_t|F_t-1, I_t), P(B_t|B_t+1, I_t)
 
  In concrete classes, can be either implemented with multinomial tables 
  or feed-forward neural networks.

 */

// forward declaration to declare Model friend
class Model;

class Parameterisation {
 public:
 Parameterisation(): NI(-1), NO(-1), NF(-1), NB(-1) {}
 Parameterisation(int idim, int odim, int fdim, int bdim):
  NI(idim), NO(odim), NF(fdim), NB(bdim) {}
  Parameterisation(std::istream&);
  virtual ~Parameterisation() { }

  virtual void read(std::istream&) = 0;
  virtual void write(std::ostream&) = 0;
  
  virtual void maximisation(Float****, Float***, Float***, Float = .1, Float = .05) = 0;

  virtual Float p_OFBI(int, int, int, int) = 0;
  virtual Float p_FFI(int, int, int) = 0;
  virtual Float p_BBI(int, int, int) = 0;

  class BadParameterisationCreation: public std::logic_error {
  public:
  BadParameterisationCreation(std::string type):
    logic_error("Cannot create type " + type) {}
  };

  // factory methods
  static Parameterisation* factory(int, int, int, int, const std::string&)
    throw(BadParameterisationCreation);
  static Parameterisation* factory(std::istream&, const std::string&)
    throw(BadParameterisationCreation);
  
  // to allow it to access the various dimensions
  // and the distribution of the boundary variables
  friend class Model;

 protected:
  int NI, NO, NF, NB;
  
  // To be complete, the model needs to specify the distributions of the boundary
  // varialbles, i.e. P(F_0), P(B_T+1)
  // Being protected, they can be accessed by concrete subclasses
  Float* P_F;
  Float* P_B;
  
};

/*

  Implements BIOHMMs conditional probabilities 
  with multinomial tables

 */
class MTParameterisation: public Parameterisation {
  MTParameterisation();
  MTParameterisation(int, int, int, int);
  MTParameterisation(std::istream&);
  friend class Parameterisation;
  
 public:
  ~MTParameterisation();

  void read(std::istream&);
  void write(std::ostream&);

  Float p_OFBI(int, int, int, int);
  Float p_FFI(int, int, int);
  Float p_BBI(int, int, int);

  void maximisation(Float****, Float***, Float***, Float = .1, Float = .05);
  
 private:
  /*
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

  void alloc();
  void initTables();
  void normaliseTables();
  void readTables(std::istream&);
  void writeTables(std::ostream&);

};

/*

  Implements BIOHMMs conditional probabilities 
  with feed-forward neural networks

 */
class NNParameterisation;

#endif // PARAMETERISATION_H
