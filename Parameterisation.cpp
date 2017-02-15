#include "Parameterisation.h"
using namespace std;

Parameterisation* Parameterisation::factory(int NI, int NO, int NF, int NU, const string& type)
  throw(BadParameterisationCreation) {
  if(type == "MT") { return new MTParameterisation(NI, NO, NF, NU); }

  throw BadParameterisationCreation(type);

}

Parameterisation* Parameterisation::factory(istream& is, const string& type) 
  throw(BadParameterisationCreation) {
  if(type == "MT") { return new MTParameterisation(is); }

  throw BadParameterisationCreation(type);

}
