#ifndef INSTANCE_H
#define INSTANCE_H

#include "General.h"
#include "Alphabet.h"
#include <iostream>
#include <cmath>

#define MAX_T 8196

class Instance {
 public:
  char name[256];

  int* u;
  int* y;
  int* y_pred;
  Float* app;

  int length;

  Alphabet *inputSymbols, *outputSymbols;
  
  int alignments_loaded;
  int** ALIGNMENTS;

  Float* HeP;
  int HePl;

  Instance(std::istream&, Alphabet*, Alphabet*, int = 0);
  ~Instance();
  
  int load_alignments(char* alidir, int dir = 0);
  int unload_alignments();
  Float profile_entropy();
  int generate_profile(char* alidir, int dir = 0);
  void unload_profile();

  // TODO: use Alphabet decode
  void write(std::ostream& os);
  void write_predictions(std::ostream& os);
};

#endif // INSTANCE_H
