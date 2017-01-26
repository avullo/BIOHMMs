#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <iostream>
#include <cmath>
#include "General.h"

#define MAX_T 8196

static int translateY[26] = {
  -1,	// A
  -1,	// B
  2,	// C
  -1,	// D
  1,	// E
  -1,	// F
  -1,	// G
  0,	// H
  -1,	// I
  -1,	// J
  -1,	// K
  -1,	// L
  -1,	// M
  -1,	// N
  -1,	// O
  -1,	// P
  -1,	// Q
  -1,	// R
  -1,	// S
  -1,	// T
  -1,	// U
  -1,	// V
  -1,	// W
  -1,	// X
  -1,	// Y
  -1	// Z
};

static char Ytranslate[26] = "HEC";

class Sequence {
 public:

  char name[256];

  int* u;
  int* y;
  int* y_pred;
  Float* app;

  int length;

  int alignments_loaded;
  int** ALIGNMENTS;

  Float* HeP;
  int HePl;

  Sequence(std::istream& is, int quot = 0);

  int load_alignments(char* alidir, int dir = 0);
  int unload_alignments();
  Float profile_entropy();
  int generate_profile(char* alidir, int dir = 0);
  void unload_profile();
  void write(std::ostream& os);
  void write_predictions(std::ostream& os);
};

#endif // SEQUENCE_H
