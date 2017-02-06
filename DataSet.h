#ifndef DATASET_H
#define DATASET_H

#include "Alphabet.h"
#include "Instance.h"
#include <iostream>

class DataSet {
 public:
  int length;
  Alphabet* inputSymbols;
  Alphabet* outputSymbols;
  Instance** seq;
  int totSize;

  DataSet() {};
  DataSet(int the_length);
  DataSet(std::istream& is, int quot = 0);
  ~DataSet();

  int getInputDim();
  int getOutputDim();
  void write(std::ostream& os);
  void write(char* fname);
  void write_predictions(std::ostream& os);
  void write_predictions(char* fname);
};


#endif // DATASET_H
