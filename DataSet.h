#ifndef DATASET_H
#define DATASET_H

#include "Alphabet.h"
#include "Instance.h"
#include <cassert>
#include <iostream>

class DataSet {
  int length;

  Alphabet* inputSymbols;
  Alphabet* outputSymbols;

  Instance** instances;
  int* pos;
  int totSize;

 public:
  DataSet() {};
  DataSet(int the_length);
  DataSet(std::istream& is, int quot = 0);
  ~DataSet();

  int size() const { return length; }
  int getInputDim() const;
  int getOutputDim() const;
  Instance* operator[](int i) {
    assert(i >= 0 && i < length);
    return instances[pos[i]];
  }

  void shuffle();
  
  void write(std::ostream& os);
  void write(char* fname);
  void write_predictions(std::ostream& os);
  void write_predictions(char* fname);
};


#endif // DATASET_H
