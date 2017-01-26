#ifndef DATASET_H
#define DATASET_H

#include "Sequence.h"
#include <iostream>

class DataSet {
 public:
  int length;
  Sequence** seq;
  int totSize;

  DataSet() {};
  DataSet(int the_length);
  DataSet(std::istream& is, int quot = 0);

  void write(std::ostream& os);
  void write(char* fname);
  void write_predictions(std::ostream& os);
  void write_predictions(char* fname);
};


#endif // DATASET_H
