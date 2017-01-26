#include "General.h"
#include "DataSet.h"
using namespace std;

DataSet::DataSet(int the_length): totSize(0), length(the_length), seq(new Sequence*[length]) {}

DataSet::DataSet(istream& is, int quot) {
  totSize = 0;
  is >> length;
  
  int foo;
  is >> foo >> foo;
  seq = new Sequence*[length];
  for (int p=0; p<length; ++p) {
    seq[p] = new Sequence(is, quot);
    totSize += seq[p]->length;
  }
};

void DataSet::write(ostream& os) {
  os << length << endl;
  for(int p=0; p<length; ++p) {
    seq[p]->write(os);
  }
}

void DataSet::write(char* fname) {
  filebuf outbuf;
  if (outbuf.open(fname, ios::out) != 0) {
    ostream os(&outbuf);
    this->write(os);
  } else {
    FAULT("Failed to write to file " << fname);
  }
  
  outbuf.close();
};


void DataSet::write_predictions(ostream& os) {
  os << length << endl;
  for (int p=0; p<length; ++p) {
    seq[p]->write_predictions(os);
  }
}

void DataSet::write_predictions(char* fname) {
  filebuf outbuf;
  if (outbuf.open(fname, ios::out) != 0) {
    ostream os(&outbuf);
    this->write_predictions(os);
  } else {
    FAULT("Failed to write to file " << fname);
  }
  
  outbuf.close();
}


