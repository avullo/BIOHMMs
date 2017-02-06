#include "General.h"
#include "DataSet.h"
#include <cassert>
#include <cstdlib>
#include <string>
using namespace std;

DataSet::DataSet(int the_length): totSize(0), length(the_length), seq(new Instance*[length]) {}

DataSet::DataSet(istream& is, int quot): totSize(0), length(0) {
  is >> length;
  assert(length > 0);
  
  string inputAlphabet, outputAlphabet; 
  is >> inputAlphabet >> outputAlphabet;

  try {
    inputSymbols = Alphabet::factory(inputAlphabet);
    outputSymbols = Alphabet::factory(outputAlphabet);
  } catch(Alphabet::BadAlphabetCreation e) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }
  
  seq = new Instance*[length];
  for (int p=0; p<length; ++p) {
    seq[p] = new Instance(is, inputSymbols, outputSymbols, quot);
    totSize += seq[p]->length;
  }
  
  assert(totSize > 0);
};

DataSet::~DataSet() {
  for(int p=0; p<length; ++p)
    if(seq[p] != NULL)
      delete seq[p];
  if(seq != NULL)
    delete[] seq;

  delete inputSymbols;
  delete outputSymbols;
}

int DataSet::getInputDim() { return inputSymbols->size(); }

int DataSet::getOutputDim() { return outputSymbols->size(); }

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


