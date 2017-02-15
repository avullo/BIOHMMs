#include "General.h"
#include "DataSet.h"
#include <cassert>
#include <cstdlib>
#include <string>
using namespace std;

DataSet::DataSet(int the_length): totSize(0), length(the_length), instances(new Instance*[length]), pos(new int[length]) {}

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
  
  instances = new Instance*[length];
  pos = new int[length];
  for (int p=0; p<length; ++p) {
    instances[p] = new Instance(is, inputSymbols, outputSymbols, quot);
    pos[p] = p;
    
    totSize += instances[p]->length;
  }
  
  assert(totSize > 0);
};

DataSet::~DataSet() {
  for(int p=0; p<length; ++p)
    if(instances[p] != NULL)
      delete instances[p];
  if(instances != NULL)
    delete[] instances;
  if(pos != NULL)
    delete[] pos;

  delete inputSymbols;
  delete outputSymbols;
}

int DataSet::getInputDim() const { return inputSymbols->size(); }

int DataSet::getOutputDim() const { return outputSymbols->size(); }

void DataSet::shuffle() {
  // Shuffle training set positions
  for (int k=0; k<size(); ++k) {
    int f1 = rand();
    int f2 = rand();
    int p1 = (int)((double)f1/(1.0+(double)(RAND_MAX))*size());
    int p2 = (int)((double)f2/(1.0+(double)(RAND_MAX))*size());
    
    int tmp = pos[p1];
    pos[p1] = pos[p2];
    pos[p2] = tmp;
  }
}

void DataSet::write(ostream& os) {
  os << length << endl;
  for(int p=0; p<length; ++p) {
    instances[p]->write(os);
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
    instances[p]->write_predictions(os);
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


