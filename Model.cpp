#include "Model.h"
using namespace std;

void Model::alloc() {
  // TODO initialise?!
  counted = new int[NY];
  nerrors_ = new int[NY];
	
  Conf = new int*[NY];
  for(int y=0; y<NY; ++y)
    Conf[y] = new int[NY];
  // TODO: initialise Conf?
}

Model::Model(int the_NF, int the_NB): NF(the_NF), NB(the_NB) {
  // these are hardwired for a secondary structure prediction problem
  // where the input is the AA or AA profile at each position and the
  // output is one-hot encoding of three classes
  // TODO: make it general
  NU = 27;
  NY = 3;
  
  model = new CPTParameterisation(NU, NY, NF, NB);

  alloc();
}

Model::Model(istream& is) {
  is >> NU >> NY;
  is >> NF >> NB;

  model = new CPTParameterisation(is);

  alloc();
}

Model::~Model() {
  for(int y=0; y<NY; ++y)
    delete[] Conf[y];
  delete[] Conf;

  delete model;
}

void Model::read(istream& is) {
  is >> NU >> NY;
  is >> NF >> NB;

  model->read(is);

  threshold = .0;
}

void Model::write(ostream& os) {
  os << NU <<  " " << NY << endl;
  os << NF << " " << NB << endl;

  model->write(os);
}

void Model::randomize(int seed) {
  //Net->initWeights(seed);
  //NetF->initWeights(seed);
}

void Model::extimation(Sequence *seq) {
  int* I;
  int* O;

  I = new int[seq->length+1];
  O = new int[seq->length+1];

  for(int t=1; t<=seq->length; ++t) {
    I[t] = seq->u[t];
    O[t] = seq->y[t];
  }

  model->extimation(I, O, seq->length);

  delete[] I;
  delete[] O;
}

void Model::maximization(Float att, Float prior) {
  model->maximization(att, prior);
}

void Model::predict(Sequence* seq) {
  int* I;
  int* O;

  I = new int[seq->length+1];
  O = new int[seq->length+1];

  for(int t=1; t<=seq->length; ++t) {
    I[t] = seq->u[t];
    O[t] = seq->y[t];
  }

  model->predict(I, seq->length);

  for(int t=1; t<=seq->length; ++t) {
    Float pred = .0;
    int arg = -1;

    for(int c=0; c<NY; ++c) {
      if (model->out()[NY*t+c] > pred) {
	pred = model->out()[NY*t+c];
	arg = c;
      }
    }
    seq->y_pred[t] = arg;
  }

  for(int t=1; t<=seq->length; ++t) {
    // cout << t << ":\t" << seq->y[t] << ',' << seq->y_pred[t] << endl;
    if (seq->y[t] != -1 && seq->y_pred[t] != -1) {
      Conf[seq->y_pred[t]][seq->y[t]]++;
      counted[seq->y[t]]++;
    } else { 
      cout << seq->y_pred[t] << " " << seq->y[t] << endl << flush;
    }

    if (seq->y[t] != seq->y_pred[t]) {
      ++nerrors;
      nerrors_[seq->y[t]]++;
    }
  }
  
  delete[] I;
  delete[] O;
}

void Model::resetNErrors() { 
  nerrors = 0;
  memset(nerrors_, 0, NY*sizeof(int));
  memset(counted, 0, NY*sizeof(int));
  
  for(int p=0; p<NY; ++p)
    for(int y=0; y<NY; ++y)
      Conf[p][y] = 0;

  model->resetError();
}
