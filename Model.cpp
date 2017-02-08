#include "Model.h"
using namespace std;

Model::Model(int the_NU, int the_NY, int the_NF, int the_NB): NU(the_NU), NY(the_NY), NF(the_NF), NB(the_NB) {
  model = new CPTParameterisation(NU, NY, NF, NB);
}

Model::Model(istream& is) {
  is >> NU >> NY;
  is >> NF >> NB;

  model = new CPTParameterisation(is);
}

Model::~Model() {
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

void Model::extimation(Instance *seq) {
  int* I;
  int* O;

  // NOTE
  // Here it makes a copy of the input/output sequence
  // which seems a waste of time/space
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

void Model::predict(Instance* seq) {
  int* I;
  int* O;

  // NOTE
  // Here it makes a copy of the input/output sequence
  // which seems a waste of time/space
  // It doesn't also seem to be using the output sequence O
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
  
  delete[] I;
  delete[] O;
}
