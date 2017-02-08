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

void Model::extimation(Instance* instance) {
  model->extimation(instance->u, instance->y, instance->length);
}

void Model::maximization(Float att, Float prior) {
  model->maximization(att, prior);
}

void Model::predict(Instance* instance) {
  model->predict(instance->u, instance->length);

  for(int t=1; t<=instance->length; ++t) {
    Float pred = .0;
    int arg = -1;

    for(int c=0; c<NY; ++c) {
      if (model->out()[NY*t+c] > pred) {
	pred = model->out()[NY*t+c];
	arg = c;
      }
    }
    instance->y_pred[t] = arg;
  }
}
