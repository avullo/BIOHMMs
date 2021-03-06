#include "Model.h"
#include "DataSet.h"

#include <cassert>
#include <cmath>
#include <cfloat>
#include <iostream>
using namespace std;

class Options {
public:
  Float epsilon;
  int NF;
  int NB;
  int seed;
  int shuffle;
  int batch_blocks;
  int readModel;
  int readEpoch;
  int nEpochs;

  int adaptive;
  int reload;
  Float belief;
  Float attenuation;

  void write(std::ostream& os = cout) {
    os << "epsilon " << epsilon << endl;
    os << "NF " << NF << endl;
    os << "NB " << NB << endl;
    os << "seed " << seed << endl;
    os << "shuffle " << shuffle << endl;
    os << "batch_blocks " << batch_blocks << endl;
    os << "readModel " << readModel << endl;
    os << "readEpoch " << readEpoch << endl;
    os << "nEpochs " << nEpochs << endl;
    os << "adaptive " << adaptive << endl;
    os << "reload " << reload << endl;
    os << "belief " << belief << endl;
    os << "attenuation " << attenuation << endl;
  }
  
  Options(std::istream& is) {
    char str[128];
    while(is) {
      is >> str;
      if(!strcmp(str, "epsilon")) is >> epsilon;
      else if(!strcmp(str, "NF")) is >> NF;
      else if(!strcmp(str, "NB")) is >> NB;
      else if(!strcmp(str, "seed")) is >> seed;
      else if(!strcmp(str, "shuffle")) is >> shuffle;
      else if(!strcmp(str, "batch_blocks")) is >> batch_blocks;
      else if(!strcmp(str, "readModel")) is >> readModel;
      else if(!strcmp(str, "readEpoch")) is >> readEpoch;
      else if(!strcmp(str, "nEpochs")) is >> nEpochs;
      else if(!strcmp(str, "adaptive")) is >> adaptive;
      else if(!strcmp(str, "reload")) is >> reload; 
      else if(!strcmp(str, "belief")) is >> belief; 
      else if(!strcmp(str, "attenuation")) is >> attenuation; 
    }
  }
};

class Performance {
  int NY;
  int** _cf; // confusion matrix

  string msg;
public:
  Performance(int nclasses, const string& message = ""): NY(nclasses), msg(message) {
    _cf = new int*[NY];
    for(int y=0; y<NY; ++y)
      _cf[y] = new int[NY];
  }
  ~Performance() {
    for(int y=0; y<NY; ++y)
      delete[] _cf[y];
    delete[] _cf;
  }

  void update(Instance* instance) {
    for(int t=1; t<=instance->length; ++t) {
      assert(instance->y[t] != -1 && instance->y_pred[t] != -1);
      _cf[instance->y_pred[t]][instance->y[t]]++;
    }
  }
  
  void reset() {
    for(int y1=0; y1<NY; ++y1)
      for(int y2=0; y2<NY; ++y2)
	_cf[y1][y2] = 0;
  }

  int getNumInstancesOfClass(int c) const {
    assert(c >= 0 && c < NY);

    int count = 0;
    for(int y=0; y<NY; ++y)
      count += _cf[y][c];

    return count;
  }

  int getNumErrorsForClass(int c) const {
    assert(c >= 0 && c < NY);
    
    int count = 0;
    for(int y=0; y<NY; ++y)
      if(y != c)
	count += _cf[y][c];

    return count;
  }
  
  int getTotalNumErrors() const {
    int count = 0;
    for(int y1=0; y1<NY; ++y1)
      for(int y2=0; y2!=NY; ++y2) 
	if(y1 != y2)
	  count += _cf[y1][y2];

    return count;
    
  }

  friend ostream& operator<<(ostream& os, const Performance& p) {
    double a[128];
    double all = 0;
    for(int y=0; y<p.NY; ++y) {
      a[y] = p.getNumInstancesOfClass(y);
      all += a[y];
    }

    int tot_class_errors = p.getTotalNumErrors();
    os << endl << endl << p.msg << "_NErrors= " << tot_class_errors << "/" << all;
    os << " " << (double)tot_class_errors/(double)(all) * 100.0;

    for(int y=0; y<p.NY; ++y) {
      int tot_errors_per_class = p.getNumErrorsForClass(y);
      os << endl << "Class" << y <<"= " << tot_errors_per_class << "/" << a[y];
      os << "\t" << (double)tot_errors_per_class/(double)a[y] * 100.0;
    }
  }
  
};

int Errcomp=100000;

void save(int epoch, Model* M) {
  filebuf outbuf;
  char fname[1024];
  sprintf(fname, "trained-%d.model", epoch);
  
  if(outbuf.open(fname, ios::out) != 0) {
    ostream os(&outbuf);
    M->write(os);
  } else {
    FAULT("Failed to write to file " << fname);
  }
  
  outbuf.close();
}

void load(int epoch, Model* M) {
  filebuf inbuf;
  char fname[1024];
  sprintf(fname, "trained-%d.model", epoch);
  
  if(inbuf.open(fname, ios::in) != 0) {
    istream is(&inbuf);
    M->read(is);
  } else {
    FAULT("Failed to read file " << fname);
  }
  
  inbuf.close();
}

void evaluate(Model* M, DataSet& D, char* which) {
  cout << endl << " counting_" << which << "_errors" << flush;

  Performance perf(M->getClasses(), which);
  perf.reset();
  for(int p=0; p<D.size(); ++p) {
    M->predict(D[p]);
    perf.update(D[p]);
    
    if(p%20==0) cout << "." << flush;
  }

  cout << perf;

  int tot_errors = perf.getTotalNumErrors();
  if((strncmp(which, "test", 4)==0) && tot_errors < Errcomp) {
    save(-10,M);
    Errcomp = tot_errors;
  }

  cout << endl;
}

void train(Model* M, DataSet& D, DataSet& T, Options& Opt) {
  int Gui = Opt.adaptive;
  //Number of steps at increasing error before
  //rescaling the learning rate.
  int gui=0;

  /*
    Float ep=Opt.epsilon/(Float)D.totSize;
    Float ep0=ep;
    if (Opt.batch_blocks>1) {
    ep *= (Float)(Opt.batch_blocks-1);
    }
    cout << "Actual lrate= " << ep << "\n";
    M->setEpsilon(ep);
  */

  cout << "Start learning" << endl;
  srand(Opt.seed);

  Float previous_error = FLT_MAX;  
  for(int epoch=Opt.readEpoch+1; epoch<=Opt.readEpoch + Opt.nEpochs; ++epoch) {
    if(Opt.shuffle)
      D.shuffle();

    M->resetError();
    int batch_cnt = 0;
    for(int pp=0; pp<D.size(); ++pp) {
      M->e_step(D[pp]);

      batch_cnt++;
      if(batch_cnt >= D.size()/Opt.batch_blocks && D.size()-pp >= D.size()/Opt.batch_blocks) {
	M->m_step(Opt.attenuation, Opt.belief);
	batch_cnt = 0;
      }
      if(pp%20 == 0) cout << "." << flush;
    }
    
    if(batch_cnt > 0)
      M->m_step(Opt.attenuation, Opt.belief);

    Float current_error = M->getError();
    cout << "\nEpoch " << epoch << " Error= " << current_error;

    save(0,M);
    if (current_error < previous_error) {
      gui=0;
      save(0,M);
      if (Gui>0) {
	//      ep += ep0*0.01;
	//      M->setEpsilon(ep);
      }
      previous_error = current_error;
    }

    if (epoch && epoch%5==0) {
      save(epoch, M);
      evaluate(M, D, "train");
      evaluate(M, T, "test");
      D.write("train-predictions");
      T.write("test-predictions");
    }
    cout << endl << flush;
  }
}

int main(int argc, char** argv) {
  if (argc<2) {
    cerr << "Usage: " << argv[0] << " option-file\n";
    exit(1);
  }
  ifstream optstream(argv[1]);
  Options Opt(optstream);
  Opt.write();

  // TODO
  // training/test sets file names should be given as arguments,
  // consider default as well
  cout << "Reading train dataset" << endl << flush;
  ifstream dstream("train.dataset");
  DataSet trainingSet(dstream);
  // trainingSet.set_belief(Opt.belief);
  
  cout << "Reading test dataset " << endl << flush;
  ifstream tstream("test.dataset");
  DataSet testSet(tstream);
  // testSet.set_belief(Opt.belief);

  int inputDim = trainingSet.getInputDim();
  int outputDim = trainingSet.getOutputDim();
  assert(inputDim == testSet.getInputDim());
  assert(outputDim == testSet.getOutputDim());
  
  Model* M;
  if (Opt.readModel) {
    char tmp[1024];
    sprintf(tmp, "trained-%d.model", Opt.readEpoch);
    cout << "Reading model from " << tmp << "\n";
    ifstream mstream(tmp);
    M = new Model(mstream);
  } else {
    cout << "Creating model" << endl << flush;

    M = new Model(inputDim, outputDim, Opt.NF, Opt.NB);

    cout << "Generating random parameters" << endl << flush;
    M->randomize(Opt.seed);
    save(-1, M);
    Opt.readEpoch = 0;
  }
  
  train(M, trainingSet, testSet, Opt);

  delete M;
  
  return 0;
}


