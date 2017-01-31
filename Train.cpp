#include "Model.h"
#include "DataSet.h"

#include <cmath>
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

void shuffle(DataSet& D, int* pos) {
  // Shuffle training set positions
  for (int k=0; k<D.length; k++) {
    //    int p1= (int)(drand48()*D.length);
    //    int p2= (int)(drand48()*D.length);
    int f1 = rand();
    int f2 = rand();
    int p1 = (int)((double)f1/(1.0+(double)(RAND_MAX))*D.length);
    int p2 = (int)((double)f2/(1.0+(double)(RAND_MAX))*D.length);
    int tmp = pos[p1];
    pos[p1] = pos[p2];
    pos[p2] = tmp;
  }
}

void evaluate(Model* M, DataSet& D, char* which) {
  cout << endl << " counting_" << which << "_errors" << flush;
  
  M->resetNErrors();
  for(int p=0; p<D.length; ++p) {
    M->predict(D.seq[p]);
    if(p%20==0) cout << "." << flush;
  }

  double a[128];
  double all = 0;
  for(int y=0; y<M->getClasses(); ++y) {
    a[y] = M->getCounted()[y];
    all += a[y];
  }

  cout << endl << endl << which << "_NErrors= " << M->getNErrors() << "/" << all;
  cout << " " << (double)M->getNErrors()/(double)(all) * 100.0;

  for(int y=0; y<M->getClasses(); ++y) {
    cout << endl << "Class" << y <<"= " << M->getNErrors_(y) << "/" << a[y];
    cout << "\t" << (double)M->getNErrors_(y)/(double)a[y] * 100.0;
  }

  if((strncmp(which, "test", 4)==0) && (M->getNErrors()<Errcomp)) {
    save(-10,M);
    Errcomp = M->getNErrors();
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
  
  int* wait = new int[D.length];
  for(int p=0; p<D.length; ++p)
    wait[p]=0;
  
  //  srand48(9199298);
  srand(Opt.seed);

  int* pos = new int[D.length];
  for (int pp=0; pp<D.length; ++pp)
    pos[pp] = pp;
  Float previous_squared_error = 1e35; // FLT_MAX?!
  
  for (int epoch=Opt.readEpoch+1; epoch<=Opt.readEpoch + Opt.nEpochs; ++epoch) {
    if (Opt.shuffle)
      shuffle(D, pos);

    M->reset_squared_error();
    Float sofar=0;
    int batch_cnt=0;
    for (int pp=0; pp<D.length; pp++) {
      int p= pos[pp]; // drand48()*D.length;
      wait[p]++;
      if (wait[p]<=0) {
	cout << p << ",";
	continue;
      }

      M->extimation(D.seq[p]);
      sofar += (D.seq[p])->length;

      batch_cnt++;
      if (batch_cnt >= D.length/Opt.batch_blocks && 
	  D.length-pp >= D.length/Opt.batch_blocks) {
	M->maximization(Opt.attenuation, Opt.belief);
	//		  	cout << "ha"<<flush;

	//	cout << -M->get_squared_error() << "\n";
	//	cout << "(" << pp << ")" << flush;
	batch_cnt=0;
      }
      if (pp%20==0) {
	//		  	for (int cy=0;cy<Opt.cycles;cy++) {
	//				cout << sqrt(M->getdcycles()[cy]/sofar) << " ";
	//			}
	//			cout << "\n" << flush;
	//			cout << M->get_squared_error() << "." << flush;
	cout << "." << flush;
      }
    }
    if (batch_cnt>0)
      M->maximization(Opt.attenuation, Opt.belief);

    Float current_squared_error = M->get_squared_error();
    cout << "\nEpoch " << epoch << " Error= " << current_squared_error;

    save(0,M);
    if (current_squared_error < previous_squared_error) {
      gui=0;
      save(0,M);
      if (Gui>0) {
	//      ep += ep0*0.01;
	//      M->setEpsilon(ep);
      }
      previous_squared_error=current_squared_error;
    }

    if (epoch && epoch%5==0) {
      save(epoch, M);
      evaluate(M, D, "train");
      evaluate(M, T, "test");
      D.write("train-predictions");
      T.write("test-predictions");
    }
    cout << "\n"<<flush;
  }
}



int
main(int argc, char** argv)
{
  if (argc<2) {
    cerr << "Usage: " << argv[0] << " option-file\n";
    exit(1);
  }
  ifstream optstream(argv[1]);
  Options Opt(optstream);
  Opt.write();

  Model* M;
  if (Opt.readModel) {
    char tmp[1024];
    sprintf(tmp, "trained-%d.model", Opt.readEpoch);
    cout << "Reading model from " << tmp << "\n";
    ifstream mstream(tmp);
    M = new Model(mstream);
  } else {
    cout << "Creating model\n"<<flush;

    M = new Model(Opt.NF, Opt.NB);

    cout << "Generating random parameters\n"<<flush;
    M->randomize(Opt.seed);
    save(-1, M);
    Opt.readEpoch = 0;
  }

  cout << "Reading train dataset\n"<<flush;
  ifstream dstream("train.dataset");
  DataSet D(dstream);
  //  D.set_belief(Opt.belief);
  cout << "Reading test dataset\n"<<flush;
  ifstream tstream("test.dataset");
  DataSet T(tstream);
  //  T.set_belief(Opt.belief);
  
  train(M, D, T, Opt);

  delete M;
  
  return 0;
}


