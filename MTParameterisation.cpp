#include "Parameterisation.h"
using namespace std;

MTParameterisation::MTParameterisation() {
  NI = NO = NF = NB = -1;
}

MTParameterisation::MTParameterisation(int idim, int odim, int fdim, int bdim):
  Parameterisation(idim, odim, fdim, bdim) {
  
  alloc();
  initTables();
}

MTParameterisation::MTParameterisation(istream& is) {
  read(is);
}

MTParameterisation::~MTParameterisation() {
  delete[] P_F; delete[] P_B;
  
  for(int i=0; i<NI; ++i) {
    for(int f=0; f<NF; ++f)
      delete[] P_FFI[i][f];

    for(int b=0; b<NB; ++b) {
      delete[] P_BBI[i][b];
      
      for(int f=0; f<NF; ++f)
	delete[] P_OFBI[i][b][f];
      delete[] P_OFBI[i][b];
    }

    delete[] P_FFI[i];
    delete[] P_BBI[i];
    delete[] P_OFBI[i];
  }

  delete[] P_FFI;
  delete[] P_BBI;
  delete[] P_OFBI;
}

void MTParameterisation::read(istream& is) {
  is >> NI >> NO >> NF >> NB;

  alloc();
  readTables(is);  
}

void MTParameterisation::write(ostream& os) {
  os << NI << ' ' << NO << ' ' << NF << ' ' << NB << endl;
  writeTables(os);
}

void MTParameterisation::maximisation(Float**** P_OFBIss, Float*** P_FFIss, Float*** P_BBIss, Float att, Float prior) {
  int i,o,f,f1,b,b1;
  Float tot,totF,totB;
  tot=0;totF=0;totB=0;

  //printTables();

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (f1=0;f1<NF;f1++) {
	totF += P_FFIss[i][f][f1];
      }
    }
    for (b=0;b<NB;b++) {
      for (b1=0;b1<NB;b1++) {
	totB += P_BBIss[i][b][b1];
      }
    }
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  tot += P_OFBIss[i][b][f][o];
	}
      }
    }
  }
  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (f1=0;f1<NF;f1++) {
	P_FFI[i][f][f1] = att*P_FFI[i][f][f1]+(1-att)*P_FFIss[i][f][f1]/totF +prior/double(NI*NF);
      }
    }
    for (b=0;b<NB;b++) {
      for (b1=0;b1<NB;b1++) {
	P_BBI[i][b][b1] = att*P_BBI[i][b][b1]+(1-att)*P_BBIss[i][b][b1]/totB +prior/double(NI*NB);
      }
    }
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  P_OFBI[i][b][f][o] = att*P_OFBI[i][b][f][o]+(1-att)*P_OFBIss[i][b][f][o]/tot +prior/double(NI*NB*NF);
	}
      }
    }
  }

  normaliseTables();
  
}

Float MTParameterisation::p_OFBI(int i, int b, int f, int o) {
  return P_OFBI[i][b][f][o];
}

Float MTParameterisation::p_FFI(int i, int f, int f1) {
  return P_FFI[i][f][f1];
}

Float MTParameterisation::p_BBI(int i, int b, int b1) {
  return P_BBI[i][b][b1];
}


// NOTE: SS initialisation should occur at the model level
void MTParameterisation::alloc() {
  P_F = new Float[NF];
  P_B = new Float[NB];
  
  
  P_FFI = new Float**[NI];
  P_BBI = new Float**[NI];
  P_OFBI = new Float***[NI];
  
  for(int i=0; i<NI; ++i) {
    P_FFI[i] = new Float*[NF];
    
    for(int f=0; f<NF; ++f)
      P_FFI[i][f] = new Float[NF];

    P_BBI[i] = new Float*[NB];
    P_OFBI[i] = new Float**[NB];
    for(int b=0; b<NB; ++b) {
      P_BBI[i][b] = new Float[NB];

      P_OFBI[i][b] = new Float*[NF];
      for(int f=0; f<NF; ++f)
	P_OFBI[i][b][f] = new Float[NO];
    }
  }

}

void MTParameterisation::initTables() {
  for(int f=0; f<NF; ++f)
    P_F[f] = 1.0;
  
  for(int b=0; b<NB; ++b)
    P_B[b] = 1.0;

  for(int i=0; i<NI; ++i) {
    for(int f=0; f<NF; ++f)
      for(int f1=0; f1<NF; ++f1)
	P_FFI[i][f][f1] = double(rand());
      
    for(int b=0; b<NB; ++b) {
      for(int b1=0; b1<NB; ++b1)
	P_BBI[i][b][b1] = double(rand());

      for (int f=0; f<NF; ++f) 
	for(int o=0; o<NO; ++o) 
	  P_OFBI[i][b][f][o] = double(rand());
    }	
  }
  
  normaliseTables();
}

void MTParameterisation::normaliseTables() {
  int i,o,f,f1,b,b1;
  Float tot,totF,totB;

  tot=0;
  for (f=0;f<NF;f++) {
    tot+=P_F[f];
  }
  for (f=0;f<NF;f++) {
    P_F[f]/=tot;
  }
  tot=0;
  for (b=0;b<NB;b++) {
    tot += P_B[b];
  }
  for (b=0;b<NB;b++) {
    P_B[b]/=tot;
  }

  tot=0;totF=0;totB=0;
  for (i=0;i<NI;i++) {
    for (f1=0;f1<NF;f1++) {
      totF=0;
      for (f=0;f<NF;f++) {
	totF += P_FFI[i][f][f1];
      }
      for (f=0;f<NF;f++) {
	P_FFI[i][f][f1] /= totF;
      }
    }
    for (b1=0;b1<NB;b1++) {
      totB=0;
      for (b=0;b<NB;b++) {
	totB+=P_BBI[i][b][b1];
      }
      for (b=0;b<NB;b++) {
	P_BBI[i][b][b1] /= totB;
      }
    }
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	tot=0;
	for (o=0;o<NO;o++) {
	  tot+=P_OFBI[i][b][f][o];
	}
	for (o=0;o<NO;o++) {
	  P_OFBI[i][b][f][o] /= tot;
	}
      }
    }
  }

}

void MTParameterisation::readTables(istream& is) {
  int i,o,f,f1,b,b1;

  for (f=0;f<NF;f++) {
    is>> P_F[f];
  }
  for (b=0;b<NB;b++) {
    is>>P_B[b];
  }

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (f1=0;f1<NF;f1++) {
	is>>P_FFI[i][f][f1];
      }
    }
    for (b=0;b<NB;b++) {
      for (b1=0;b1<NB;b1++) {
	is>>P_BBI[i][b][b1];
      }
    }
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  is>>P_OFBI[i][b][f][o];
	}
      }
    }
  }
}

void MTParameterisation::writeTables(ostream& os) {
  for(int f=0; f<NF; ++f)
    os << P_F[f] << ' ';
  os << endl;
  
  for(int b=0; b<NB; ++b)
    os << P_B[b] << ' ';
  os << endl;

  for(int i=0; i<NI; ++i) {
    for(int f=0; f<NF; ++f)
      for (int f1=0; f1<NF; ++f1)
	os << P_FFI[i][f][f1] << ' ';
    os << endl;
    
    for(int b=0; b<NB; ++b)
      for(int b1=0; b1<NB; ++b1)
	os << P_BBI[i][b][b1] << ' ';
    os << endl;
    
    for(int b=0; b<NB; ++b)
      for(int f=0; f<NF; ++f)
	for(int o=0; o<NO; ++o)
	  os << P_OFBI[i][b][f][o] << ' ';
    os << endl; 
  }
}
