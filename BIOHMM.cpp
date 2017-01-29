
// BIOHMM ver. 1.0 (16/9/2005)
// Copyright (C) Gianluca Pollastri 2005


#include "General.h"
#include "BIOHMM.h"
using namespace std;

/*****************
  Public methods
******************/

BIOHMM::BIOHMM(int the_NI, int the_NO, int the_NF, int the_NB): NI(the_NI), NO(the_NO), NF(the_NF), NB(the_NB) {
  alloc();
  resetStats();
  initialiseTables();
}

BIOHMM::BIOHMM(istream& is) {
  is >> NI >> NO >> NF >> NB;
  alloc();
  resetStats();
  readTables(is);
}

BIOHMM::~BIOHMM() {
  delete[] P_F;
  delete[] P_B;
  delete[] Y;

  for (int i=0; i<NI; ++i) {
    for (int f=0; f<NF; ++f) {
      delete[] P_FFI[i][f];
      delete[] P_FFIss[i][f];
    }
    for (int b=0; b<NB; ++b) {
      delete[] P_BBIss[i][b];
      
      for (int f=0; f<NF; ++f) {
	delete[] P_OFBI[i][b][f];
	delete[] P_OFBIss[i][b][f];
      }
    }
  }
}

void BIOHMM::read(istream& is) {
  is >> NI >> NO >> NF >> NB;
  // NOTE: why doesn't it allocate?
  readTables(is);
  resetStats();
}


void BIOHMM::write(ostream& os) {
  os << NI << ' ' << NO << ' ' << NF << ' ' << NB << endl;
  writeTables(os);
}


void BIOHMM::propagate(int length) {
  int t,i,f,f1,b,b1,o;

  for (t=1;t<=length;t++) {
    //marginalise F1FBI down to FBI
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (f1=0;f1<NF;f1++) {
	    FBIa[t][i][b][f] += FFBI[t][i][b][f][f1];
	  }
	  mua[t] += FBIa[t][i][b][f];
	}
      }
    }
    mua[t] = 1.0/mua[t];
    //multiply OFBI by FBI
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (o=0;o<NO;o++) {
	    OFBI[t][i][b][f][o] *= FBIa[t][i][b][f]*mua[t];
	  }
	}
      }
    }

    //marginalise OFBI down to FBI
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (o=0;o<NO;o++) {
	    FBIb[t][i][b][f] += OFBI[t][i][b][f][o];
	  }
	  mub[t] += FBIb[t][i][b][f];
	}
      }
    }
    mub[t] = 1.0/mub[t];

    //multiply BB1FI by FBI
    for (i=0;i<NI;i++) {
      for (f=0;f<NF;f++) {
	for (b=0;b<NB;b++) {
	  for (b1=0;b1<NB;b1++) {
	    BBFI[t][i][f][b][b1] *= FBIb[t][i][b][f]*mub[t];
	  }
	}
      }
    }

    //marginalise BB1FI down to FB1
    for (b1=0;b1<NB;b1++) {
      for (f=0;f<NF;f++) {
	for (i=0;i<NI;i++) {
	  for (b=0;b<NB;b++) {
	    // are we sure this isn't BBFI[t][i][f][b1][b]?
	    FB1[t][b1][f] += BBFI[t][i][f][b][b1];
	  }
	}
	mu1[t] += FB1[t][b1][f];
      }
    }
    mu1[t] = 1.0/mu1[t];

    //multiply F1FBI by FB1
    if (t<length) {
      for (i=0;i<NI;i++) {
	for (b=0;b<NB;b++) {
	  for (f=0;f<NF;f++) {
	    for (f1=0;f1<NF;f1++) {
	      FFBI[t+1][i][b][f][f1] *= FB1[t][b][f1]*mu1[t];  //tricky: check
	    }
	  }
	}
      }
    }
  }

  /*
  // reset separators
  for (i=0;i<NI;i++) {
  for (b=0;b<NB;b++) {
  memset(FBIa[i][b],0,NF*sizeof(Float));
  memset(FBIb[i][b],0,NF*sizeof(Float));
  }
  }
  for (b=0;b<NB;b++) {
  memset(FB1[b],0,NF*sizeof(Float));
  }
  */

  for (t=1;t<=length;t++) {
    //	error -= log(1/mua[t])+log(1/mub[t])+log(1/mu1[t]);
    mua[t]=mub[t]=mu1[t]=0;
  }

  Float temp=0;
  // do it the other way around
  // here it seems it's doing like the known message passing,
  // dividing the new separator potential values by the old
  // perhaps it should be done in the previous passage as well?!
  for (t=1;t<=length;t++) {
    //marginalise BB1FI down to FBI
    for (i=0;i<NI;i++) {
      for (f=0;f<NF;f++) {
	for (b=0;b<NB;b++) {
	  temp=0;
	  for (b1=0;b1<NB;b1++) {
	    temp += BBFI[t][i][f][b][b1];
	  }
	  if (FBIb[t][i][b][f]!=0)
	    FBIb[t][i][b][f] = temp/FBIb[t][i][b][f];
	  mub[t] += FBIb[t][i][b][f];
	}
      }
    }
    mub[t] = 1.0/mub[t];
    //multiply OFBI by FBI
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (o=0;o<NO;o++) {
	    OFBI[t][i][b][f][o] *= FBIb[t][i][b][f]*mub[t];
	    //cout << OFBI[t][i][b][f][o] << ' ';
	  }
	}
      }
    }

    //marginalise OFBI down to FBI
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  temp=0;
	  for (o=0;o<NO;o++) {
	    temp += OFBI[t][i][b][f][o];
	  }
	  if (FBIa[t][i][b][f]!=0)
	    FBIa[t][i][b][f] = temp/FBIa[t][i][b][f];
	  mua[t] += FBIa[t][i][b][f];
	}
      }
    }
    mua[t] = 1.0/mua[t];
    //multiply F1FBI by FBI
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (f1=0;f1<NF;f1++) {
	    FFBI[t][i][b][f][f1] *= FBIa[t][i][b][f]*mua[t];
	  }
	}
      }
    }
    //marginalise F1FBI down to FB1
    for (b=0;b<NB;b++) {
      for (f1=0;f1<NF;f1++) {
	temp=0;
	for (f=0;f<NF;f++) {
	  for (i=0;i<NI;i++) {
	    temp += FFBI[t][i][b][f][f1];
	  }
	}
	//				cout << temp << " " << flush;
	if (FB1[t][b][f1]!=0)
	  FB1[t][b][f1] = temp/FB1[t][b][f1];
	mu1[t] += FB1[t][b][f1];
      }
    }
    mu1[t] = 1.0/mu1[t];
    //multiply BB1FI by FB1
    if (t>1) {
      for (i=0;i<NI;i++) {
	for (f=0;f<NF;f++) {
	  for (b=0;b<NB;b++) {
	    for (b1=0;b1<NB;b1++) {
	      BBFI[t-1][i][f][b][b1] *= FB1[t][b1][f]*mu1[t];
	    }
	  }
	}
      }
    }
    for (t=1;t<=length;t++) {
      //	error -= log(1/mua[t])+log(1/mub[t])+log(1/mu1[t]);
      mua[t]=mub[t]=mu1[t]=0;
    }
  }

}



int max(Float*a,int l) {
  Float m=0;
  int ml=-1;
  for (int g=0;g<l;g++) {
    if (a[g]>m) {
      m=a[g];
      ml=g;
    }
  }
  return ml;
};

void BIOHMM::sufficientStats(int length) {
  int t,i,f,f1,b,b1,o;
  //int imax,bmax,b1max,fmax,f1max,omax;
  
  Float* F;Float* F1;Float* B;Float *B1;Float* O;Float* I;
  F=new Float[NF];
  F1=new Float[NF];
  B=new Float[NB];
  B1=new Float[NB];
  O=new Float[NO];
  I=new Float[NI];

  //resetStats();

  for (t=1;t<=length;t++) {
    memset(F,0,NF*sizeof(Float));
    memset(F1,0,NF*sizeof(Float));
    memset(B,0,NB*sizeof(Float));
    memset(B1,0,NB*sizeof(Float));
    memset(I,0,NI*sizeof(Float));
    memset(O,0,NO*sizeof(Float));
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (o=0;o<NO;o++) {
	    //					if (OFBI[t][i][b][f][o]>max) {
	    //						max = OFBI[t][i][b][f][o];
	    F[f] += OFBI[t][i][b][f][o];
	    B[b] += OFBI[t][i][b][f][o];
	    I[i] += OFBI[t][i][b][f][o];
	    O[o] += OFBI[t][i][b][f][o];
	    //						imax=i;fmax=f;bmax=b;omax=o;
	    //					}
	  }
	}
      }
    }
    //	cout << max << " " << flush;
    // WARNING, POSSIBLE ERROR
    // should second index be max(B,NB) and third max(F,NF)?!
    P_OFBIss[max(I,NI)][max(B,NB)][max(F,NF)][max(O,NO)] += 1.0;
    //	cout << "r"<<flush;
    for (i=0;i<NI;i++) {
      for (f=0;f<NF;f++) {
	for (f1=0;f1<NF;f1++) {
	  for (b=1;b<NB;b++) {
	    F1[f1] += FFBI[t][i][b][f][f1];
	  }
	}
      }
    }
    P_FFIss[max(I,NI)][max(F,NF)][max(F1,NF)] += 1.0;
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (b1=0;b1<NB;b1++) {
	  for (f=1;f<NF;f++) {
	    B1[b1] += BBFI[t][i][f][b][b1];
	  }
	}
      }
    }
    P_BBIss[max(I,NI)][max(B,NB)][max(B1,NB)] += 1.0;
  }

  delete[] F;
  delete[] F1;
  delete[] B;
  delete[] B1;
  delete[] O;
  delete[] I;

}

void BIOHMM::extimation(int *seq, int* y, int length) {
  // initialise tree
  tree_alloc(length);
  //	attach CPTs to all tables
  attachCPT(length);

  // inject evidence
  injectOut(y, length); // multiply P_OFBI by outputs
  injectIn(seq,length); // multiply P_OFBI by inputs

  // propagate evidence
  propagate(length);
  
  // collect and store sufficient stats for all tables
  sufficientStats(length);

  //save the outputs into Y (to check the error)
  saveOutput(length);

  // deallocate tree
  tree_dealloc(length);
}

void BIOHMM::maximization(Float att, Float prior) {
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
  resetStats();
  normaliseTables();

  //printTables();
}

void BIOHMM::printStats() {
  int i,o,f,f1,b,b1;
  Float tot,totF,totB;
  tot=0;totF=0;totB=0;

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (f1=0;f1<NF;f1++) {
	cout << P_FFIss[i][f][f1] << ' ';
      }
    }
    for (b=0;b<NB;b++) {
      for (b1=0;b1<NB;b1++) {
	cout << P_BBIss[i][b][b1] << ' ';
      }
    }
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  cout << P_OFBIss[i][b][f][o] << ' ';
	}
      }
    }
    cout << endl << flush;
  }
}

void BIOHMM::printTables() {
  int i,o,f,f1,b,b1;
  Float tot,totF,totB;
  tot=0;totF=0;totB=0;

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (f1=0;f1<NF;f1++) {
	cout << P_FFI[i][f][f1] << ' ';
      }
    }
    for (b=0;b<NB;b++) {
      for (b1=0;b1<NB;b1++) {
	cout << P_BBI[i][b][b1] << ' ';
      }
    }
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  cout << P_OFBI[i][b][f][o] << ' ';
	}
      }
    }
    cout << endl << flush;
  }
}

void BIOHMM::Feed(int* seq, int length) {
  // initialise tree
  tree_alloc(length);
  attachCPT(length); // attach CPTs to all tables

  // inject evidence
  injectIn(seq,length); // multiply P_OFBI by inputs

  // propagate evidence
  propagate(length);

  //save the outputs into Y
  saveOutput(length);

  // deallocate tree
  tree_dealloc(length);
}


void BIOHMM::predict(int* seq, int length) {
  Feed(seq,length);
}


/*****************
  Private methods
******************/

void BIOHMM::alloc() {
  int i,f,b;

  P_F = new Float[NF];
  P_B = new Float[NB];
  P_OFBI = new Float***[NI];
  P_FFI = new Float**[NI];
  P_BBI = new Float**[NI];
  P_OFBIss = new Float***[NI];
  P_FFIss = new Float**[NI];
  P_BBIss = new Float**[NI];
  Y=new Float[NO*MAX];

  for (i=0;i<NI;i++) {
    P_OFBI[i]=new Float**[NB];
    P_FFI[i]=new Float*[NF];
    P_BBI[i]=new Float*[NB];
    P_OFBIss[i]=new Float**[NB];
    P_FFIss[i]=new Float*[NF];
    P_BBIss[i]=new Float*[NB];
    for (f=0;f<NF;f++) {
      P_FFI[i][f]=new Float[NF];
      P_FFIss[i][f]=new Float[NF];
    }
    for (b=0;b<NB;b++) {
      P_OFBI[i][b]=new Float*[NF];
      P_BBI[i][b]=new Float[NB];
      P_OFBIss[i][b]=new Float*[NF];
      P_BBIss[i][b]=new Float[NB];
      for (f=0;f<NF;f++) {
	P_OFBI[i][b][f]=new Float[NO];
	P_OFBIss[i][b][f]=new Float[NO];
      }
    }
  }
}


void BIOHMM::initialiseTables() {
  int i,o,f,f1,b,b1;

  for (f=0;f<NF;f++) {
    P_F[f]=1.0;
  }
  for (b=0;b<NB;b++) {
    P_B[b]=1.0;
  }

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (f1=0;f1<NF;f1++) {
	P_FFI[i][f][f1]=double(rand());
      }
    }
    for (b=0;b<NB;b++) {
      for (b1=0;b1<NB;b1++) {
	P_BBI[i][b][b1]=double(rand());
      }
    }
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  P_OFBI[i][b][f][o]=double(rand());
	}
      }
    }
  }
  normaliseTables();
}

void BIOHMM::normaliseTables() {
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

void BIOHMM::resetTables() {
  int i,o,f,f1,b,b1;

  for (f=0;f<NF;f++) {
    P_F[f]=1.0;
  }
  for (b=0;b<NB;b++) {
    P_B[b]=1.0;
  }

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (f1=0;f1<NF;f1++) {
	P_FFI[i][f][f1]=0.0;
      }
    }
    for (b=0;b<NB;b++) {
      for (b1=0;b1<NB;b1++) {
	P_BBI[i][b][b1]=0.0;
      }
    }
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  P_OFBI[i][b][f][o]=0.0;
	}
      }
    }
  }
}


void BIOHMM::resetStats() {
  int i,o,f,f1,b,b1;

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (f1=0;f1<NF;f1++) {
	P_FFIss[i][f][f1]=0.0;
      }
    }
    for (b=0;b<NB;b++) {
      for (b1=0;b1<NB;b1++) {
	P_BBIss[i][b][b1]=0.0;
      }
    }
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  P_OFBIss[i][b][f][o]=0.0;
	}
      }
    }
  }
}

void BIOHMM::readTables(istream& is) {
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

void BIOHMM::writeTables(ostream& os) {
  int i,o,f,f1,b,b1;

  for (f=0;f<NF;f++) {
    os << P_F[f] << ' ';
  }
  os << endl;
  
  for (b=0;b<NB;b++) {
    os << P_B[b] << ' ';
  }
  os << endl;

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (f1=0;f1<NF;f1++) {
	os << P_FFI[i][f][f1] << ' ';
      }
    }
    os << endl;
    for (b=0;b<NB;b++) {
      for (b1=0;b1<NB;b1++) {
	os << P_BBI[i][b][b1] << ' ';
      }
    }
    os << endl;
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  os << P_OFBI[i][b][f][o] << ' ';
	}
      }
    }
    os << endl;
  }
}

void BIOHMM::tree_alloc(int length) {
  int t,i,o,f,f1,b,b1;

  mua=new Float[length+1];
  mub=new Float[length+1];
  mu1=new Float[length+1];
  memset(mua,0,(length+1)*sizeof(Float));
  memset(mub,0,(length+1)*sizeof(Float));
  memset(mu1,0,(length+1)*sizeof(Float));

  OFBI = new Float****[length+1];
  FFBI = new Float****[length+1];
  BBFI = new Float****[length+1];

  // junction tree inference initilisation, step 1
  // for each cluster and separator set, sets its corresponding
  // potential to 1
  for (t=0;t<=length;t++) {
    OFBI[t]=new Float***[NI];
    FFBI[t]=new Float***[NI];
    BBFI[t]=new Float***[NI];
    for (i=0;i<NI;i++) {
      OFBI[t][i]=new Float**[NB];
      FFBI[t][i]=new Float**[NB];
      BBFI[t][i]=new Float**[NF];
      for (b=0;b<NB;b++) {
	OFBI[t][i][b]=new Float*[NF];
	FFBI[t][i][b]=new Float*[NF];
	for (f=0;f<NF;f++) {
	  OFBI[t][i][b][f]=new Float[NO];
	  for (o=0;o<NO;o++) {
	    OFBI[t][i][b][f][o]=1.0;
	  }
	  FFBI[t][i][b][f]=new Float[NF];
	  for (f1=0;f1<NF;f1++) {
	    FFBI[t][i][b][f][f1]=1.0;
	  }
	}
      }
      for (f=0;f<NF;f++) {
	BBFI[t][i][f]=new Float*[NB];
	for (b=0;b<NB;b++) {
	  BBFI[t][i][f][b]=new Float[NB];
	  for (b1=0;b1<NB;b1++) {
	    BBFI[t][i][f][b][b1]=1.0;
	  }
	}
      }
    }
  }

  // initialise separators
  // NOTE: I don't understand why it's setting them to 0
  // they should be set to 1 according to step 1 of the
  // initialisation procedure
  FBIa=new Float***[length+1];
  FBIb=new Float***[length+1];
  FB1=new Float**[length+1];
  for (t=0;t<=length;t++) {
    FBIa[t]=new Float**[NI];
    FBIb[t]=new Float**[NI];
    FB1[t]=new Float*[NB];
    for (i=0;i<NI;i++) {
      FBIa[t][i]=new Float*[NB];
      FBIb[t][i]=new Float*[NB];
      for (b=0;b<NB;b++) {
	FBIa[t][i][b]=new Float[NF];
	memset(FBIa[t][i][b],0,NF*sizeof(Float));
	FBIb[t][i][b]=new Float[NF];
	memset(FBIb[t][i][b],0,NF*sizeof(Float));
      }
    }
    for (b=0;b<NB;b++) {
      FB1[t][b] = new Float[NF];
      memset(FB1[t][b],0,NF*sizeof(Float));
    }
  }


} // tree_alloc


void BIOHMM::tree_dealloc(int length) {
  int t,i,f,b;

  delete[] mua;
  delete[] mub;
  delete[] mu1;

  for (t=0;t<=length;t++) {
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  delete[] OFBI[t][i][b][f];
	  delete[] FFBI[t][i][b][f];
	}
	delete[] OFBI[t][i][b];
	delete[] FFBI[t][i][b];
      }
      for (f=0;f<NF;f++) {
	for (b=0;b<NB;b++) {
	  delete[] BBFI[t][i][f][b];
	}
	delete[] BBFI[t][i][f];
      }
      delete[] OFBI[t][i];
      delete[] FFBI[t][i];
      delete[] BBFI[t][i];
    }
    delete[] OFBI[t];
    delete[] FFBI[t];
    delete[] BBFI[t];
  }

  delete[] OFBI;
  delete[] FFBI;
  delete[] BBFI;


  for (t=0;t<=length;t++) {
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	delete[] FBIa[t][i][b];
	delete[] FBIb[t][i][b];
      }
      delete[] FBIa[t][i];
      delete[] FBIb[t][i];
    }
    for (b=0;b<NB;b++) {
      delete[] FB1[t][b];
    }
    delete[] FBIa[t];
    delete[] FBIb[t];
    delete[] FB1[t];
  }
  delete[] FBIa;
  delete[] FBIb;
  delete[] FB1;

} // tree_dealloc


void BIOHMM::attachPOFBI(int t) {
  int i,o,f,b;

  for (i=0;i<NI;i++) {
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  //cout << P_OFBI[i][b][f][o] << " ";
	  OFBI[t][i][b][f][o] *= P_OFBI[i][b][f][o];
	}
      }
    }
  }
  //cout << endl;
} // attachPOFBI

void BIOHMM::attachPFFI(int t) {
  int i,f,f1,b;

  for (i=0;i<NI;i++) {
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (f1=0;f1<NF;f1++) {
	  //cout << P_FFI[i][f][f1] << " ";
	  FFBI[t][i][b][f][f1] *= P_FFI[i][f][f1];
	  //cout << FFBI[t][i][b][f][f1] << " ";
	}
      }
    }
  }
  //cout << endl;
} // attachPFFI

void BIOHMM::attachPBBI(int t) {
  int i,f,b,b1;

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (b=0;b<NB;b++) {
	for (b1=0;b1<NB;b1++) {
	  BBFI[t][i][f][b][b1] *= P_BBI[i][b][b1];
	}
      }
    }
  }
} // attachPBBI

void BIOHMM::attachPF() {
  int i,f,f1,b;

  for (i=0;i<NI;i++) {
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (f1=0;f1<NF;f1++) {
	  FFBI[1][i][b][f][f1] *= P_F[f1];
	}
      }
    }
  }
} // attachPF

void BIOHMM::attachPB(int T) {
  int i,f,b,b1;

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (b=0;b<NB;b++) {
	for (b1=0;b1<NB;b1++) {
	  BBFI[T][i][f][b][b1] *= P_B[b1];
	}
      }
    }
  }
} // attachPB

// this should correspond to the second step of the initialisation
// procedure of the junction tree for inference
void BIOHMM::attachCPT(int length) {
  int t;
  for (t=1;t<=length;t++) {
    attachPOFBI(t);
    attachPFFI(t);
    attachPBBI(t);
  }
  attachPF();
  attachPB(length);
} // attachCPT




void
BIOHMM::injectOut(int* y, int length) {
  int t,i,o,f,b;

  for (t=1;t<=length;t++) {
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (o=0;o<NO;o++) {
	    if (o != y[t])
	      OFBI[t][i][b][f][o] = 0.0;
	  }
	}
      }
    }
  }
} // injectOut

void
BIOHMM::injectIn(int* x, int length) {
  int t,i,o,f,b;

  for (t=1;t<=length;t++) {
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (o=0;o<NO;o++) {
	    if (i != x[t])
	      OFBI[t][i][b][f][o] = 0.0;
	  }
	}
      }
    }
  }
} // injectIn




void BIOHMM::saveOutput(int length) {
  int t,i,o,f,b;
  Float tot=0;

  memset(Y,0,NO*MAX*sizeof(Float));

  for (t=1;t<=length;t++) {
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (o=0;o<NO;o++) {
	    Y[NO*t+o] += OFBI[t][i][b][f][o];
	    //				cout << t << " " << i << " " << b << " " << f << " " << o << " " << OFBI[t][i][b][f][o] << "\n"<<flush;
	  }
	}
      }
    }
    tot=0;
    for (o=0;o<NO;o++) {
      tot += Y[NO*t+o];
      if (Y[NO*t+o]>0) error -= log(Y[NO*t+o]);
    }
    for (o=0;o<NO;o++) {
      Y[NO*t+o]/=tot;
    }
  }
}


