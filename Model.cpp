#include "Model.h"
using namespace std;

Model::Model() {
  NI = NO = NF = NB = -1;
}

Model::Model(int idim, int odim, int fdim, int bdim, const string& type):
  NI(idim), NO(odim), NF(fdim), NB(bdim), ptype(type) {

  try {
    model = Parameterisation::factory(NI, NO, NF, NB, ptype);
  } catch(Parameterisation::BadParameterisationCreation e) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  Y = new Float[NO*MAX];
  
  allocStats();
  resetStats();
}

Model::Model(istream& is) {
  read(is);
}

Model::~Model() {
  delete model;
  delete[] Y;
  
  deallocStats();
}

void Model::read(istream& is) {
  is >> ptype;
  
  try {
    model = Parameterisation::factory(is, ptype);
  } catch(Parameterisation::BadParameterisationCreation e) {
    cerr << e.what() << endl;
    exit(EXIT_FAILURE);
  }

  NI = model->NI; // inputDim();
  NO = model->NO; // outputDim();
  NF = model->NF; // forwardStateDim();
  NB = model->NB; // backwardStateDim();

  Y = new Float[NO*MAX];
  
  allocStats();
  resetStats();  
}

void Model::write(ostream& os) {
  os << ptype << endl;
  model->write(os);
}

void Model::randomize(int seed) {
  //Net->initWeights(seed);
  //NetF->initWeights(seed);
}

void Model::e_step(Instance* instance) {
  int length = instance->length;
  
  tree_alloc(length); // initialise tree
  attachCPT(length); // attach CPTs to all clique potentials

  // inject evidence
  injectOut(instance->y, length); // multiply P_OFBI by outputs
  injectIn(instance->u, length); // multiply P_OFBI by inputs

  // propagate evidence
  propagate(length);
  
  // collect and store sufficient stats for all tables
  sufficientStats(length);

  //save the outputs into Y (to check the error)
  saveOutput(length);

  tree_dealloc(length);
  
}

void Model::m_step(Float att, Float prior) {
  model->maximisation(P_OFBIss, P_FFIss, P_BBIss, att, prior);
  resetStats();
}

void Model::predict(Instance* instance) {
  int length = instance->length;
  
  // initialise tree
  tree_alloc(length);
  attachCPT(length); // attach CPTs to all tables

  // inject evidence
  injectIn(instance->u, length); // multiply P_OFBI by inputs

  // propagate evidence
  propagate(length);
  
  //save the outputs into Y
  saveOutput(length);

  // deallocate tree
  tree_dealloc(length);

  for(int t=1; t<=instance->length; ++t) {
    Float pred = .0;
    int arg = -1;

    for(int c=0; c<NO; ++c) {
      Float out_c = Y[NO*t + c];
      if (out_c > pred) {
	pred = out_c;
	arg = c;
      }
    }
    instance->y_pred[t] = arg;
  }
}

void Model::allocStats() {
  P_FFIss = new Float**[NI];
  P_BBIss = new Float**[NI];
  P_OFBIss = new Float***[NI];
  
  for(int i=0; i<NI; ++i) {    
    P_FFIss[i] = new Float*[NF];
    for(int f=0; f<NF; ++f)
      P_FFIss[i][f] = new Float[NF];

    P_BBIss[i] = new Float*[NB];
    P_OFBIss[i] = new Float**[NB];
    for(int b=0; b<NB; ++b) {
      P_BBIss[i][b] = new Float[NB];

      P_OFBIss[i][b] = new Float*[NF];
      for(int f=0; f<NF; ++f)
	P_OFBIss[i][b][f] = new Float[NO]; 
    }
  }
}

void Model::deallocStats() {
  for (int i=0; i<NI; ++i) {
    for (int f=0; f<NF; ++f)
      delete[] P_FFIss[i][f];
    
    for (int b=0; b<NB; ++b) {
      delete[] P_BBIss[i][b];
      
      for (int f=0; f<NF; ++f)
	delete[] P_OFBIss[i][b][f];

      delete[] P_OFBIss[i][b];
    }
    
    delete[] P_FFIss[i];
    delete[] P_BBIss[i];
    delete[] P_OFBIss[i];
  }

  delete[] P_FFIss;
  delete[] P_BBIss;
  delete[] P_OFBIss;
}

void Model::resetStats() {
  int i,o,f,f1,b,b1;

  for(int i=0; i<NI; ++i) {
    for(int f=0; f<NF; ++f)
      for(int f1=0; f1<NF; ++f1)
	P_FFIss[i][f][f1] = .0;
   
    for(int b=0; b<NB; ++b) {
      for(int b1=0; b1<NB; ++b1)
	P_BBIss[i][b][b1] =.0;

      for(int f=0; f<NF; ++f)
	for(int o=0; o<NO; ++o)
	  P_OFBIss[i][b][f][o] = .0;
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

void Model::sufficientStats(int length) {
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

void Model::tree_alloc(int length) {
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

void Model::tree_dealloc(int length) {
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

// this should correspond to the second step of the initialisation
// procedure of the junction tree for inference
void Model::attachCPT(int length) {
  int t;
  for (t=1;t<=length;t++) {
    attachPOFBI(t);
    attachPFFI(t);
    attachPBBI(t);
  }
  attachPF();
  attachPB(length);
} // attachCPT

void Model::attachPOFBI(int t) {
  int i,o,f,b;

  for (i=0;i<NI;i++) {
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (o=0;o<NO;o++) {
	  //cout << P_OFBI[i][b][f][o] << " ";
	  OFBI[t][i][b][f][o] *= model->p_OFBI(i, b, f, o);
	}
      }
    }
  }
  //cout << endl;
} // attachPOFBI

void Model::attachPFFI(int t) {
  int i,f,f1,b;

  for (i=0;i<NI;i++) {
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (f1=0;f1<NF;f1++) {
	  //cout << P_FFI[i][f][f1] << " ";
	  FFBI[t][i][b][f][f1] *= model->p_FFI(i, f, f1);
	  //cout << FFBI[t][i][b][f][f1] << " ";
	}
      }
    }
  }
  //cout << endl;
} // attachPFFI

void Model::attachPBBI(int t) {
  int i,f,b,b1;

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (b=0;b<NB;b++) {
	for (b1=0;b1<NB;b1++) {
	  BBFI[t][i][f][b][b1] *= model->p_BBI(i, b, b1);
	}
      }
    }
  }
} // attachPBBI

void Model::attachPF() {
  int i,f,f1,b;

  for (i=0;i<NI;i++) {
    for (b=0;b<NB;b++) {
      for (f=0;f<NF;f++) {
	for (f1=0;f1<NF;f1++) {
	  FFBI[1][i][b][f][f1] *= model->P_F[f1];
	}
      }
    }
  }
} // attachPF

void Model::attachPB(int T) {
  int i,f,b,b1;

  for (i=0;i<NI;i++) {
    for (f=0;f<NF;f++) {
      for (b=0;b<NB;b++) {
	for (b1=0;b1<NB;b1++) {
	  BBFI[T][i][f][b][b1] *= model->P_B[b1];
	}
      }
    }
  }
} // attachPB

void Model::injectIn(int* x, int length) {
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

void Model::injectOut(int* y, int length) {
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

void Model::propagate(int length) {
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

void Model::saveOutput(int length) {
  int t,i,o,f,b;
  Float tot=0;

  memset(Y,0,NO*MAX*sizeof(Float));

  for (t=1;t<=length;t++) {
    for (i=0;i<NI;i++) {
      for (b=0;b<NB;b++) {
	for (f=0;f<NF;f++) {
	  for (o=0;o<NO;o++)
	    Y[NO*t+o] += OFBI[t][i][b][f][o];
	}
      }
    }

    tot=0;
    for (o=0;o<NO;o++) {
      tot += Y[NO*t+o];
      if (Y[NO*t+o]>0) error -= log(Y[NO*t+o]);
    }
    for (o=0;o<NO;o++)
      Y[NO*t+o] /= tot;
  }
}
