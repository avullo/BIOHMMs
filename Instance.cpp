#include "Instance.h"
using namespace std;

Instance::Instance(istream& is, int quot): alignments_loaded(0), HePl(0) {
  char temp[MAX_T];

  if (quot == 0) is >> name;
  is >> temp;

  length = strlen(temp);

  u = new int[length+1];
  y = new int[length+1];
  y_pred = new int[length+1];
  memset(y_pred, 0, (length+1)*sizeof(int));
  
  app = new Float[3*(length+1)];
  memset(app, 0, 3*(length+1)*sizeof(Float));

  for(int i=1; i<=length; ++i) {
    u[i] = temp[i-1] - 'A';
    if (u[i]<0 || u[i]>26) u[i]=-1;
  }

  char c;
  if(quot == 0)
    for (int i=1; i<=length; ++i) {
      is >> c;
      y[i] = translateY[c-'A'];
      if (y[i]<0 || y[i]>26) y[i]=-1;
    }
}

Instance::~Instance() {
  unload_alignments();
  
  if(HePl)
    unload_profile();

  delete[] app;
  delete[] y_pred;
  delete[] y;
  delete[] u;
}

int Instance::load_alignments(char* alidir, int dir) {
  filebuf inbuf;
  char temp[MAX_T];
  char fname[1024];

  if(dir == 0) {
    strcpy(fname, alidir);
    strncat(fname, name, 5);
    strcat(fname, ".flatblast");
  } else {
    strcpy(fname, alidir);
  }

  if(inbuf.open(fname, ios::in) != 0) {
    istream is(&inbuf);
    is >> alignments_loaded;

    if (alignments_loaded) {
      ALIGNMENTS = new int*[alignments_loaded];
    }

    for (int quanti = 0; quanti<alignments_loaded; ++quanti) {
      is >> temp;

      int ali_length = strlen(temp);
      int maxa = length;

      if(ali_length > maxa) maxa = ali_length;
      ALIGNMENTS[quanti] = new int[maxa+1];
      ALIGNMENTS[quanti][0] = ali_length;
      for(int i=1; i<=ali_length; ++i) {
	if(temp[i-1] == '.' || temp[i-1] == '-')
	  ALIGNMENTS[quanti][i] = 26;
	else
	  ALIGNMENTS[quanti][i] = temp[i-1] - 'A';
	if(ALIGNMENTS[quanti][i]<0 || ALIGNMENTS[quanti][i]>26) {
	  ALIGNMENTS[quanti][i] = 26;
	}
      }
    }
  } else {
    alignments_loaded = 0; 
  }
  
  inbuf.close();

  return alignments_loaded;

}

int Instance::unload_alignments() {
  if(alignments_loaded) {
    for(int quanti=0; quanti<alignments_loaded; ++quanti)
      delete[] ALIGNMENTS[quanti];
    delete[] ALIGNMENTS;
    
    alignments_loaded = 0;
    return 1;
  } else
    return 0;
}

Float Instance::profile_entropy() {
  Float sum=0;
  for(int t=1; t<=length; ++t) {
    for (int aa=0; aa<27; ++aa) {
      if (HeP[27*t+aa])
	sum -= HeP[27*t+aa]*(Float)log(HeP[27*t+aa]);
      if (HeP[27*t+aa]<0 || HeP[27*t+aa]>1) {
	cout << HeP[27*t+aa] << " " << aa << " " << t <<" "<< length<<"\n" << flush;
      }
    }
  }
  
  return sum;
}

int Instance::generate_profile(char* alidir, int dir) {
  if (HePl) return 1;

  Float* Pres = new Float[length+1];
  HeP = new Float[27*(length+1)];

  memset(HeP,0,27*(length+1)*sizeof(Float));
  memset(Pres,0,(length+1)*sizeof(Float));
  
  // Add the main sequence to HeP
  for (int t=1; t<=length; ++t) {
    int aa = u[t];
    if (aa >=0 && aa <27) {
      HeP[27*t+aa]++;
      Pres[t]++;
    }
  }

  load_alignments(alidir, dir);
  int offset = 0;
  int AL = alignments_loaded;

  if(alignments_loaded > 0) {
    for(int a=0; a<AL; ++a) {
      for(int t=1; t<=length; ++t) {
	int aa = ALIGNMENTS[a][t];
	if (aa >=0 && aa <27) {
	  HeP[27*t+aa]++;
	  Pres[t]++;
	}
      }
    }

    for(int t=1; t<=length; ++t) {
      Float total = 0.0;
      Float totalg = 0.0;
      int aa;
      for(int aa=0; aa<26; ++aa) {
	// HeP[27*t+aa] += belief*frequencies[aa];
	total += HeP[27*t+aa];
	totalg += HeP[27*t+aa];
      }
      totalg += HeP[27*t+26];
      for(int aa=0; aa<27; ++aa) {
	if(total)
	  HeP[27*t+aa] /= total;
      }

      HeP[27*t+26] /= totalg;
    }

    Float* weights = new Float[AL];

    //cout << profile_entropy() << "\n" << flush;
    //cout << "\n" << flush;

    for(int iter=0; iter<1; iter++) {
      // We can now think to weight each sequence based on its information content
      memset(weights, 0, AL*sizeof(Float));

      for(int a=0; a<AL; ++a) {
	Float sum = 0;
	for(int t=1; t<=length; ++t) {
	  int aa = ALIGNMENTS[a][t];
	  if(aa >=0 && aa <26) {
	    sum -= (Float)log(HeP[27*t+aa]);
	  }
	}
	sum /= (Float)length;
	weights[a] = sum;
      }

      // Now we recompute the profile matrix
      memset(HeP, 0, 27*(length+1)*sizeof(Float));

      for(int t=1; t<=length; ++t) {
	int aa = u[t];
	if (aa >=0 && aa <27) {
	  HeP[27*t+aa]++;
	}
      }

      for(int a=0; a<AL; ++a) {
	if(1) {
	  for(int t=1; t<=length; ++t) {
	    int aa = ALIGNMENTS[a][t];
	    if (aa >=0 && aa <27) {
	      HeP[27*t+aa] += weights[a];
	    }
	  }
	}
      }

      for(int t=1; t<=length; ++t) {
	Float total = 0.0;
	Float totalg = 0.0;
	int aa;
	for(int aa=0; aa<26; ++aa) {
	  // HeP[27*t+aa] += belief*frequencies[aa];
	  total += HeP[27*t+aa];
	  totalg += HeP[27*t+aa];
	}
	totalg += HeP[27*t+26];
	for(int aa=0; aa<27; ++aa) {
	  if(total)
	    HeP[27*t+aa] /= total;
	}
	HeP[27*t+26] /= totalg;
      }
    }

    delete[] weights;
  }
  
  unload_alignments();
  HePl=1;

  delete[] Pres;

  return 0;

}

inline void Instance::unload_profile() {
  HePl = 0;
  delete[] HeP;
}

void Instance::write(ostream& os) {
  os << name << endl;

  for(int i=1; i<=length; ++i)
    os << (char)('A'+u[i]);
  os << endl;

  for(int i=1; i<=length; ++i)
    os << (char)(Ytranslate[y[i]]);
  os << endl;

  for(int i=1; i<=length; ++i)
    os << (char)(Ytranslate[y_pred[i]]);
  os << endl << endl;
}

void Instance::write_predictions(ostream& os) {
  for(int i=1; i<=length; ++i)
    os << (char)('A'+u[i]);
  os << endl;

  for(int i=1; i<=length; ++i)
    os << (char)(Ytranslate[y_pred[i]]);
  os << endl;

  // the following can be compressed in a cycle of three steps
  for(int t=1; t<=length; ++t) {
    if(y_pred[t] == -1)
      os << "0.0000\t";
    else {
      char num[16];
      sprintf(num, "%.4f", app[3*t]);
      os << num << "\t";
    }
  }
  os << endl;
  
  for(int t=1; t<=length; ++t) {
    if(y_pred[t] == -1)
      os << "0.0000\t";
    else {
      char num[16];
      sprintf(num, "%.4f", app[3*t+1]);
      os << num << "\t";
    }
  }
  os << endl;
  
  for(int t=1; t<=length; ++t) {
    if (y_pred[t] == -1)
      os << "0.0000\t";
    else {
      char num[16];
      sprintf(num, "%.4f", app[3*t+2]);
      os << num << "\t";
    }
  }
  os << endl;
};
