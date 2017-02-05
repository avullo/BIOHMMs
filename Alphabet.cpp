#include "Alphabet.h"
#include <cctype>
#include <cstdlib>
using namespace std;

string AminoAcidSymbols::alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";

Alphabet* Alphabet::factory(const string& type)
  throw(BadAlphabetCreation) {
  if(type == "AA") { return new AminoAcidSymbols(); }
  if(type == "SS") { return new SecondaryStructureSymbols(); }

  throw BadAlphabetCreation(type);
}

// TODO: could get encoding from string::find returned value
int AminoAcidSymbols::encode(char c) {
  if(!isupper(c)) // check character is upper case, convert if necessary
     c = toupper(c);
     
  // check character is in valid range
  if(alphabet.find(c) == string::npos) {
    cerr << "Character " << c << " in sequence is not allowed" << endl;
    exit(EXIT_FAILURE);
  }
  
  // return integer encoding as the distance wrt
  // to first character in the alphabet
  return c - alphabet[0];
}

char AminoAcidSymbols::decode(int encoding) {
  if(encoding < 0 || encoding > alphabet.size() - 1) {
    cerr << "AA encoding " << encoding << " out of range" << endl;
    exit(EXIT_FAILURE);
  }

  return alphabet[encoding];
}

string SecondaryStructureSymbols::alphabet = "HEC";

// TODO: could get encoding from string::find returned value
int SecondaryStructureSymbols::encode(char c) {
  if(!isupper(c)) // check character is upper case, convert if necessary
    c = toupper(c);
    
  // check character is in valid range
  if(alphabet.find(c) == string::npos) {
    cerr << "Character " << c << " in sequence is not allowed" << endl;
    exit(EXIT_FAILURE);
  }
    
  // return integer encoding
  if(c == 'H') return 0;
  else if(c == 'E') return 1;
  else if(c == 'C') return 2;
  else {
    cerr << "Shouldn't have this: " << c << endl;
    exit(EXIT_FAILURE);
  }
}

char SecondaryStructureSymbols::decode(int encoding) {
  if(encoding < 0 || encoding > alphabet.size() - 1) {
    cerr << "SS encoding " << encoding << " out of range" << endl;
    exit(EXIT_FAILURE);
  }

  return alphabet[encoding];
}
