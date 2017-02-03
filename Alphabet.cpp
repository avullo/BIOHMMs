#include "Alphabet.h"
using namespace std;

Alphabet* Alphabet::factory(const string& type)
  throw(BadAlphabetCreation) {
  if(type == "AA") { return new AminoAcidSymbols(); }
  if(type == "SS") { return new SecondaryStructureSymbols(); }

  throw BadAlphabetCreation(type);
}

int AminoAcidSymbols::encode(char c) {
  // check character is upper case, convert if necessary
  // check character is in valid range
  // return integer encoding
}

int SecondaryStructureSymbols::encode(char c) {
  // check character is upper case, convert if necessary
  // check character is in valid range
  // return integer encoding  
}
