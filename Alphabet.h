#ifndef ALPHABET_H
#define ALPHABET_H

#include <iostream>
#include <stdexcept>
#include <cstddef>
#include <string>

class Alphabet {
 public:
  virtual int encode(char) = 0;
  virtual ~Alphabet() {}
  
  class BadAlphabetCreation: public std::logic_error {
  public:
  BadAlphabetCreation(std::string type):
    logic_error("Cannot create type " + type) {}
  };
  
  // factory method
  static Alphabet* factory(const std::string&)
    throw(BadAlphabetCreation);
  
};

class AminoAcidSymbols: public Alphabet {
  AminoAcidSymbols() {}
  friend class Alphabet;
  
 public:
  int encode(char);
  ~AminoAcidSymbols() {}

};

class SecondaryStructureSymbols: public Alphabet {
  SecondaryStructureSymbols() {}
  friend class Alphabet;
  
 public:
  int encode(char);  
  ~SecondaryStructureSymbols() {}

};

#endif // ALPHABET_H
