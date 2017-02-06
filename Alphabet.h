#ifndef ALPHABET_H
#define ALPHABET_H

#include <iostream>
#include <stdexcept>
#include <cstddef>
#include <string>

// typedef char Symbol;

class Alphabet {
 public:
  
  class BadAlphabetCreation: public std::logic_error {
  public:
  BadAlphabetCreation(std::string type):
    logic_error("Cannot create type " + type) {}
 };

  /* class BadSymbol: public std::logic_error { */
  /* public: */
  /* BadSymbol(char symbol): */
  /*   logic_error("Cannot encode symbol " + symbol) {} */
  /* }; */
 virtual int size() = 0;
 virtual int encode(char) = 0;
 virtual char decode(int) = 0;
 virtual ~Alphabet() {}
 
 // factory method
 static Alphabet* factory(const std::string&)
   throw(BadAlphabetCreation);
};

class AminoAcidSymbols: public Alphabet {
  AminoAcidSymbols() {}
  friend class Alphabet;
  static std::string alphabet;
  
 public:
  int size() { return alphabet.size(); }
  int encode(char);
  char decode(int);
  ~AminoAcidSymbols() {}
  
};

class SecondaryStructureSymbols: public Alphabet {
  SecondaryStructureSymbols() {}
  friend class Alphabet;
  static std::string alphabet;
  
 public:
  int size() { return alphabet.size(); }
  int encode(char);
  char decode(int);
  ~SecondaryStructureSymbols() {}

};

#endif // ALPHABET_H
