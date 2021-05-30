#ifndef LDSC_H
#define LDSC_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unordered_map>
#include "Eigen/Core"

using namespace std;

class LDSC{
public:
  int M_all;
  int NumBlock;
  unordered_map<string, double> _lds;
  unordered_map<string, int>  _block;
  unordered_map<string, bool> _found;
  string *SNP;

  // Many different constructors
  // will be considered
  LDSC(string ldscfile, bool verbose);
  ~LDSC(){};
};
#endif
