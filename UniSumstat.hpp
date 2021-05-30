#ifndef UniSumstat_H
#define UniSumstat_H

#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include "Eigen/Core"
#include "LDSC.hpp"

using namespace std;
using namespace Eigen;

class UniSumstat{
public:
  // Data
  int M;            // Number of SNPs
  int NumBlock;     // Number of blocks
  VectorXd Z;       // -Z-score
  VectorXd L2;      // LD scores
  VectorXd scL2;    // Scaled LD scores: scL2 = L2 sqrt(N) / M;
  VectorXd N;
  VectorXi block;
  string  *SNP;

  // Parameters to estimate
  double b_ID;
  double i_ID;
  double c_ID; // slope from constrained regression

  UniSumstat(string sumstatfile, LDSC ldscores, bool verbose);
  ~UniSumstat(){};
  
  void fitLDSCdom(bool verbose);
};
#endif
