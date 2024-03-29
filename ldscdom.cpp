#include <string.h>
#include <math.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <random>

#include "LDSC.hpp"
#include "UniSumstat.hpp"
using namespace std;

// Main funtion used on cov_files.txt with a sequential access
int main(int argc, char *argv[]){
  if(argc==1){
    cerr<<"# Arguments must be specified. Type --help for more details."<<endl;
    exit(1);
  }
  
  // Read arguments
  string sw = argv[1];
  if (sw == "--help"){
    cerr<<"\t--sumstat    : GWAS summary statistics file with 3 columns. Format: [SNPID] [Dominance Z-score] [N: Sample size]."<<endl;
    cerr<<"\t--ld-score   : LD scores file with 5 columns CHR SNPID BP LDSCORE."<<endl;
    cerr<<"\t--out        : Specifies a prefix for output file: [prefix].ldscdom"<<endl;
    cerr<<"\t--silent     : Different steps of the calculations will NOT be displayed on screen."<<endl;
    exit(1);
  }else{
    if (argc == 1) {
      cerr<<"# Arguments must be specified. Type --help for more details."<<endl;
      exit(1);
    }
  }
  
  // Input parameters
  string sumstatfile = "";
  string ldscfile = "";
  string prefix = "";
  bool verbose = true;
  
  // Indices
  int i;
  
  for(i = 1; i<argc;i++){
    sw = argv[i]; // check if field is not recognised and throw an error.
    if (sw == "--sumstat"){
      sumstatfile = argv[i + 1];
    }
    if (sw == "--ld-score"){
      ldscfile = argv[i + 1];
    }
    if (sw == "--silent"){
      verbose = false;
    } 
    if (sw == "--out"){
      prefix = argv[i + 1];
    }
  }
  
  if(prefix==""){
    cerr<<"\tA prefix for output files must be specified."<<endl;
    exit(1);
  }
  if(sumstatfile==""){
    cerr<<"\tA GWAS summary statistics file must be specified. Use [--sumstat] or check help [--help]."<<endl;
    exit(1);
  }
  if(ldscfile==""){
    cerr<<"\tA LD score file must be specified. Use [--ld-score] or check help [--help]."<<endl;
    exit(1);
  }
  
  if(verbose){
    cout << "#################################################################\n";
    cout << "# ldscdom 0.1 Beta                                              #\n";
    cout << "# LD score regression estimation of Inbreeding Depression       #\n";
    cout << "# Author: Loic Yengo                                            #\n";
    cout << "# License (TBD)                                                 #\n";
    cout << "#################################################################\n";
  }  
  clock_t tic = clock();
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  
  if(verbose){
    cout <<"# Analysis starts: ";
    cout << (now->tm_year + 1900) << '-' 
         << (now->tm_mon + 1) << '-'
         <<  now->tm_mday << " at "
         <<  now->tm_hour <<":"<<now->tm_min<<":"<<now->tm_sec
         <<  endl;
  }
  
  string logFile = prefix+".ldscdom.log";
  ofstream fileLog(logFile.c_str());
  fileLog << "#################################################################\n";
  fileLog << "# ldscdom 0.1 Beta                                              #\n";
  fileLog << "# LD score regression estimation of Inbreeding Depression       #\n";
  fileLog << "# Author: Loic Yengo                                            #\n";
  fileLog << "# License (TBD)                                                 #\n";
  fileLog << "#################################################################\n";
  fileLog <<"# Analysis starts: ";
  fileLog << (now->tm_year + 1900) << '-' 
          << (now->tm_mon + 1) << '-'
          <<  now->tm_mday << " at "
          <<  now->tm_hour <<":"<<now->tm_min<<":"<<now->tm_sec
          << endl;
  
  LDSC ldscores   = LDSC(ldscfile,verbose,fileLog);
  UniSumstat gwas = UniSumstat(sumstatfile,ldscores,verbose,fileLog);
  gwas.fitLDSCdom(verbose,fileLog);

  time_t t2 = time(0);   // get time now
  struct tm * now2 = localtime( & t2 );
  if(verbose){
    cout <<"# Analysis ends: ";
    cout << (now2->tm_year + 1900) << '-' 
         << (now2->tm_mon + 1) << '-'
         <<  now2->tm_mday << " at "
         <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
         << endl;
  }
  clock_t toc = clock();
  
  float timeElapsed = (float)(toc - tic) / CLOCKS_PER_SEC;
  if(verbose){
    printf("# Time elapsed: %f seconds\n", timeElapsed);
  }
  
  fileLog <<"# Analysis ends: ";
  fileLog << (now2->tm_year + 1900) << '-' 
          << (now2->tm_mon + 1) << '-'
          <<  now2->tm_mday << " at "
          <<  now2->tm_hour <<":"<<now2->tm_min<<":"<<now2->tm_sec
          << endl;
  fileLog<<"# Time elapsed: "<<timeElapsed<<" seconds\n";
  fileLog.close();
  return EXIT_SUCCESS;
}


