#include "LDSC.hpp"
LDSC::LDSC(string ldscfile, bool verbose, ofstream &fileLog){
  int j;
  this->M_all = 0; // starts at 0 becasue there is a title line
  string line = "";
  string tok  = ""; 
  ifstream tmpStream;
  
  tmpStream.open(ldscfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    this->M_all++;
  }
  tmpStream.close();
  if(verbose){
    cout<<"# Found "<<this->M_all<<" SNPs in LD score file: "<<ldscfile<<endl;
  }
  fileLog<<"# Found "<<this->M_all<<" SNPs in LD score file: "<<ldscfile<<endl;
  
  tmpStream.open(ldscfile.c_str());
  int nColScoreFile = 0;
  getline(tmpStream,line);
  stringstream ss;
  ss << line;
  while( ss >> tok )
  {
    ++nColScoreFile;
  }
  if(verbose){
    cout<<"# Found "<<nColScoreFile<<" columns in LD score file.\n";
  }
  fileLog<<"# Found "<<nColScoreFile<<" columns in LD score file.\n";
  
  string snp, lds, block;
  this->SNP = new string[this->M_all];
  // Read title
  int idblock;
  this->NumBlock = 0;
  for(j=0;j<this->M_all;j++){
    tmpStream >> snp;
    tmpStream >> lds;
    tmpStream >> block;
    
    SNP[j] = snp;
    this->_lds.insert({snp,atof(lds.c_str())});
    
    idblock = atoi(block.c_str());
    if(idblock>this->NumBlock){
      this->NumBlock = idblock;
    }
    this->_block.insert({snp,idblock});
    this->_found.insert({snp,true});
  }
  if(verbose){
    cout<<"# Found "<<this->NumBlock<<" blocks in LD score file.\n";
  }
  fileLog<<"# Found "<<this->NumBlock<<" blocks in LD score file.\n";
  // for(j=0;j<10;j++){
  //   cout<<SNP[j]<<"\t"<<_chr[SNP[j]]<<"\t"<<_pos[SNP[j]]<<"\t"<<_lds[SNP[j]]<<endl;
  // }
}
