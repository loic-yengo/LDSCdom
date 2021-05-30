#include "UniSumstat.hpp"

UniSumstat::UniSumstat(string sumstatfile, LDSC ldscores, bool verbose){
  int j;
  int Mi = 0; // starts at 0 becasue there is a title line
  string line = "";
  string tok  = ""; 
  ifstream tmpStream;
  
  tmpStream.open(sumstatfile.c_str());
  while(tmpStream){
    getline(tmpStream,line);
    Mi++;
  }
  tmpStream.close();
  if(verbose){
    cout<<"# Found "<<Mi<<" SNPs in GWAS summary statistics file: "<<sumstatfile<<endl;
  }
  
  tmpStream.open(sumstatfile.c_str());
  int nColScoreFile = 0;
  getline(tmpStream,line);
  stringstream ss;
  ss << line;
  while( ss >> tok )
  {
    ++nColScoreFile;
  }
  if(verbose){
    cout<<"# Found "<<nColScoreFile<<" columns in summary statistics file.\n";
  }
  
  // Format is SNP A1 A2 A1FREQ BETA SE P N
  string z, n;
  string *SNPi   = new string[Mi];
  double *Zi     = new double[Mi];
  double *Ni     = new double[Mi];
  bool*found     = new bool[Mi];
  
  // Read title
  this->M = 0;
  for(j=0;j<Mi;j++){
    tmpStream >> SNPi[j];
    tmpStream >> z;
    tmpStream >> n;
    
    Zi[j]     = -atof(z.c_str()); // ADDED A MINUS SIGN HERE
    Ni[j]     =  atof(n.c_str());
    
    if(ldscores._found[SNPi[j]]){
      // DEAL HERE WITH MISSING DATA
      //if(!isnan(Zi[j]) and !isnan(Ni[j])){
      if(true){  
        if(Ni[j]>0.){
          found[j] = true;
        }
      }
      this->M++;
    }else{
      found[j] = false;
    }
  }
  
  // Allocate object
  this->Z.resize(this->M);    
  this->scL2.resize(this->M);
  this->L2.resize(this->M);  
  this->N.resize(this->M);
  this->block.resize(this->M);
  this->SNP = new string[this->M];
  
  cout<<"# "<<this->M<<" SNPs included in analysis.\n";
  
  int i = -1;
  this->NumBlock = ldscores.NumBlock;
  for(j=0;j<Mi;j++){
    if(found[j]){
      i++;
      this->SNP[i]   = SNPi[j];
      this->N(i)     = Ni[j];
      this->Z(i)     = Zi[j];
      this->block(i) = ldscores._block[SNP[i]];
      this->L2(i)    = ldscores._lds[SNP[i]];
      this->scL2(i)  = this->L2(i) * sqrt( this->N(i) ) / ldscores.M_all;
    }
  }

  // this->b_OLS   = (mxy-mx*my)/(mx2-mx*mx);
  // this->Int_OLS = my - this->b_OLS * mx;
  // this->meanZ   = this->meanZ / this->M;
  // 
  // if(verbose){
  //   cout<<"# Mean of -Z-score  = "<<this->meanZ<<"\n";
  //   cout<<"# OLS ID: Intercept = "<<this->Int_OLS<<" - b = "<<this->b_OLS<<"\n";
  // }
};

void UniSumstat::fitLDSCdom(bool verbose){
  
  int B = this->NumBlock + 1;
  VectorXd Sw   = VectorXd::Zero(B);
  VectorXd Swx  = VectorXd::Zero(B);
  VectorXd Swy  = VectorXd::Zero(B);
  VectorXd Swxy = VectorXd::Zero(B);
  VectorXd Swx2 = VectorXd::Zero(B);
 
  int j,k;
  double w, wx, wx2, wxy, wy;

  for(j=0;j<this->M;j++){
    if(this->L2(j) > 1.){
      w = 1. / L2(j);
    }else{
      w = 1.;
    }

    wx  = w * this->scL2(j);
    wy  = w * this->Z(j);
    wxy = w * this->scL2(j) * this->Z(j);
    wx2 = wx * this->scL2(j);

    Sw(0)   += w;
    Swx(0)  += wx;
    Swy(0)  += wy;
    Swxy(0) += wxy;
    Swx2(0) += wx2;
    
    k        = this->block(j); // get block
    Sw(k)   += w;
    Swx(k)  += wx;
    Swy(k)  += wy;
    Swxy(k) += wxy;
    Swx2(k) += wx2;
  }
  double b_ref = ( (Swxy(0)/Sw(0))  - (Swx(0)/Sw(0)) * (Swy(0)/Sw(0)) ) / ( (Swx2(0)/Sw(0))  - (Swx(0)/Sw(0)) * (Swx(0)/Sw(0)) );
  double i_ref = (Swy(0)/Sw(0)) - b_ref * (Swx(0)/Sw(0));
  double c_ref = Swxy(0) / Swx2(0);

  int Beff = 0; // Effective number of blocks
  double mx, my, mxy, mx2, mw;
  double b_tmp, i_tmp, c_tmp;
  double se_b = 0., se_i = 0., se_c = 0.;
  for(k=1;k<B;k++){
    if(Sw(k)>0.){
      Beff++;
      mxy = Swxy(0) - Swxy(k);
      mx2 = Swx2(0) - Swx2(k);
      mx  = Swx(0)  - Swx(k);
      my  = Swy(0)  - Swy(k);
      mw  = Sw(0)   - Sw(k);
      
      b_tmp = ( (mxy/mw)  - (mx/mw) * (my/mw) ) / ( (mx2/mw)  - (mx/mw) * (mx/mw) );
      i_tmp = (my/mw) - b_tmp * (mx/mw);
      c_tmp = mxy / mx2;
      
      se_b += (b_tmp - b_ref) * (b_tmp - b_ref);
      se_i += (i_tmp - i_ref) * (i_tmp - i_ref);
      se_c += (c_tmp - c_ref) * (c_tmp - c_ref);
    }
  }

  se_b = sqrt( (1. - 1./Beff) * se_b ); 
  se_i = sqrt( (1. - 1./Beff) * se_i ); 
  se_c = sqrt( (1. - 1./Beff) * se_c );
  
  this->b_ID   = b_ref;
  this->i_ID   = i_ref;
  this->c_ID   = c_ref;
  
  if(verbose){
    cout<<"# Effective number of blocks : "<<Beff<<".\n";
    cout<<"# LDSC weighted Estimator : Intercept = "<<this->i_ID<<" ("<<se_i<<") - b = "<<this->b_ID<<" ("<<se_b<<")\n";
    //cout<<"# Intercept-constrained   : b = "<<this->c_ID<<" ("<<se_b<<")\n";
  }
};



