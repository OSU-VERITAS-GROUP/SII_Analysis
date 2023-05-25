/*
// Here is when #4 moves but without the middle "data point" -> That shelf was not well defined
const int npt=12;
//Double_t theta[npt+1] = {0.0    , 1.0   , 2.0   , 3.0   , 4.0   ,  6.0  , 8.0   , 6.0   , 4.0   , 3.0   , 2.0};// deg
Double_t xipoint[npt]   = {-35.900,-35.700,-35.500,-35.300,-35.100,-34.701,-35.900,-31.979,-34.304,-34.701,-35.100,-35.300};
Double_t xfpoint[npt]   = {-35.700,-35.500,-35.300,-35.100,-34.701,-34.304,-31.979,-35.900,-34.701,-35.100,-35.300,-35.500};
Double_t yipoint[npt]   = {22.768 , 22.766, 22.761, 22.752, 22.740, 22.705, 23.988, 23.367, 22.657, 22.705, 22.740, 22.752};
Double_t yfpoint[npt]   = {22.766 , 22.761, 22.752, 22.740, 22.705, 22.657, 23.367, 23.988, 22.705, 22.740, 22.752, 22.761};
Double_t zipoint[npt]   = {15.03  , 15.03 , 15.03 , 15.03 , 15.03 , 15.03 , 12.917, 12.917, 15.03 , 15.03 , 15.03 , 15.03 };
Double_t zfpoint[npt]   = {15.03  , 15.03 , 15.03 , 15.03 , 15.03 , 15.03 , 12.917, 12.917, 15.03 , 15.03 , 15.03 , 15.03 };
Double_t Dphi[npt]      = { 0.7010, 0.7251, 0.6549, 0.5745, 1.0450, 0.8133, 9.2046,-9.2321,-0.8064,-0.9231,-0.5903,-0.5475};
Double_t dDphi[npt]     = {0.17   , 0.17  , 0.011 , 0.030 , 0.03  , 0.01  , 0.06  , 0.05  , 0.02  , 0.02  , 0.01  , 0.3   };


// Here is when #3 moves
const int npt=14;
//Double_t theta[npt+1] = {0.0   , 1.0   , 2.0   , 3.0   , 4.0   , 6.0   , 8.0  ,  6.0  , 4.0    , 3.0    , 2.0    , 1.0   , 0.0};//deg
Double_t xipoint[npt]   = {29.400, 29.600, 29.800, 30.000, 30.200, 30.599, 29.400, 33.321, 30.996, 30.599, 30.200 , 30.000 , 29.800 , 29.600};
Double_t xfpoint[npt]   = {29.600, 29.800, 30.000, 30.200, 30.599, 30.996, 33.321, 29.400, 30.599, 30.200, 30.000 , 29.800 , 29.600 , 29.400};
Double_t yipoint[npt]   = {71.568, 71.566, 71.561, 71.552, 71.540, 71.505, 72.788, 72.167, 71.456, 71.505, 71.540 , 71.552 , 71.561 , 71.566}; 
Double_t yfpoint[npt]   = {71.566, 71.561, 71.552, 71.540, 71.505, 71.456, 72.167, 72.788, 71.505, 71.540, 71.552 , 71.561 , 71.566 , 71.568};
Double_t zipoint[npt]   = {17.83 , 17.83 , 17.83 , 17.83 , 17.83 , 17.83 , 15.717, 15.717, 17.83 , 17.83 , 17.83  , 17.83  , 17.83  , 17.83 };
Double_t zfpoint[npt]   = {17.83 , 17.83 , 17.83 , 17.83 , 17.83 , 17.83 , 15.717, 15.717, 17.83 , 17.83 , 17.83  , 17.83  , 17.83  , 17.83 };
Double_t Dphi[npt]      = {0.3393, 0.2984, 0.2773, 0.2701, 0.6073, 0.7692, 7.3872,-7.3994,-0.7454,-0.5915,-0.2595 ,-0.2527 ,-0.3127 ,-0.3240}; 
Double_t dDphi[npt]     = {0.016 , 0.013 , 0.009 , 0.011 , 0.010 , 0.006 , 0.5417, 0.0195, 0.010 , 0.013 , 0.011  ,  0.012 , 0.009  , 0.012 }; 
*/

//T1 T2 T3 T4 
const int npt=30;
Double_t xipoint[npt]   = {135.480, 139.401,44.100, 48.021, 29.400, 29.600, 29.800, 30.000, 30.200, 30.599, 29.400, 33.321, 30.996, 30.599, 30.200 , 30.000 , 29.800 , 29.600, -35.900,-35.700,-35.500,-35.300,-35.100,-34.701,-35.900,-31.979,-34.304,-34.701,-35.100,-35.300};
Double_t xfpoint[npt]   = {139.401, 135.480,48.021, 44.100, 29.600, 29.800, 30.000, 30.200, 30.599, 30.996, 33.321, 29.400, 30.599, 30.200, 30.000 , 29.800 , 29.600 , 29.400, -35.700,-35.500,-35.300,-35.100,-34.701,-34.304,-31.979,-35.900,-34.701,-35.100,-35.300,-35.500};
Double_t yipoint[npt]   = {4.078  , 3.457 , -35.012, -35.633,71.568, 71.566, 71.561, 71.552, 71.540, 71.505, 72.788, 72.167, 71.456, 71.505, 71.540 , 71.552 , 71.561 , 71.566, 22.768 , 22.766, 22.761, 22.752, 22.740, 22.705, 23.988, 23.367, 22.657, 22.705, 22.740, 22.752}; 
Double_t yfpoint[npt]   = {3.457  , 4.078, -35.633, -35.012, 71.566, 71.561, 71.552, 71.540, 71.505, 71.456, 72.167, 72.788, 71.505, 71.540, 71.552 , 71.561 , 71.566 , 71.568, 22.766 , 22.761, 22.752, 22.740, 22.705, 22.657, 23.367, 23.988, 22.705, 22.740, 22.752, 22.761};
Double_t zipoint[npt]   = { 18.15, 18.15,10.317, 10.317,17.83 , 17.83 , 17.83 , 17.83 , 17.83 , 17.83 , 15.717, 15.717, 17.83 , 17.83 , 17.83  , 17.83  , 17.83  , 17.83 ,15.03  , 15.03 , 15.03 , 15.03 , 15.03 , 15.03 , 12.917, 12.917, 15.03 , 15.03 , 15.03 , 15.03};
Double_t zfpoint[npt]   = {18.15, 18.15, 10.317, 10.317,17.83 , 17.83 , 17.83 , 17.83 , 17.83 , 17.83 , 15.717, 15.717, 17.83 , 17.83 , 17.83  , 17.83  , 17.83  , 17.83 ,15.03  , 15.03 , 15.03 , 15.03 , 15.03 , 15.03 , 12.917, 12.917, 15.03 , 15.03 , 15.03 , 15.03 };
Double_t Dphi[npt]      = {-12.3129, 12.2977,-4.0229, 4.0331,0.3393, 0.2984, 0.2773, 0.2701, 0.6073, 0.7692, 7.3872, -7.3994, -0.7454, -0.5915, -0.2595 , -0.2527 , -0.3127 ,-0.3240, 0.7010, 0.7251, 0.6549, 0.5745, 1.0450, 0.8133, 9.2046,-9.2321,-0.8064,-0.9231,-0.5903,-0.5475 }; 
Double_t dDphi[npt]     = {0.024 , 0.04 ,0.03, 0.03,0.016 , 0.013 , 0.009 , 0.011 , 0.010 , 0.006 , 0.5417, 0.0195, 0.010 , 0.013 , 0.011  ,  0.012 , 0.009  , 0.012 ,0.17   , 0.17  , 0.011 , 0.030 , 0.03  , 0.01  , 0.06  , 0.05  , 0.02  , 0.02  , 0.01  , 0.3   }; 



//Here is the Theoretical Equation for Spherical Waves
double Theory(Double_t* par, Double_t xf, Double_t xi, Double_t yf, Double_t yi, Double_t zf, Double_t zi){
  double x      = par[0];
  double y      = par[1];
  double z      = par[2];
  double lambda = par[3];
  double value  = 2*TMath::Pi()*(TMath::Sqrt(TMath::Power(xf-x,2.0)+TMath::Power(yf-y,2.0)+TMath::Power(zf-z,2.0)) - TMath::Sqrt(TMath::Power(xi-x,2.0)+TMath::Power(yi-y,2.0)+TMath::Power(zi-z,2.0)))/(lambda);
  value *=-1;  // note this is an unjustified (for now) sign flip!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  return value;
}

//Here is where we calulate the chi squared
void fcn(Int_t &npar, Double_t* gin, Double_t &f, Double_t* par, Int_t iflag){
  double chi2=0; 
  double TwoPi = 2.0*TMath::Pi();
  
  for (int i=0; i<npt; i++){
    double Dphi_TH = Theory(par, xfpoint[i], xipoint[i], yfpoint[i], yipoint[i], zfpoint[i], zipoint[i]);
   
    while(abs(Dphi_TH-TwoPi)>TwoPi){
      if (Dphi_TH<0.0){Dphi_TH+=TwoPi;}
      else {Dphi_TH-=TwoPi;}
    }
    if (Dphi_TH>TMath::Pi()) Dphi_TH -= TwoPi;
    chi2 += pow(Dphi[i] - Dphi_TH,2) / pow(dDphi[i],2);
  }
  f = chi2;
}

void GeometricNoiseAnalysis(){

  TMinuit* gMinuit = new TMinuit(4); //Here is where we define the number of paramters for minuit
  gMinuit->SetFCN(fcn); // Here we declare the function for minuit as fcn
  
  Double_t arglist[10];
  Int_t ierflg=0;

  arglist[0] = 1;
  gMinuit->mnparm(0,"x"     ,  44.0  , 0.05, 0, 0, ierflg);
  gMinuit->mnparm(1,"y"     ,  55.0  , 0.05, 0, 0, ierflg);
  gMinuit->mnparm(2,"z"     ,  8.5   , 0.05, 0, 0, ierflg);
  gMinuit->mnparm(3,"lambda", 1.75765, 0.05, 0, 0, ierflg);
  
  //////Fix a parameter if you like 
  //gMinuit->FixParameter(0);
  //gMinuit->FixParameter(1);
  //gMinuit->FixParameter(2);
  gMinuit->FixParameter(3); //We fix this Parameter as we are fairly certain we know what the phase is

  arglist[0]=500;
  arglist[1]=1.;
  gMinuit->mnexcm("MIGRAD",arglist,2,ierflg);

  // -----------------------------------------
  // From here on, it is just drawing
  double params[4];
  double temp,dtemp;
  for (int i=0;i<4;i++){
    gMinuit->GetParameter(i,temp,dtemp);
    params[i] = temp;
  }

  double pointNumber[npt];
  double pointNumberNudge[npt];
  double Zero[npt];
  double TheoryValue[npt];
  for (int i=0; i< npt; i++){
    pointNumber[i] = (double)i;
    pointNumberNudge[i] = (double)i+0.05;
    Zero[i] = 0.0;
    TheoryValue[i] = Theory(params, xfpoint[i], xipoint[i], yfpoint[i], yipoint[i], zfpoint[i], zipoint[i]);
  }
 
 //-------------------------
 // Root Graph 
  TGraphErrors* tgData = new TGraphErrors(npt,pointNumber,Dphi,Zero,dDphi);
  tgData->SetTitle(";Increment of Relative Phase Points; Relative Phase (Radians)");
  tgData->SetMarkerStyle(20);
  TGraph* tgTheory = new TGraph(npt,pointNumber,TheoryValue);
  tgTheory->SetMarkerColor(2);
  tgTheory->SetMarkerStyle(21);
  tgTheory->SetLineColor(2);
  
  tgData->Draw();
  tgTheory->Draw("same");
  TLine* zeroline = new TLine(0,0,(double)npt,0.0);
  zeroline->SetLineColor(4);
  zeroline->Draw();
}


