// josie 19 july 2024
// off peak error dist for the same area only as the fitted area

TH1D* simpeak(double area, double tau, double sigma, int windowsize);
double peakval(double x, double area, double tau, double sigma);

const double pi = TMath::Pi();

void offPeakErrors(TString filename, TString pairname){//offPeakErrorsOneArea(TString filename, TString pairname){

  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);

  // parameters for simulated peaks/fits
  int npsims(50); //31
  double simsigma(4.0), simtau(0.0);
  double taumin(-10), taumax(10), sigmamin(2), sigmamax(7);
  int nwindows(18), windowsize(40); // ??? discuss noise bias grouping

  // define canvases/dist plots
  TFile* file = new TFile(filename.Data(), "READONLY");
  TProfile* cf = new TProfile;     file->GetDirectory(pairname.Data())->GetObject("projectedCF",cf);
  TProfile* cffit = new TProfile;  file->GetDirectory(pairname.Data())->GetObject("fittedCF",cffit);

  TCanvas* cwindows = new TCanvas("","",1600,700); cwindows->Divide(6,3);
  cwindows->Print("offpeakwindows.pdf[");

  //TCanvas* cdists = new TCanvas("","",600,1600); cdists->Divide(1,3);
  //cdists->Print("fitdists.pdf[");

  TCanvas* crms = new TCanvas("","",1200,1100);   crms->Divide(2,2);
  TGraph* RMSvsArea = new TGraph;   RMSvsArea->SetMarkerStyle(20);     RMSvsArea->SetTitle("RMS vs simulated area;peak area input to simulation;RMS of \"good\" peaks only");
  TGraph* nUncut = new TGraph;      nUncut->SetMarkerStyle(20);        nUncut->SetTitle("number of peaks that pass cuts vs simulated area;peak area input to simulation;number of good peaks");
  TGraph* avgAreas = new TGraph;    avgAreas->SetMarkerStyle(20);      avgAreas->SetTitle("avg area returned vs input area;peak area input to simulation;avg peak area from fits");
  TGraph* tauRMS = new TGraph;      tauRMS->SetMarkerStyle(20);        tauRMS->SetTitle("RMS of tau parameter - input tau;peak area input to simulation;RMS of tau parameter error");
  
  // standard things to draw on plots
  TLine* left = new TLine(-10,-0.6e-6,-10,0.6e-6);
  TLine* right = new TLine(10,-0.6e-6,10,0.6e-6);

  TF1* hbtFit = new TF1("hbtFit","([0]*exp(-pow((x-[1])/[2],2)/2.0))/([2]*sqrt(2*pi))",-300,300);
  hbtFit->SetParName(0,"A");          hbtFit->SetParameter(0,0.0);
  hbtFit->SetParName(1,"#tau_{0}");   hbtFit->SetParameter(1, simtau);        hbtFit->SetParLimits(1, taumin, taumax);
  hbtFit->SetParName(2,"#sigma");     hbtFit->SetParameter(2, simsigma);      hbtFit->SetParLimits(2, sigmamin, sigmamax);

  double simarea = cffit->GetFunction("hbtfit")->GetParameter(0);

  TCanvas* cpdist = new TCanvas("","",1600,600);    cpdist->Divide(2,1);
  cpdist->cd(1);   cffit->GetXaxis()->SetRangeUser(-128,128);  cffit->Draw();
  TH1D* areaInOutDist = new TH1D("areaInOutDist",Form("simulated peaks for real peak area = %f ns;input area - fitted area;counts", simarea), 50,-5e-6,5e-6);


  // --------------------------------- start of loop over simulated peak areas ------------------------------------------------------
  for(int ip=0; ip<npsims; ip++){

  int startbin(1), endbin(40);

    //TGraphErrors* areadist = new TGraphErrors;   areadist->SetMarkerStyle(20);     areadist->SetTitle("areas;window #;area (ns)");
    //TGraphErrors* widthdist = new TGraphErrors;  widthdist->SetMarkerStyle(20);    widthdist->SetTitle("widths;window #;width (ns)");
    //TGraphErrors* taudist = new TGraphErrors;    taudist->SetMarkerStyle(20);      taudist->SetTitle("taus;window #;tau (ns)");

    //double simarea = ip*0.5e-6;

  double sumRMS(0.0);
  int nRMS(0);
  double sumArea(0.0);
  double sumTauRMS(0.0);

  // loop over making windows from OFF peak region
  for(int i=0; i<nwindows; i++){

    // make the fake peak, at a different tau shift each time
    //TRandom3* tr = new TRandom3;
    //tr->SetSeed(0);
    //double taushift;
    TH1D* fakepeak = simpeak(simarea, simtau, simsigma, windowsize);

    int flag(0);

    cwindows->cd(i+1);

    int thisbin = startbin;
    TH1D* thiswindow = new TH1D("thiswindow",Form("window %d",i+1),40,-20,20);

    for(int ix=1; ix<=windowsize; ix++){
      thiswindow->SetBinContent(ix, cf->GetBinContent(thisbin)+fakepeak->GetBinContent(ix));
      thiswindow->SetBinError(ix, cf->GetBinError(thisbin));
      thisbin++;
    }
    thiswindow->GetYaxis()->SetRangeUser(-0.6e-6,1.2e-6);
    hbtFit->SetParameter(0,0.0);   hbtFit->SetParameter(1,simtau);   hbtFit->SetParameter(2,simsigma);
    thiswindow->Fit("hbtFit");
    thiswindow->Draw();
    left->Draw("SAME");
    right->Draw("SAME");

    //if(abs(hbtFit->GetParameter(1) - taumin) < 0.01 || abs(hbtFit->GetParameter(1) - taumax) < 0.01){ flag = 2; }
    //if(abs(hbtFit->GetParameter(2) - sigmamin) < 0.01 || abs(hbtFit->GetParameter(2) - sigmamax) < 0.01){ flag = 3; }

    if(flag < 1){
      areaInOutDist->Fill(simarea - hbtFit->GetParameter(0));
      sumRMS += pow((hbtFit->GetParameter(0) - simarea), 2.0);
      sumArea += hbtFit->GetParameter(0);
      sumTauRMS += pow((hbtFit->GetParameter(1) - simtau), 2.0);  // how to pass in the true random peak tau?
      nRMS++;
    }

    // fill fit parameter distributions
    //areadist->AddPoint(i+1,hbtFit->GetParameter(0));    areadist->SetPointError(i,hbtFit->GetParError(0));
    //taudist->AddPoint(i+1,hbtFit->GetParameter(1));     taudist->SetPointError(i,hbtFit->GetParError(1));
    //widthdist->AddPoint(i+1,hbtFit->GetParameter(2));   widthdist->SetPointError(i,hbtFit->GetParError(2));

    if(i!=nwindows/2){ startbin = thisbin - 20; } // startbin +20??
    else{ startbin = thisbin + 90; }// how big do we want to skip for peak region?? 
    
  } // end of loop over windows

  double rms = sqrt(sumRMS/double(nRMS));
  double avgarea = sumArea/double(nRMS);
  RMSvsArea->AddPoint(simarea, rms);
  nUncut->AddPoint(simarea, nRMS);
  avgAreas->AddPoint(simarea, avgarea);

  cwindows->Print("offpeakwindows.pdf");
  //cdists->cd(1);   areadist->Draw("AP");
  //cdists->cd(2);   taudist->Draw("AP");
  //cdists->cd(3);   widthdist->Draw("AP");
  //cdists->Print("fitdists.pdf");


    } // end of loop over sim peaks

    cwindows->Print("offpeakwindows.pdf]");
    //cdists->Print("fitdists.pdf]");

  crms->cd(1);   RMSvsArea->Draw("AP");
  crms->cd(2);   nUncut->Draw("AP");
  crms->cd(3);   avgAreas->SetMinimum(0);    avgAreas->Draw("AP");
  TLine* one = new TLine(0,0,1,1);    one->Draw("SAME");

  crms->cd(4);
  areaInOutDist->Draw();
  cpdist->cd(2);
  areaInOutDist->Draw();

  cout << "this peak area is " << simarea << endl;
  cout << "areas dist, std dev sigma: " << areaInOutDist->GetStdDev() << endl;
  cout << "areas dist, RMS: " << areaInOutDist->GetRMS() << "    err: " << areaInOutDist->GetRMSError() << endl;


  /*cout << "just a check: " << endl;
  TCanvas* temp = new TCanvas;
  TH1D* test = new TH1D("test","test",40,-20,20);
  for(int i=1; i<=40; i++){
    test->SetBinContent(i,peakval(
    }*/
  
} // end of macro

// =================================================================================================================================

// function to return math value at each bin for the simulated gaussian peak
double peakval(double x, double area, double tau, double sigma){
  double val = (area*exp(-pow((x-tau)/sigma,2)/2.0))/(sigma*sqrt(2*pi));
  return val;
}

// function to return simulated peak to add to data
TH1D* simpeak(double area, double tau, double sigma, int windowsize){
  // generate random tau +/- 5 from 0 rel time
  TRandom3* tr = new TRandom3();
  tr->SetSeed(0);
  double taushift = tau - ((tr->Rndm()*10.0) - 5);
  
  TH1D* thepeak = new TH1D("thepeak",Form("peak area %f",area),windowsize,-(windowsize/2),windowsize/2);
  for(int i=1; i<=windowsize; i++){
    thepeak->SetBinContent(i,peakval(i-1-(windowsize/2), area, taushift, sigma)); // this should probably have some assoc errorbar? 
  }
  return thepeak;
}

