// josie 19 aug 2024                                                                                                                                                                               
// macro to draw each peak, formatted nicely, on 6 panel canvas                                                                                                                                    

void drawone(TString pair, TFile* file){
  TProfile* peak = new TProfile;
  file->GetDirectory(pair)->GetObject("fittedCF",peak);
  peak->SetTitle(Form("%s",pair.Data()));
  peak->RebinX(4);
  peak->GetXaxis()->SetRangeUser(-128,128);
  peak->Draw();
}

void drawPeaksFast(TString filename){

  TFile* file = new TFile(filename.Data(),"READONLY");

  TCanvas* c1 = new TCanvas;
  c1->Divide(3,2);

  c1->cd(1);
  drawone("T1T2",file);
  c1->cd(2);
  drawone("T1T3",file);
  c1->cd(3);
  drawone("T1T4",file);
  c1->cd(4);
  drawone("T2T3",file);
  c1->cd(5);
  drawone("T2T4",file);
  c1->cd(6);
  drawone("T3T4",file);

}
