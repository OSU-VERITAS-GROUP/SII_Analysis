// josie 26 jan 2024
// draw objects from analysis root file to popup canvases, since pdf from analysis script is too large

void seeRootPlots(TString filename){
    
    // open file from input
    TFile* analysisfile = new TFile(filename.Data(),"READONLY");
    
    // loop over all directories
    TList* keyList = gDirectory->GetListOfKeys();
    for (int ip=0; ip<keyList->GetSize(); ip++){
      TKey* key = (TKey*)keyList->At(ip);
      TString keyName = key->GetName();
      TString className = key->GetClassName();
      if ((keyName.Data()[0]=='T')&&(className=="TDirectoryFile")){
          cout << "drawing pair " << keyName.Data() << endl;
          TString tee1 = keyName[1];    int t1 = atoi(tee1.Data()); 
          TString tee2 = keyName[3];    int t2 = atoi(tee2.Data());
          
          // open canvas
          TCanvas* c1 = new TCanvas("","",1400,800);
          c1->Divide(3,1);
          
          // get objects from root file
          TProfile2D* heatmap = new TProfile2D;     analysisfile->GetDirectory(keyName)->GetObject("heatmap",heatmap); 
          TProfile* fittedcf = new TProfile;        analysisfile->GetDirectory(keyName)->GetObject("fittedCF",fittedcf);
          TNtuple* geometry = new TNtuple;          analysisfile->GetDirectory(keyName)->GetObject("PairGeometryInfo",geometry);
          TH1D* fft1 = new TH1D;                    analysisfile->GetDirectory(keyName)->GetObject("FFTProjection1",fft1);  fft1->SetTitle("FFT before noise removal");
          TH1D* fft2 = new TH1D;                    analysisfile->GetDirectory(keyName)->GetObject("FFTProjection2",fft2);  fft2->SetTitle("FFT after noise removal");
          TH2D* adc1 = new TH2D;                    analysisfile->GetDirectory("Singles")->GetObject(Form("ADCT%s", tee1.Data()), adc1);
          TH2D* adc2 = new TH2D;                    analysisfile->GetDirectory("Singles")->GetObject(Form("ADCT%s", tee2.Data()), adc2);
          TProfile2D* ps1 = new TProfile2D;         analysisfile->GetDirectory("Singles")->GetObject(Form("PowerSpectrumT%s", tee1.Data()), ps1);
          TProfile2D* ps2 = new TProfile2D;         analysisfile->GetDirectory("Singles")->GetObject(Form("PowerSpectrumT%s", tee2.Data()), ps2);
          
          // draw plots on canvas
          c1->cd(1);    gPad->Divide(2,3);
          c1->cd(1);    gPad->cd(1);    adc1->Draw("COLZ");
          c1->cd(1);    gPad->cd(2);    adc2->Draw("COLZ");
          c1->cd(1);    gPad->cd(3);    ps1->Draw("COLZ");
          c1->cd(1);    gPad->cd(4);    ps2->Draw("COLZ");
          c1->cd(1);    gPad->cd(5);    fft1->Draw();
          c1->cd(1);    gPad->cd(6);    fft2->Draw();
          
          c1->cd(2);    gPad->Divide(1,2);
          c1->cd(2);    gPad->cd(1);    heatmap->Draw("COLZ");
          c1->cd(2);    gPad->cd(2);    fittedcf->Draw();
          
          c1->cd(3);    gPad->Divide(1,3);
          c1->cd(3);    gPad->cd(1);    geometry->Draw("baseline:time");
          c1->cd(3);    gPad->cd(2);    geometry->Draw("OPD:time");
          c1->cd(3);    gPad->cd(3);    geometry->Draw("v:u");
      }
    }
    
}
