// josie 21 feb 2024, 23 june 2024
// combine old FPGA correlated ANALYSIS root files and python info into one root file per run, adding hour angle to py files
// run in directory of the run, above the pair directories
// run as: "root convertOldFiles.C\(\"gam\ cas\"\)"

// ++ how to add in off runs?!?! (usually saved in a separate text file - this is fine I guess to save the whole text file again here, IF I can find it)

void convertOldFiles(TString sourcein){
    
    string source; //sourcein; //= "gam cas";
    source = sourcein.Data();
    //cout << "enter the star name, with underscores instead of spaces" << endl;  // this is great if I am running files by hand but BAD if I am running a bash script
    //cin >> sourcein;
    //source = sourcein.substr(0,sourcein.find("_")) + " " + sourcein.substr(sourcein.find("_")+1,sourcein.length());
    cout << "source is " << source << endl;
    
    // get date, info from directories to name new output file
    string dirname, year, month, filename;
    dirname = gSystem->pwd();
    dirname = dirname.substr(dirname.find("UTC"), dirname.length());
    year = dirname.substr(dirname.find("UTC")+3, 4);
    month = dirname.substr(dirname.find("UTC")+7, 2);
    dirname = dirname.substr(dirname.find("/")+1, dirname.length());
    filename = source + "_y" + year + "m" + month + dirname + "_OldAnalysis.root";
    cout << "combining each pair root file into one new root file named: " << filename << endl;
    TString newfilename = filename;
    TFile* newfile = new TFile(newfilename, "RECREATE");
    newfile->mkdir("Singles");
    
    TString pairlist[6] = {"T1T2","T1T3","T1T4","T2T3","T2T4","T3T4"};
    
    bool ifadc[4] = {false, false, false, false};  // this isn't an ideal solution but I think better than defining new histograms within loop if not needed
    
    // beginning of loop over possible directories to copy data into new file
    for (int i=0; i<6; i++){
        
        // get root file from old directories
        TFile* oldfile = new TFile(pairlist[i]+"/Analysis.root","READONLY");
        
        if(oldfile->IsOpen()){
            
            cout << "found file for " << pairlist[i] << endl;
            newfile->cd();
            gDirectory->mkdir(pairlist[i]);    gDirectory->cd(pairlist[i]);
            
            string pairstring = pairlist[i].Data();
            
            // get existing objects from analysis file
            TTree* pairhead = new TTree;                oldfile->GetObject("Header",pairhead);              newfile->cd(pairlist[i]);      pairhead->CloneTree()->Write("pairheader");
            TProfile2D* cforig = new TProfile2D;        oldfile->GetObject("CFOriginal", cforig);           cforig->Write("originalCF");
            TProfile2D* cfnorm = new TProfile2D;        oldfile->GetObject("CFNormalized", cfnorm);         cfnorm->Write("normalizedCF");
            TH1D* fftbefore = new TH1D;                 oldfile->GetObject("FFTProjection1", fftbefore);    fftbefore->Write("FFTProjection1");
            TProfile2D* cfnoise1 = new TProfile2D;      oldfile->GetObject("CFMinusNoise1", cfnoise1);      cfnoise1->Write("CFMinusNoise1");
            TProfile2D* cfnoise2 = new TProfile2D;      oldfile->GetObject("CFMinusNoise2", cfnoise2);      cfnoise2->Write("CFMinusNoise2");
            TProfile2D* cfnoise3 = new TProfile2D;      oldfile->GetObject("CFMinusNoise3", cfnoise3);      cfnoise3->Write("CFMinusNoise3");
            TH1D* fftafter = new TH1D;                  oldfile->GetObject("FFTProjection2", fftafter);     fftafter->Write("FFTProjection2");
            TProfile2D* heatmap = new TProfile2D;       oldfile->GetObject("HeatMap", heatmap);             heatmap->Write("heatmap");
            TProfile2D* cfshifted = new TProfile2D;     oldfile->GetObject("JustShift", cfshifted);         cfshifted->Write("shiftedCF");
            TProfile* hbtpeak = new TProfile;           oldfile->GetObject("HBTPeakWithFit", hbtpeak);      hbtpeak->Write("fittedCF");
            
            // get which telescope is 1/2 from header info to save ADCs identified correctly in a separate singles directory (parallel structure to versii file format)
            int t1, t2;
            string tel1, tel2;
            pairhead->SetBranchAddress("T1Br",&t1);     pairhead->SetBranchAddress("T2Br",&t2);     pairhead->GetEvent(0);
            tel1 = "T"+to_string(t1);
            tel2 = "T"+to_string(t2);
            
            newfile->cd("Singles");
            if(!ifadc[t1-1]){   TH1D* adcA = new TH1D;      oldfile->GetObject("ADC1", adcA);   adcA->Write(Form("ADC%s", tel1.c_str())); }
            if(!ifadc[t2-1]){   TH1D* adcB = new TH1D;      oldfile->GetObject("ADC2", adcB);   adcB->Write(Form("ADC%s", tel2.c_str())); }
            newfile->cd(pairlist[i]);
            
            if(t1 == 1 || t2 == 1){ ifadc[0] = true; }
            if(t1 == 2 || t2 == 2){ ifadc[1] = true; }
            if(t1 == 3 || t2 == 3){ ifadc[2] = true; }
            if(t1 == 4 || t2 == 4){ ifadc[3] = true; }
            
            // save python info IN root file, like for versii version
            int nframes = cfshifted->GetYaxis()->GetNbins();
            //TNtuple* geom = new TNtuple("PairGeometryInfo","T1T2 geometry","frame:time:u:v:baseline:OPD:hourAngle");
            TNtuple* geom = new TNtuple("PairGeometryInfo",Form("%s geometry",pairstring.c_str()),"frame:time:u:v:baseline:OPD:hourAngle");
            float vals[7], junk;
            
            ifstream pyinfoFile;
            pyinfoFile.open(Form("delays/pyinfo%s.txt", pairstring.c_str()));       if (!pyinfoFile.is_open()){cout << "Problem opening python file!!!\n";  return 0;}
            
            pyinfoFile >> junk >> junk >> junk >> junk >> junk >> junk >> junk; 
            for (int i=0; i<nframes; i++){
                pyinfoFile >> vals[0] >> vals[1] >> vals[2] >> vals[3] >> vals[4] >> vals[5] >> vals[6]; // frame, time, u, v, baseline, opd, hour angle
                geom->Fill(vals);
            }
            pyinfoFile.close();
            geom->Write("PairGeometryInfo");
        
        }
    }
        
} // end of main
