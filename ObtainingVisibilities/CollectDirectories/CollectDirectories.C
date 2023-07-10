void CollectDirectories(){
  
  //Runs on command line to pull out all analysis files 
  gSystem->Exec("ls -1 */*/*/Analysis.root | awk -F \"/Anal\" '{print $1}' > Directories.txt");
  ifstream inFile("Directories.txt");
  
  
  ofstream outFile;
  outFile.open("Subroutine.C");
  outFile << "int Setup(int irun, TString& DirName, int& HumanFlag, int& splitNum){\n";
  outFile << "  splitNum = 1;";
  outFile << "  switch (irun) {\n";

  string line;
  int nn=0;
  //  while (inFile.getline(line,256) != EOF){
  while (std::getline(inFile,line)){
    outFile << "  case " << nn << ":\n";
    outFile << "    DirName = \"" << line << "\";\n";
    outFile << "    HumanFlag = 0;\n";
    outFile << "    break;\n";
    nn++;
  }
  outFile <<   "   default:\n";
  outFile <<   "     return 0;\n";
  outFile <<   "     break;\n";
  outFile <<   "   }\n";
  outFile <<   "  return 1;\n";
  outFile <<   "}\n";
  inFile.close();
  outFile.close();
  
}
