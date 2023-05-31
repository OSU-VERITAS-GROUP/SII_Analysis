/*
Mackenzie Scott -- Summer 2021 based on SimpleReader.cpp by Mike Lisa 

Read and Perform FFT on values from a TDMS file. Prints out results as a PDF
   * Specifically is looking for the 79.5 MHz noise 

Compile with g++ 'root-config --cflags' -o GoodFFT.exe GoodFFT.cpp 'root-config --libs' -Wno-deprecated -lfftw3 
Execute with path/file on the execution line

*/
#include <iostream>
#include <fstream>
#include <cmath>
#include <fftw3.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>

using namespace std;

int main(int argc, char *argv[]){

	//Make sure user puts in a files
	cout << endl;
        cout << "Alright... well. Let's try to open these files \n\n";
        if (argc < 2) {
                cout << "Please include path to a tdms file to read in\n";
                return 1;
        }
	
	//Take in file from terminal and open it
	char *InputFile = argv [1];
        ifstream MyFile;
        MyFile.open(InputFile, std::ifstream::binary);

	//Read in Header and output it in Hex
	
	const int HeaderSize = 4096;
        char Header[HeaderSize];
        unsigned char* HeaderData;
        HeaderData = (unsigned char*)(Header);
        MyFile.read(Header,HeaderSize);

        printf("----------------- Output Header (Hexadecimal) ------------------\n");
        for (int line=0; line<9; line++){
                for (int i=0; i<16; i++){
                        printf("%02X  ",HeaderData[i+16*line]);
                }
                printf("\n");
        }
        printf("----------------------------------------------------------------\n\n");
	

	int nPointsToTransform = 250000;                                       //CHANGE THIS PARAMETER AROUND TO GET BIGGER FFT!!!
	double SecondsPerReading = 0.000000004;
	double MaxForGoodData = (1.0/SecondsPerReading)/2.0;
	double nGoodDataPoints = 0;
	double Freq[nPointsToTransform];
	double sum[nPointsToTransform] = {0};
	double HighestPoint = 0;
	double LocationHP = 0;

	for (int i=0; i < 9000; ++i){                                          //CHANGE NUMBER OF ITERATIONS TO GET BIGGER FFT!!!
		//Read in data after header and output first part in decimal 
        	char RawData[nPointsToTransform];
      		unsigned char* RawDataUnsigned = (unsigned char*)RawData;
        	signed char* RawDataSigned = (signed char*)RawData;
       		MyFile.read(RawData,nPointsToTransform);

		//Perform Fourier Transform
		fftw_complex in[nPointsToTransform];//Data before transfrom 
		fftw_complex out[nPointsToTransform];//Data after transfrom
		fftw_plan plan = fftw_plan_dft_1d(nPointsToTransform, in, out, FFTW_FORWARD, FFTW_ESTIMATE);//Plan to transform 

		//Assign data to our before array  
		for(int i=0; i < nPointsToTransform; i++){ 
			in[i][0] = RawDataSigned[i];
			in[i][1] = 0.0;
		}	
	
		//execute and destroy plan
		fftw_execute(plan);
		fftw_destroy_plan(plan);
	
		double mag[nPointsToTransform];

		for (int i=1; i < nPointsToTransform; ++i) {
			Freq[i] = (1.0/(nPointsToTransform*SecondsPerReading))*i;
			if (Freq[i] > MaxForGoodData){
				nGoodDataPoints = i;
				break;
			}	
        		mag[i] = sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1]);		
			sum[i] += mag[i];
		}	
    
	}

	for(int i=1; i < nGoodDataPoints; ++i){ 
	 	if (Freq[i] > 60000000.0 && Freq[i] < 80000000.0){
                	if(sum[i] > HighestPoint){
                        HighestPoint = sum[i];
                        LocationHP   = Freq[i];
                	}
		}
	}	

	cout << "(" << LocationHP << "," << HighestPoint << ")" <<  endl;
	cout << "Alright now that we have transformed let's put this into root so we can graph\n\n";

	TCanvas *c =new TCanvas("c", "Fast Fourier Transform", 800, 600);             
	TGraph* gr = new TGraph(nGoodDataPoints, Freq, sum);

        gr->Draw("AC*");
        gr->SetTitle("Fast Fourier Transform");
        gr->GetXaxis()->SetTitle("Frequency (Hz)");
        gr->GetYaxis()->SetTitle("Magnitude");
        c->SaveAs("FastFourierTransform.pdf"); 
	MyFile.close();

return 0;
}


