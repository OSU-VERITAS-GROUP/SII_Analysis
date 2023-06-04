# VERITAS SII Analysis

A repository of information regarding the SII Analysis for the VERITAS telescope performed at Ohio State. There are a number of technical investiations as well as what is needed to analyze a stellar source. 

## Analyzing Stars

There are a number of steps that are needed to convert raw data to visibility measurements. We will outline some of those here. 

### Creating Correlations 

Currently, there are two methods used for creating correlations a FPGA based method and a C++ software method. Both take the raw TDMS files from the telescopes and correlate them together. The result of the FPGA correlations is text file and the result of C++ software is a root object. The FPGA text files must be converted to root objects and we can do this by running ***ConvertNolanFrames.cpp***. 

### Analyzing Correlations

The result of the both correlation programs is a 2D object which is a series (in time) of correlation functions (in relative time). We must "analyze" and time average these series of correlations to obtain an HBT peak (psudo-visibility). The main steps to the analysis are 

1. Normalize each correlation function
2. Correct for ambient signal with *off* data
3. Remove extraneous noise
4. Optical Path Delay Shift 
5. Time Average

The result of this program is ***Analysis.root*** with root objects that can be open in drawn in root.  

### Combining Visibilities 
