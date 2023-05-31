# VERITAS SII Analysis

A repository of information regarding the SII Analysis for the VERITAS telescope performed at Ohio State. THere are a number of technical investiations as well as what is needed to analyze a stellar source. 

## Analyzing Stars

There are a number of steps that are needed to convert raw data to visibility measurements. We will outline some of those here. 

### Creating Correlations 

Currently, there are two methods used for creating correlations a FPGA based method and a C++ software method. Both take the raw TDMS files from the telescopes and correlate them together. The result of the FPGA correlations is text file and the result of C++ software is a root object. The FPGA text files must be converted to root objects.  
