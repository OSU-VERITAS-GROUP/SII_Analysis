# Noise Analysis 

As you may or may not be aware of there is a 79.5 MHz Noise that frequently interferes with observations at the VERITAS site. It appears sinusoidally in the correlation functions and there are a number of things we can do in order to remove this noise. 

    $A*sin(\omega\tau +\phi)$

## Geometric Noise Analysis 
Using the Phase Parameter we can extract important information about the location of our potential beacon. GeometericNoiseAnalysis.C contains callibation data of incramenetenting the telescopes at fitting with a point source noise model to find the location of the noise beacon.  