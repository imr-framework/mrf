<p align="center"> <a>
    <img title="Open Source MRF LOGO" src="https://github.com/imr-framework/imr-framework.github.io/blob/master/img/portfolio/mrf.png" width="225">
  </a></p>
<h1 align="center"> Open Source MRF </h1> <br>

**Please note:** You need to install Pulseq and Michigan Image Reconstruction Toolbox (MIRT) in your system for the package to run. Visit [Pulseq](http://pulseq.github.io/) and [MIRT](https://web.eecs.umich.edu/~fessler/code/) for more information.   

This package contains five different parts. To use the package, please pull the complete mrf repository. 

## How to use the package
#### Generate MRF sequence 
1. Call **Spiral_Design/gen_MRF_sequence_pulseq.m**. The function returns parameters of the sequence that are used in other parts of the    package
2. A MRF sequence based on Jiang's paper [1](https://www.ncbi.nlm.nih.gov/pubmed/25491018) is implemented using Pulseq framework.  
3. The user should be able play the generated .seq file on a scanner with Pulseq interpreter installed. This sequence has been verified    on a scanner by the author. 

#### Image Reconstruction
1. Run **Image_Reconstruction/MRF_recon.m**
2. The script lets user select raw data file (.dat) from scanner and trajectory file (.mat) and reconstruct the image using sliding         window method.  

#### Generate EPG dictionary
1. Call **EPG_dict_sim/EPGsim_MRF.m**.
2. A dictionary is generated based on the input parameters range using Extend Phase Graph(EPG)[2](https://www.ncbi.nlm.nih.gov/pubmed/24737382) method.

#### Dictionary matching
1. Call **Dictionary_Matching/MRF_Dic_Matching.m**.
2. The script uses vector dot product to find the best matching dictionary entry and retrieves all parameters for that entry. 

#### ROI analysis tool
1. Call **ROI_analysis_tool/NIST_ROI_ana.m** .
2. The function takes in the parameter maps of ISMRM/NIST phantom and return locations of spheres in correct orders.

## Demo
To test the package, a demo script can be found at **mrf/demo_MRF.m**. This script generates a sequence, then read a sample data,       reconstruct the data. A dictionary is generated and dictionary matching is performed. Then ROI analysis is done to provide T1 and T2 values for each spheres.  

## Reference
[1] Jiang Y, Ma D, Seiberlich N, Gulani V, Griswold M. MR fingerprinting using fast imaging with steady state precession (FISP) with spiral readout. Magn Reson Med. 2015;74(6):spcone-spcone. doi:10.1002/mrm.26048

[2] Weigel M. Extended phase graphs: Dephasing, RF pulses, and echoes - pure and simple. Journal of Magnetic Resonance Imaging. 2014;41(2):266-295. doi:10.1002/jmri.24619



