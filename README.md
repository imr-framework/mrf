<p align="center"> <a>
    <img title="Open Source MRF LOGO" src="https://github.com/imr-framework/imr-framework.github.io/blob/master/img/portfolio/mrf.png" width="225">
  </a></p>
<h1 align="center"> Open Source MRF </h1> <br>

**Please note:** You need to install Pulseq and Michigan Image Reconstruction Toolbox (MIRT) in your system for the package to run. Visit [Pulseq](http://pulseq.github.io/) and [MIRT](https://web.eecs.umich.edu/~fessler/code/) for more information.   

This package contains five different parts. To use the package, please pull the complete mrf repository. 

## How to use the package
#### Generate MRF sequence 
1. Run **Spiral_Design/gen_MRF_sequence_pulseq.m**. 
2. A MRF sequence based on Jiang's paper [1](https://www.ncbi.nlm.nih.gov/pubmed/25491018) is implemented using Pulseq framework, with flip angles and TR shown below. 
3. The user should be able play the generated .seq file on a scanner with Pulseq interpreter installed. This sequence has been verified on a scanner by the author. 

#### Generate EPG dictionary
1. Run **EPG_simulation/gen_EPG_dict.m**.
2. A dictionary is generated based on the input parameters range using Extend Phase Graph(EPG)[2](https://www.ncbi.nlm.nih.gov/pubmed/24737382) method.

#### Dictionary matching
1. Run **Dictionary_Matching/MRF_Dic_Matching_v2.m**.
2. The script uses vector dot product to find the best matching dictionary entry and retrieves all parameters for that entry. 

#### Image Reconstruction
1. Run **Image_Reconstruction/MRF_Image_Reconstruction_Complex_Chan.m**
2. The script takes in .dat raw data and trajectory file and reconstruct the image using sliding window method. 

#### ROI analysis tool
1. Run **ROI_analysis_tool/NIST_ROI_ana.m** 
2. The function takes in the parameter maps and return locations of spheres in correct order of ISMRM/NIST phantom.
3. To demo the tool, ROI_ana_demo.m can be ran with data included the folder.

## Demo
To test the package, a demo script can be found at **mrf/demo_MRF.m**. 

## Reference
[1] Jiang Y, Ma D, Seiberlich N, Gulani V, Griswold M. MR fingerprinting using fast imaging with steady state precession (FISP) with spiral readout. Magn Reson Med. 2015;74(6):spcone-spcone. doi:10.1002/mrm.26048

[2] Weigel M. Extended phase graphs: Dephasing, RF pulses, and echoes - pure and simple. Journal of Magnetic Resonance Imaging. 2014;41(2):266-295. doi:10.1002/jmri.24619



