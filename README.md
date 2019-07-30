<p align="center"> <a>
    <img title="Open Source MRF LOGO" src="https://github.com/imr-framework/imr-framework.github.io/blob/master/img/portfolio/mrf.png" width="225">
  </a></p>
<h1 align="center"> Open Source MRF </h1> <br>

**Please note:** You need to install Pulseq in your system for the package to run. Visit their [GitHub page](http://pulseq.github.io/) for more information

This package contains several different parts. 

#### Generate MRF sequence 
In MRF_with_Original_FA_TR.m, a MRF sequene based on Jiang's paper [1] is implemented. This sequence has been verified on the scanner. User should be able to play with the sequence or implement their own sequence using the Pulseq framework. 

#### Dictionary matching
In MRF_Dic_Matching_v2.m, a dictionary matching method using vector dot product is implemented. User should be able to match the the dictionary with their signal evolutions and the most similar entry will be returned. 

#### Image Reconstruction
To run MRF_Image_Reconstruction_Complex_Chan.m, Michigan Image Reconstruction Toolbox (MIRT) is required. The script takes in .dat raw data and trajectory file and reconstruct the image using sliding window method. The user should be able to reconstruct their own raw data file and trajectory. 

#### ROI analysis tool
NIST_ROI_ana.m is the main function which takes in the parameter maps and return locations of spheres in correct order of ISMRM/NIST phantom. To test the tool, ROI_ana_demo.m can be ran with data included the folder. 


## Reference
[1] Jiang Y, Ma D, Seiberlich N, Gulani V, Griswold M. MR fingerprinting using fast imaging with steady state precession (FISP) with spiral readout. Magn Reson Med. 2015;74(6):spcone-spcone. doi:10.1002/mrm.26048





