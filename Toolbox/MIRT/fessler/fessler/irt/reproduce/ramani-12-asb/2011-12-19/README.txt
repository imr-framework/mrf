This toolbox contains Matlab code that implements the 
alternating direction method of multipliers (ADMM)
algorithm from 

[1] S. Ramani and J. A. Fessler,
"A splitting-based iterative algorithm
for statistical X-ray CT reconstruction,"
IEEE Trans. Medical Imaging, 
doi: 10.1109/TMI.2011.2175233.

This algorithm performs statistical 2-D X-ray CT reconstruction. 
For comparison, we also include code for 
monotone fast iterative shrinkage-thresholding algorithm (MFISTA)
as described in [1] and the following paper:

[2] A. Beck and M. Teboulle,
"Fast gradient-based algorithms for constrained total variation image 
denoising and deblurring problems,"
IEEE Trans. Image Processing, vol. 18, no. 11, 2009.


Our implementation uses objects from J. Fessler's 
image reconstruction toolbox (IRT) that is available at

http://eecs.umich.edu/~fessler/code/index.html

IRT should be added to Matlab's path for proper working of our code.
Please start from "example1.m" which presents a simple simulation
with a 2-D NCAT phantom (also available in IRT).
Our code was tested using Matlab 2009a running on a PC with LINUX. 

Terms of Use:
The user is free to utilize our code for research purposes, 
but is not authorized to redistribute or commercialize 
this product. It would also be courteous to include a citation 
to a published version of [1] in your papers and presentations 
whenever you use our code. University of Michigan and the 
authors cannot be held liable for damages of any kind incurred 
in using this software.

Contact: 
For further information, to send comments / report bugs, 
please email Sathish Ramani (sramani at umich period edu), 
University of Michigan.
