Matlab code by Antonis Matakos for his image deblurring paper
and MRI ghost correction paper.


The main code files (functions) are in toolbox/im_restore and toolbox/mr_ghost.

toolbox/im_restore contents:

image_restore: function used in paper with ADMM method for analysis formulation

image_restore_synth: Function for ADMM formulation for synthesis formulation.
Not used in paper and not certain if fully complete and/or debugged.
@BlurMatrix: Toolbox used for DCT restoration.
Used in the DCT options of image_restore.
Not essential for the proposed ADMM algorithm.


toolbox/mr_ghost contents:

ghost_correct: Main function used in paper with proposed and other ghost correction methods.


toolbox/mr_trajectories contents:

Function files to create trajectories for the scanner. Also contains modification of mri_trajectory function that uses all the new trajectories as options.


toolbox/systems contents:

block_fatrix: modified to create fatrix with only one block and allow blocks
in any diagonal apart from the main.

Gdiag2: Extension of Gdiag to allow non-square diagonal matrices.
Used as mask for image restoration code
and also for some parts of ghost correction testing.

sparse_bccb: creates a sparse bccb matrix from psf.
Unsure if used for image restoration code.

zero_fatrix: matrix of zeros.
More straightforward implementation than forcing a diagonal with all zeros.


toolbox/utilities:

Some utility functions that may already be includued in the toolbox.
Not sure where they are used.
I think nufft_sinc_nopi is from Michael and was used for the QS algorithm.


script folder:
contains some test scripts for ghost correction and image restoration
along with some trajectory creation code for the EPIs used in the scanner.
Some of the scripts are self-contained but others depend on data files
that are included in the data folder for completeness.
