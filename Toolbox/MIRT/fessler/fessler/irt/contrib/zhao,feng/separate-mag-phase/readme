This directory contains codes and data for MRI reconstruction
with separate magnitude and phase regularization via Compressed Sensing (CS).
They are based on the paper:
Feng Zhao, Douglas Noll, Jon-Fredrick Nielsen, J A Fessler.
"Separate magnitude and phase regularization via compressed sensing."
IEEE Trans. Med. Imag., 31(9):1713-23, Sep. 2012.
DOI: 10.1109/TMI.2012.2196707

Running this script requires the entire image reconstruction toolbox
by Prof. Jeffrey A. Fessler:
http://www.eecs.umich.edu/~fessler/code/index.html

codes:
mag_phs_cs_sim.m: an example for simulation experiment (MRI thermometry)
mag_phs_cs_real.m: an example for real experiment (MRI velocity mapping)
pcg_bls.m: the phase updating function with regularizer 1
pcg_bls_exp.m: the phase updating function with regularizer 2
pcg_bls_ep.m: the phase updating function with regularizer 3
pcg_bls_exp_ep.m: the phase updating function with regularizer 4

data:
(1) abd_therm_sim.mat: 
(Yoon Chun Kim and Daehyun Yoon provided magnitude images based on the data
from "ISMRM Reconstruction Challenge 2010";
the phase images are from "ISMRM Reconstruction Challenge 2010")

msk1: true low res mask
msk2: loose low res mask
mask1: true high res mask
im2: "true" high res complex image, interpolated from mt and xt
img: im2 with hot spots added in the phase 
samp1: the sampling pattern for the fully sampled low res data
im2_l: "true" low res complex image, reconstructed from fully sampled low res data (samp1) generated from img
xt2: unwrapped phase of im2_l

(2) abd_vl_real.mat:
(Provided by Jon-Frederick Nielsen)

imt26: the complex image reconstructed by taking DFT
 of the fully sampled k-space data, frame 6
imt20: the complex image reconstructed by taking DFT
 of the fully sampled k-space data, reference frame
dx26: the phase that only contains velocity information;
 it is extracted from imt26 and imt20
kdata2: the raw fully sampled k-space data, 11 frames
msk2: loose mask for image domain, frame 6
im20: the reference complex image reconstructed by the proposed method
 with regularizer 2
fc: rad*s/cm, the ratio between phase and the velocity	
