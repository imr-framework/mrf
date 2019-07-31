%% Step 1: Read Image
I = imread('cameraman.tif');
I = double(I);

%% Step 2: Simulate a Blur Kernel
PSF = fspecial('gaussian',7,10); 

%% Step 3: Compute the Eigenvalues of the blur matrix (the operator)
%H = BlurMatrix(PSF, size(I), 'cir');  %circulant boundary condition
H = BlurMatrix(PSF, size(I), 'refl'); %reflection boundary condition
EigConvMat = H.eigblurmatrix; % the eigenvalue of the blur matrix; 


%% Step 4: Compute the blurred image.
Y = H * I; 
imagesc(Y); axis image; colormap gray; axis off;

%% Step 5: Deblur
X = inv(H' * H + 1e-7) * (H' * Y);
imagesc(X); axis image; colormap gray; axis off;
