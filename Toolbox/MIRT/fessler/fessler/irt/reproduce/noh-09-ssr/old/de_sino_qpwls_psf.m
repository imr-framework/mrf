%This is a code to determine regularization parameter(gamma)
%for penalty function. The setup is useful for sect and dect with ignoring
%the functional interaction between component images or component sinograms

clear all;
%Load data
load xtrue.mat;
[n.x n.y] = size(xtrue(:,:,1));
n.b = 140; n.a = 128; f.dx = 0.16 * 2;
%mask
mask = sum(xtrue, 3) > 0;
mask = imdilate(mask, strel('disk', 5));
mask = logical(mask);
%Geometry
ig = image_geom('nx', n.x, 'ny', n.y, 'dx', f.dx, 'mask', mask);
sg = sino_geom('par', 'nb', n.b, 'na', n.a, 'dr', f.dx);
%angular spacing for parallel beam: 180 deg / # of angular samples
sino_mask = logical(ones([n.b n.a]));

%Parameters
domain_sel = 1; %0 for image domain / 1 for sinogram domain
gamma1exp = (-10:5);
gamma1 = 2.^(gamma1exp);
gamma2 = gamma1;
W = []; %weight matrix (internal default: identity)

if domain_sel == 0 
    %%%%%%%%%%%%%%%%
    % Image Domain %
    %%%%%%%%%%%%%%%%
    %System matrix
    G = Gtomo2_strip(sg, ig);   
    %Regularization
    R1 = Reg1(ig.mask);  
    C = R1.C; %differencing matrix (+1/-1)
    for i=1:size(gamma1,2)
        [psf0 var(i) fwhm(i)] = qpwls_psf(G, C, gamma1(i), ig.mask, W);
        %printm('approximate stddev = %g', sqrt(var(i)))
        %Contour plot of PSF
        con_vec = [0.99 0.75 0.5 0.25]*max(max(psf0));
        psf0imax = imax(psf0,2);
        psf0imaxx = psf0imax(1)+[-5:5];
        psf0imaxy = psf0imax(2)+[-5:5];
        psf0c = psf0(psf0imaxx,psf0imaxy);
        if (-4 < gamma1exp(i) && gamma1exp(i) < 6)
            subplot(330+i-7);contour(psf0imaxx,psf0imaxy,psf0c',con_vec);
            xlabel('x');ylabel('y');grid on;zoom on;
            tit = sprintf('Contour of PSF with gamma = 2^{%d}',gamma1exp(i));
            title(tit);
        end
    end
    tit1 = sprintf('FWHM vs Beta in the image domain');
    tit2 = sprintf('Approx. STD vs Beta in the image domain');
    
elseif domain_sel == 1 
    %%%%%%%%%%%%%%%%%%%
    % Sinogram Domain %  
    %%%%%%%%%%%%%%%%%%%
    %System matrix
    G = 1;
    %Regularization
    R1 = Reg1(sino_mask, 'offsets', [1]); %radial direction only
    C = R1.C; 
    for i=1:size(gamma1,2)
        [psf0 var(i) fwhm(i)] = qpwls_psf(G, C, gamma1(i), sino_mask, W);
        %printm('approximate stddev = %g', sqrt(var(i)))
        %Contour plot of PSF
        con_vec = [0.99 0.75 0.5 0.25]*max(max(psf0));
        psf0imax = imax(psf0,2);
        psf0imaxx = psf0imax(1)+[-5:5];
        psf0imaxy = psf0imax(2)+[-5:5];
        psf0c = psf0(psf0imaxx,psf0imaxy);
        if (-4 < gamma1exp(i) && gamma1exp(i) < 6) 
            subplot(330+i-7);contour(psf0imaxx,psf0imaxy,psf0c',con_vec);
            xlabel('radius');ylabel('angle');grid on;zoom on;
            tit = sprintf('Contour of PSF with gamma = 2^{%d}',gamma1exp(i));
            title(tit);
        end
    end
    tit1 = sprintf('FWHM vs Beta in the sinogram domain');
    tit2 = sprintf('Approx. STD vs Beta in the sinogram domain');
end

%beta vs FHWM
figure;subplot(121);plot(log2(gamma1),fwhm,'ro-');grid on;
title(tit1);xlabel('log_2 {\beta}');
ylabel('FWHM of PSF [pixels]')
subplot(122);plot(log2(gamma1),sqrt(var),'b-o');grid on;
title(tit2);xlabel('log_2 {\beta}');
ylabel('Approximate Standard Deviation')


%%When the interaction exists
%for i=1:size(gamma1,2)
%    de_G = block_fatrix({G,G});
%    de_C1 = diag_sp(sqrt(gamma1(i))) * C; %gamma is hidden in de_C
%    de_C2 = diag_sp(sqrt(gamma2(i))) * C; 
%    de_C = block_fatrix({de_C1, de_C2}); 
%    [psf0 var fwhm(i)] = qpwls_psf(de_G, de_C, beta, ig.mask, W); %?? in mask
    %printm('approximate stddev = %g', sqrt(var))
%end

% %test for Reg1 differencing matrix
% clear d;
% a = Reg1(logical(ones([4 4])),'offsets',[1]); %only one direction
% %a = Reg1(logical(ones([4 4]))); %for all dierection;
% c = a.C;
% for i=1:16
%     d(:,i) = c(:,i); %Fatrix to Matrix
% end

