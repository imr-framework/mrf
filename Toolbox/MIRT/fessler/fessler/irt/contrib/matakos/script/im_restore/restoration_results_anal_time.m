load ../results/mat/res_anal_exact_wave_bsnr40_present2.mat
err29 = err;
time29 = time;
xres29 = xres;
isnr29 = isnr;

load ../results/mat/res_anal_exact_wave_bsnr40_present.mat
err(2502:end,:) = [];
isnr(2502:end,:) = [];
time(2502:end,:,:) = [];
err(:,3) = err29(:,1);
time(:,:,3) = time29(:,:,1);
xres(:,:,3) = xres29(:,:,1);
isnr(:,3) = isnr29(:,1);


% err29 = err;
% time29 = time;
% xres29 = xres;
% isnr29 = isnr;
% load ../results/mat/res_anal_exact_tv_bsnr40.mat
% err(:,29) = err29(:,29);
% % time(:,:,29) = time29(:,:,29);
% xres(:,:,29) = xres29(:,:,29);
% isnr(:,29) = isnr29(:,29);

% Experiments for paper timing comparison
% 1) Split-Bregman
% 2) ADMM-P2
% 3) SB-PAD
% 4) SB-MIL

% TV reg xinf params:
% f.ylim = [-120 0];
% f.xlim_iter = [0 1500];
% f.xlim_cputime = [0 1000];
% f.xlim_time = [0 60];
% f.xlim_timedt = [0 60];
% f.exper = 1:4; % Keep for figures (TV reg)
% f.label = '\infty';
% f.title = 'TV';
% f.fname1 = 'fig_conv_anal_tv_bsnr40_present_iter.eps';
% f.fname2 = 'fig_conv_anal_tv_bsnr40_present_cputime.eps';
% f.fname3 = 'fig_conv_anal_tv_bsnr40_present_time.eps';
% f.fname4 = 'fig_conv_anal_tv_bsnr40_present_timedt.eps';
% f.st = 200;

% Wavelet reg xinf params:
f.ylim = [-100 0];
f.xlim_iter = [0 1000];
f.xlim_cputime = [0 500];
f.xlim_time = [0 100];
f.xlim_timedt = [0 100];
f.exper = 1:4; % Keep for figures (TV reg)
f.label = '\infty';
f.title = 'Wavelet';
f.fname1 = 'fig_conv_anal_wave_bsnr40_present_iter.eps';
f.fname2 = 'fig_conv_anal_wave_bsnr40_present_cputime.eps';
f.fname3 = 'fig_conv_anal_wave_bsnr40_present_time.eps';
f.fname4 = 'fig_conv_anal_wave_bsnr40_present_timedt.eps';
f.st = 100;


% f.st = 500;
% f.exper = 1:27; % All
% f.exper = 1:3; % NCG and ISTA
% f.exper = 8:11; % FISTA and MFISTA
% f.exper = [16 18  20 22 ]; % SALSA
% f.exper = [25 26 27]; % ALP1 and ALP2
% f.exper = [2 3 7:9 10:12 19:22];
% f.exper = [2 3 7:9 10 11 16 17 19 21 22]; % Keep for figures (wavelet reg)
% f.exper = [2 3 7:9 10 12 16 18 19 21 22]; % Keep for figures (TV reg)
% f.exper = 19:20;

len_full = length(err);
len_short = length(err(1:f.st:end,1));

symbols = {'g--d','m-o','b-.*','r--s'};

nrmsLim = -50;
[~,idx] = min(abs(err-nrmsLim));

fprintf('Iteration and time to reach %idB NRMS:\n',nrmsLim);
fprintf('Split-Bregman: iter %i - time %4.2f - %4.2f\n',idx(1), time(idx(1),2,2), err(idx(1),2));
fprintf('ADMM-P2: iter %i - time %4.2f - %4.2f\n',idx(2), time(idx(2),2,2), err(idx(2),2));
fprintf('SB-PAD :iter %i - time %4.2f - %4.2f\n',idx(3), time(idx(3),2,2), err(idx(3),2));
fprintf('SB-MIL: iter %i - time %4.2f - %4.2f\n',idx(4), time(idx(4),2,2), err(idx(4),2));

figure(1)
clf
hold on
for ii=f.exper
% 	plot(0:len_full-1,err(:,ii),symb_line{ii},'LineWidth',2,'MarkerSize',6);
	plot(0:f.st:len_full-1,err(1:f.st:end,ii),symbols{ii},'LineWidth',2,'MarkerSize',8)
	set(gca,'FontSize',13)
end
hold off
ylim(f.ylim)
xlim(f.xlim_iter)
xlabel('Iteration')
ylabel(['20 log_{10}(||x - x_{' f.label '}||/||x_{' f.label '}||)'])
title(['Analysis form convergence speed comparison for ' f.title ' regularization'])
legend('SB','ADMM-P2','SB-PAD','SB-MIL');

print('-depsc','-r300','-f1',f.fname1);

figure(2);
clf
hold on
for ii=f.exper
% 	plot(squeeze(time(:,1,ii)),err(:,ii),symb_line{ii},'LineWidth',2,'MarkerSize',6);
	plot(squeeze(time(1:f.st:end,1,ii)),err(1:f.st:end,ii),symbols{ii},'LineWidth',2,'MarkerSize',8)
	set(gca,'FontSize',13)
end
hold off
ylim(f.ylim)
xlim(f.xlim_cputime)
xlabel('CPU time (sec)')
ylabel(['20 log_{10}(||x - x_{' f.label '}||/||x_{' f.label '}||)'])
title(['Analysis form convergence speed comparison for ' f.title ' regularization'])
legend('SB','ADMM-P2','SB-PAD','SB-MIL');


print('-depsc','-r300','-f2',f.fname2);

% printm('%4.0f\n',time(2001,f.exper)/2)


figure(3);
clf
hold on
for ii=f.exper
% 	plot(squeeze(time(:,2,ii)),err(:,ii),symb_line{ii},'LineWidth',2,'MarkerSize',6);
	plot(squeeze(time(1:f.st:end,2,ii)),err(1:f.st:end,ii),symbols{ii},'LineWidth',2,'MarkerSize',8)
	set(gca,'FontSize',13)
end
hold off
ylim(f.ylim)
xlim(f.xlim_time)
xlabel('Time (sec)')
ylabel(['20 log_{10}(||x - x_{' f.label '}||/||x_{' f.label '}||)'])
title(['Analysis form convergence speed comparison for ' f.title ' regularization'])
legend1 = legend('SB','ADMM-P2','SB-PAD','SB-MIL');

set(legend1,'LineWidth',1 ,'FontSize',12);

print('-depsc','-r300','-f3',f.fname3);

figure(4);
clf
hold on
for ii=f.exper
% 	plot(squeeze(time(:,2,ii)),err(:,ii),symb_line{ii},'LineWidth',2,'MarkerSize',6);
	plot(squeeze(time(1:f.st:end,2,ii)),err(1:f.st:end,ii),symbols{ii},'LineWidth',2,'MarkerSize',8)
	set(gca,'FontSize',13)
end
hold off
ylim(f.ylim)
xlim(f.xlim_timedt)
xlabel('Time (sec)')
ylabel(['20 log_{10}(||x - x_{' f.label '}||/||x_{' f.label '}||)'])
title(['Analysis form convergence speed comparison for ' f.title ' regularization'])
legend1 = legend('SB','ADMM-P2','SB-PAD','SB-MIL');

set(legend1,'LineWidth',1 ,'FontSize',12);

print('-depsc','-r300','-f4',f.fname4);

% printm('%4.0f\n',time(2001,f.exper)/2)
