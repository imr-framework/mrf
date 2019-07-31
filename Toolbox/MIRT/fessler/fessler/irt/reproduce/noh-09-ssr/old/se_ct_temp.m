clear all; %SE-CT version
im off;
sl{1} = linspace(0, 50, 26); %For a given lth material,
sl{2} = linspace(0, 30, 31); %makes sinogram of (Nd) size
stype = 'ps1';
xray = xray_read_spectra(stype);
mtype = {'water','bone'}; 
mac = xray_read_atten(mtype, xray.en);
if im
	clf, semilogy(xray.en, mac), legend(mtype{:})
end
sll = ndgrid_jf('mat', sl{:});
fm = se_ftab_fm(sll,mac(:,1),xray.Ide(:,1)); %both sll are needed
fit = de_ftab_fit(sl{1}, fm, 'type', 'exp', 'mtype', mtype)
g = fit.fgrad(sll(:,:,1));
%fit = de_ftab_fit(sl, fm, 'type', 'poly')

if 1
	fit.show_sp(xray.en, xray.sp)
	prompt
	fit.show_fm(sl, fm)
	prompt
	fit.show_err(sl, fm)
end
