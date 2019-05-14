function Save_kshot_dcf_ind(kshot, dcf, ind)

dim1 = size(kshot);
traj = zeros(dim1(1),dim1(2),2);
traj(:,:,1) = real(kshot);
traj(:,:,2) = imag(kshot);

save ('KTraj_dcf_ind.mat','traj','dcf','ind')
end
