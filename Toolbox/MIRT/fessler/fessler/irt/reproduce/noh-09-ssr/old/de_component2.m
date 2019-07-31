% de_component2
% design 2 component materials for DE decomposition using "PCA" of MAC

if ~isvar('mac')
	kev = 40:160;
kev = 30:160;
	%elem = 1:82; % up to lead
	elem = 1:57; % all elements with k-edges below 40
	[mac kev mtype file] = xray_read_atten(elem, kev);
	amass = xray_read_atomic_mass(elem);
	if 1
		for ii=1:length(file)
			file{ii} = regexprep(file{ii}, '.*\/', '');
	%		mac(:,ii) = mac(:,ii) / max(mac(:,ii)); % normalize
		end
	end
end

if 1
	plot(kev, mac)
	xlabel 'X-ray Energy [keV]'
	ylabel 'Mass Attenuation Coefficient [cm^2/g]'
	title(sprintf('MAC for elements 1 to %d', max(elem)))
	axis([minmax(kev)' 0 20])
%	ir_savefig fig-mac-1-57
prompt
end

if 0 % try macovski formula
	kn = @(a) (1+a) ./ a.^2 .* (2*(1+a)./(1+2*a) - 1 ./ a .* log(1+2*a)) ...
		+ 1/2 ./ a .* log(1+2*a) - (1+3*a) ./ (1+2*a).^2;
	kne = @(kev) kn(kev/510.975);
	fc = @(kev) exp(-0.0028 * (kev - 30)); % * kne(30);
	plot(kev, kne(kev), '-', kev, fc(kev), '--')
return

	iz = 8;
	zz = elem(iz);

	fc = 0.3595 * exp(-0.0028 * (kev - 30));
	fp = 5.9017 * zz^3.8 ./ kev.^3.2;
	tmp = zz / amass(iz) * (1*fp + 1*fc);

	plot(kev, mac(:,iz), '-', kev, tmp, '--')
return
end

%tmp = mac * mac';
%im(kev, kev, tmp)
%[u s v] = svd(tmp);
[u s v] = svd(mac);

if 1 % fix signs
	for ii=1:2
		if all(u(:,ii) < 0)
			u(:,ii) = -u(:,ii);
			v(:,ii) = -v(:,ii);
		end
	end
	max_percent_diff(mac, u * s * v')
end

s = diag(s);
semilogy(s, '-o')
pr sum(s(1:2).^2) / sum(s.^2) * 100

u2 = u(:,1:2);
v2 = v(:,1:2);
if 1
	plot(kev, u2(:,1), '-', kev, u2(:,2:end), '--')
	xlabel 'keV'
	legend('Component 1', 'Component 2')
%	ir_savefig fig-mac-comp-1-2
prompt
end

if 0
	tmp = u2 * diag(s(1:2)) * v2';
	im([tmp- mac])
return
end

if 0
	tmp = mac' * u2;
	tmp(:,1) = tmp(:,1) / max(tmp(:,1));
	tmp(:,2) = tmp(:,2) / max(tmp(:,2));
	plot(tmp)
end
