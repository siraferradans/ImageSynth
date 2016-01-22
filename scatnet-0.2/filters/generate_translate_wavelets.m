function [filters, options, lpal]=generate_translate_wavelets(N, J, L, type,M)

if nargin < 5
    M=N;
end 

options.gpu =0 ;
options.periodinput=0;
options.border=0;%do not remove borders
options.l1renorm = 0; %Sira:?? why should we renormalize the filters?
options.mirror = 0;
options.localized = 1;
options.lr=1e-3;
options.momentum=0.9;
options.l2scatt = 0;
options.multigrid = 1;

switch type

	case 'image'

	options.J1 = J;
	options.L1 = L;
	options.spiral = 0;
    %options.stretchspiral = 1;
    options.coef_stretching= log2(N(1))/log2(M(1)); %M size audio, N size image
	
	options.J2 = J;
	options.L2 = L;	
	options.N = N;
    options.isComplex = 1; %we want the angle to be mapped to 2*pi (not pi)
    
    
	[filters,lpal] = generate_scatt_filters(options);
	options.nscattcoeffs = nsccoeffs(fwdscatt(randn(options.N), filters, options));

	case 'audio'

	options.J1 = J;
	options.L1 = L;
	options.spiral = 0;
	options.J2 = J;
	options.L2 = L;	
	options.N = N;
	options.onedim = 1;
    
	[filters,lpal] = generate_scatt_filters(options);
	options.nscattcoeffs = nsccoeffs(fwdscatt(randn(options.N,1), filters, options));


end