function filters = morlet_filter_bank_2d_spiral(size_in, options,isGPU,isComplex)
%2d filter bank of Morlet wavelets using a spiral topology, so that they can be mapped
%to uni-dimensional filters

if(nargin<2)
    options = struct;
end
if(nargin<3)
    isGPU=0;
end
if(nargin<4)
    isComplex=0;
end

white_list = {'Q', 'L', 'J', 'sigma_phi', 'sigma_psi', ...
    'xi_psi', 'slant_psi', 'min_margin', 'precision', ...
    'filter_format'};
check_options_white_list(options, white_list);

% Options
options = fill_struct(options, 'Q',1);
options = fill_struct(options, 'L',8);
options = fill_struct(options, 'J',4);
J = options.J;
Q = options.Q;
L = options.L;
options = fill_struct(options, 'sigma_phi',  0.8);
options = fill_struct(options, 'sigma_psi',  0.8);
options = fill_struct(options, 'xi_psi',  1/2*(2^(-1/Q)+1)*pi);
options = fill_struct(options, 'slant_psi',  2/3);%4/L); %TODO re-optimize this quantity now
options = fill_struct(options, 'filter_format', 'fourier_multires');
options = fill_struct(options, 'min_margin', options.sigma_phi * 2^(J/Q) );
options = fill_struct(options, 'precision', 'double');

sigma_phi  = options.sigma_phi;
sigma_psi  = options.sigma_psi;
xi_psi     = options.xi_psi;
slant_psi  = options.slant_psi;

switch options.precision
    case 'single'
        cast = @double;
    case 'double'
        cast = @double;
    otherwise
        error('precision must be either double or single');
end

% Size
res_max = floor(J/Q);
size_filter = pad_size(size_in, options.min_margin, res_max);
phi.filter.type = 'fourier_multires';

% Compute all resolution of the filters
res = 0;

N = size_filter(1);
M = size_filter(2);

% Compute low-pass filters phi
scale = 2^((J-1) / Q - res);
if(isGPU)
    filter_spatial =  gpuArray(single(gabor_2d_period(N, M, sigma_phi*scale, 1, 0, 0)));
else
    filter_spatial =  double(gabor_2d_period(N, M, sigma_phi*scale, 1, 0, 0));
end
phi.filter = cast(real(fft2(filter_spatial)));
phi.meta.J = J;

phi.filter = optimize_filter(phi.filter, 1, options);

littlewood_final = zeros(N, M);
% Compute band-pass filters psi
if(~isComplex)
    angles = (0:(L-1))  * pi / L;
else
    angles = (0:(L-1))  *2* pi / L;
end
p = 1;
for j = 0:J-1
    for theta = 1:numel(angles)
        angle = angles(theta);
        %scale = 2^(j/Q - res);
        scale = 2^(j/Q + (theta-1)/numel(angles)- res);
        filter_spatial = morlet_2d_noDC_period(N, ...
            M,...
            sigma_psi*scale,...
            slant_psi,...
            xi_psi/scale,...
            angle);
        % SLIGHT MODIF
        if(isGPU)
            psi.filter{p} = gpuArray(single(cast(real(fft2(filter_spatial)))));
        else
            psi.filter{p} = double(cast(real(fft2(filter_spatial))));
            %psi.filter{p} = single(cast((fft2(filter_spatial))));
        end
        
        if 0&(theta==6)&(j>J-2)
            keyboard
        end
        littlewood_final = littlewood_final + ...
            abs(realize_filter(psi.filter{p})).^2;
        
        psi.meta.j(p) = j;
        psi.meta.theta(p) = theta;
        psi.meta.center(p)=xi_psi/scale;
        p = p + 1;
    end
end

% Second pass : renormalize psi by max of littlewood paley to have
% an almost unitary operator
% NB : phi must not be renormalized since we want its mean to be 1
K = max(littlewood_final(:));
for p = 1:numel(psi.filter)
    psi.filter{p} = psi.filter{p}/sqrt(K/2);
    if (mod(p,L)==6)&(p>numel(psi.filter)-Q*L)&0
        keyboard
    end
    psi.filter{p} = optimize_filter(psi.filter{p}, 0, options);
end

filters.phi = phi;
filters.psi = psi;

filters.meta.Q = Q;
filters.meta.J = J;
filters.meta.L = L;
filters.meta.sigma_phi = sigma_phi;
filters.meta.sigma_psi = sigma_psi;
filters.meta.xi_psi = xi_psi;
filters.meta.slant_psi = slant_psi;
filters.meta.size_in = size_in;
filters.meta.size_filter = size_filter;
end
