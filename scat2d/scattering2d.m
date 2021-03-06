function [S, U] = scattering2d(in, filters, options)

local = getoptions(options, 'localized',1);

J = size(filters{1}.psi{1},2);
L = size(filters{1}.psi{1}{1},2);

U{1}{1}.signal = in;
U{1}{1}.asignal = in;
U{1}{1}.scale = 0;
U{1}{1}.orientation=0;

linearindx = @(j,l)(j-1).*L+l;

%zero order x*phi
S{1}{1}.l1 = apply_phi_and_subsample(in,0,filters,J);
S{1}{1}.scale = 0;
S{1}{1}.orientation=0;

%first order |x conv Psi | conv phi
Fin=fft2(in);

for j=1:J
    for l=1:L %to be optimized we can put all the angles together
        r = linearindx(j,l);
        
        S{2}{r}.scale = j;
        S{2}{r}.orientation = l;
        
        % x conv Psi
        aux = ifft2(applyfilter2d(Fin,filters{1}.psi{1}{j}{l}));
        % subsample and abs
        ids = 2^(j-1);
        U{2}{r}.signal = aux(1:ids:end,1:ids:end);
        U{2}{r}.asignal = abs(U{2}{r}.signal);
       
        U{2}{r}.res = log2(ids);
        U{2}{r}.scale = j;
        U{2}{r}.orientation = l;
        
        S{2}{r}.l1 = apply_phi_and_subsample(U{2}{r}.asignal,U{2}{r}.res,filters,J);
    end 
end 
       
%second order  D([||x conv Psi_1 | conv Psi_2 |]_j1)

rast = 1;
for i=1:size(U{2},2) %for every r in U{2}{r}
    j1 = U{2}{i}.scale;
    l1 = U{2}{i}.orientation;
    Fin = fft2(U{2}{i}.asignal);
    res = U{2}{i}.res;%associated to the size of the signals
    
    for j2=j1+1:J % j2> j1
        
        ids = 2^(j2-res-1); %for the subsampling
         
        for l2=1:L
            % |x conv Psi_lambda_1| conv Psi_lambda_2
            aux = ifft2(applyfilter2d(Fin,filters{1}.psi{1}{j2}{l2}));
           
            % subsample and abs
            U{3}{rast}.signal = aux(1:ids:end,1:ids:end);
            U{3}{rast}.asignal = abs(U{3}{rast}.signal);

            U{3}{rast}.res = res+log2(ids);
            U{3}{rast}.scale = [j1 j2];
            U{3}{rast}.orientation = [l1 l2];
        
            %compute the conv with the phi to get S
            S{3}{rast}.l1 = apply_phi_and_subsample(U{3}{rast}.asignal,U{3}{rast}.res,filters,J);

            S{3}{rast}.scale = [j1 j2];
            S{3}{rast}.orientation = [l1 l2];
        
            rast=rast+1;
        end %for l2
    end %for j2
end %for i



end


function l1 = apply_phi_and_subsample(signal,res, filters, J)

res = 1+res;
ids_f = 2^(max(0,floor(J-res)));
 
%apply filter phi
ch = ifft2(applyfilter2d(fft2(signal),filters{1}.phi{1}));
%subsample
l1 = 2^(J-1)*ch(1:ids_f:end,1:ids_f:end);

end


function S=applyfilter2d(signal,filter)

%assuming both in the fourier domain
N  = size(filter);
n = size(signal);

diff_size = (sum(abs(N-n))>0);

if diff_size
    Signal=zeros(size(filter));
    Signal((N(1)/2+(1:n(1)))-floor(n(1)/2),(N(2)/2+(1:n(2)))-floor(n(2)/2))=fftshift(signal);
    Z = Signal.*fftshift(filter);
    S = ifftshift(Z((N(1)/2-n(1)/2)+1:N(1)/2+n(1)/2,(N(2)/2-n(2)/2+1):N(2)/2+n(2)/2));
else 
    S = signal.*filter;
end 

%for optimizing at some point: 
% 
% ch = fft(ch,[], 1);
% ch = ifft(ch .* repmat(filters{1}.phi{res}(:,1), 1 , size(ch,2)));
% 
% ch = ch(1:ids_f:end,:);
% ch = permute(ch, [2 1]);
% ch = fft(ch,[], 1);
% % ch = ifft(ch .* repmat(filters{1}.phi{res}(:,1), 1 , size(ch,2)));
% ch = ifft(ch .* repmat(filters{1}.phi{res}(1,:)', 1 , size(ch,2)));



end