function [out,total_cost]=scattering2d_synthesis(Sinput,filters, options)

niters=getoptions(options,'niters',300);
periodinput=getoptions(options,'periodinput',0);
localized = getoptions(options, 'localized',0);
mgrid = getoptions(options,'multigrid',0);
lambda=getoptions(options,'lambda',1e-8);
debugging = getoptions(options,'debugging',false);
lr = getoptions(options,'lr',0.9);


[N1, N2]= size(filters{1}.psi{1}{1}{1});


%Initialization 
out=abs(randn(N1,N2));
%out=ifft2(fft2(out).*filters{1}.phi{1});
out=getoptions(options,'init',out);


%backtrack the second layer

for it=1:niters
    diffS = backtrack_scattering(out,Sinput,filters,options);

    x = x-lr*diffS;
    
end 

end 



function out=backtrack_scattering(out,Sinput,filters,options)

J = size(filters{1}.psi{1},2);
[S, U] = scattering2d(out, filters, options);

 diffS = scatdiff(S,Sinput);

    conv_filter=@(x,f)ifft2(applyfilter2d(fft2(x),f));

    indx2=1;
    sum_on_lambda1 = 0;
    for i=1:size(U{2},2) %for every lambda_1 in U{2}
        j1 = U{2}{i}.scale;
        l1 = U{2}{i}.orientation;

        sum_on_lambda2 = 0;
        while ((U{3}{indx2}.scale(1)==j1) && ... %second layer obtained from the same lambda1
               (U{3}{indx2}.orientation(1)==l1))

           j2 = U{3}{indx2}.scale(2);
           l2 = U{3}{indx2}.orientation(2);
           res= U{3}{indx2}.res-U{2}{i}.res;
           %compute (U{lambda1}{lamdba2}/abs(U{lambda1}{lamdba2})) conv conj(psi{lambda2})
           backtrack = real(upsample_and_apply_filter( U{3}{indx2}.signal./(eps+U{3}{indx2}.asignal),res, ...
                                         conj(filters{1}.psi{1}{j2}{l2}),J));

           %pointwise mult with ( v{lambda1}{lambda2} conv conj(phi) )                          
           lower_frq_component = upsample_and_apply_filter(diffS{3}{indx2}.l1,J, conj(filters{1}.phi{1}),J);

           sum_on_lambda2 = sum_on_lambda2+ real(backtrack.*lower_frq_component);
           indx2=indx2+1;
        end 

        % backtrack for lambda1: backtrack pointwise_mult U{lambda1}/abs(U{lambda1})
        deriv_sign = U{2}{i}.signal./(eps + U{2}{i}.asignal);

        accum_lambda2 = sum_on_lambda2.* deriv_sign;
        
        backtrackL1L2 = real(upsample_and_apply_filter(accum_lambda2,U{2}{i}.res, ...
                                       conj(filters{1}{j1}{l1})),J);

        backtrackL1 = upsample_and_apply_filter(diffS{2}{i}.l1,J-(U{2}{i}.res+1), conj(filters{1}.phi),J)  ; 
        
        backtrackL1 = real(upsample_and_apply_filter(backtrackL1.* deriv_sign,U{2}{i}.res,...
                                                     conj(filters{1}.psi{1}{j1}{l1})),J);

        sum_on_lambda1 = sum_on_lambda1+ backtrackL1L2+backtrackL1;
    end 

    out = conv_filter(diffS{1}{1}.l1, conj(filters{1}.phi)) + sum_on_lambda1;

end 

function ch = upsample_and_apply_filter(signal,res, filter, J)

%res = 1+res;
%ids_f = 2^(max(0,floor(J-res)));
ids_f =2^res;
%upsample
ch = zeros(size(signal,1)*ids_f,size(signal,2)*ids_f);
ch(1:ids_f:end,1:ids_f:end) = 2^(J-1)*signal;

%apply filter phi
ch = ifft2(applyfilter2d(fft2(ch),filter));

end
% function ch = upsample_and_apply_phi(signal,res, phi, J)
% 
% res = 1+res;
% ids_f = 2^(max(0,floor(J-res)));
% 
% %upsample
% ch = zeros(size(signal,1)*ids_f,size(signal,2)*ids_f);
% ch(1:ids_f:end,1:ids_f:end) = 2^(J-1)*signal;
% 
% %apply filter phi
% ch = ifft2(applyfilter2d(fft2(ch),phi));
% 
% end

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

end 




function out=blurbp(Stmp, Sref, scratch, filters, options,J)

local = getoptions(options, 'localized',1);
ninput = getoptions(options, 'ninput', numel(filters{1}.phi{end}));

res = 1+scratch.res;

ids_f = 2^(max(0,J-res-1));


    if local
        in = 2^(J-1)*(Stmp.l1 - Sref.l1);
        in = permute(in, [2 1]);
        ch = zeros(size(in,1)*ids_f,size(in,2));
        ch(1:ids_f:end,:) = in;
        ch = ifft(fft(ch).*repmat(conj(filters{1}.phi{res}(1,:)'), 1 , size(ch,2)));
        
        in = permute(ch, [2 1]);
        ch = zeros(size(in,1)*ids_f,size(in,2));
        ch(1:ids_f:end,:)=in;
        out = ifft(fft(ch).*repmat(conj(filters{1}.phi{res}(:,1)), 1 , size(ch,2)));
       
    else %delocalized
        out = ninput^2*ones(size(scratch.asignal))*(Stmp.l1 - Sref.l1)/numel(scratch.asignal) ;
       
    end
    

end


















function diffS = scatdiff(S1,S2)

diffS = S1;

for i=1:length(S1)
    for j=1:length(S1{i})
        diffS{i}{j}.l1 = (S1{i}{j}.l1-S2{i}{j}.l1);
    end 
end 

end 


function [dout,curdist, refe] = bkwscatt_mask(in, S0, filters, options)

osf = getoptions(options,'osf',0);
mgrid=getoptions(options,'multigrid',1);

%perform a forward transform first
[Stmp, scratch] = fwdscatt_dicts(in, filters, options);

[curdist,refe] = scdist(Stmp,S0,options);

[N1, N2]=size(in);

dout=zeros(N1,N2);

Q=filters{1}.Q;
J=size(filters{1}.psi{1},2);

mmax=size(Stmp,2);


for m=mmax:-1:1
    if mgrid
        [codesc_tmp, codeor_tmp] = getcode(Stmp{m});
        [codesc_ref, codeor_ref] = getcode(S0{m});
    end
    for r=1:size(Stmp{m},2)
        %init gradients at their corresponding resolution and backpropagate phi_J output
        if mgrid
            Ir = find( (codesc_tmp(r)== codesc_ref) & ( codeor_tmp(r) == codeor_ref));
            din{m}{r}=blurbp(Stmp{m}{r}, S0{m}{Ir(1)}, scratch.X{m}{r}, filters, options, J,osf);
        else
            din{m}{r}=blurbp(Stmp{m}{r}, S0{m}{r}, scratch.X{m}{r}, filters, options, J,osf);
        end
        
        %backpropagate upper layer
        if m < mmax
            aux = 0;
            res = scratch.X{m}{r}.res;
            for rr=scratch.X{m}{r}.childbot:scratch.X{m}{r}.childtop
                aux0 = 0*din{m}{r};
                
                if m ==2 
                    aux1 = din{m+1}{rr}.*(D{rr}* (scratch.X{m+1}{rr}.signal./(eps+scratch.X{m+1}{rr}.asignal)));
                else %m==1
                    aux1 = din{m+1}{rr}.*scratch.X{m+1}{rr}.signal./(eps+scratch.X{m+1}{rr}.asignal);
                end 
                
                ids = 2^(scratch.X{m+1}{rr}.res - scratch.X{m}{r}.res);
                if onedim
                    aux0(1:ids:end)=aux1;
                else
                    aux0(1:ids:end,1:ids:end)=aux1;
                end
                
                sc = scratch.X{m+1}{rr}.scale(end);
                or = scratch.X{m+1}{rr}.orientation(end);
                aux = aux + real(ifft2( fourier(aux0).*filters{m}.dpsi{1+res}{sc}{or}));
               % din{m}{r} = din{m}{r} + real(ifourier( fourier(aux0).*conj(filters{m}.psi{1+res}{sc}{or})));
                %din{m}{r} = din{m}{r} + real(ifft2( fft2(aux0).*(filters{m}.dpsi{1+res}{sc}{or})));
            end
            din{m}{r} = din{m}{r} + aux;
        end
    end
end

dout=din{1}{1};

end







