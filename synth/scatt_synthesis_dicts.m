function [out,total_cost]=scatt_synthesis_mask(S,Mask, filters, options, inmax)

niters=getoptions(options,'niters',300);
rho=getoptions(options,'rhotol',0.02);
periodinput=getoptions(options,'periodinput',0);
localized = getoptions(options, 'localized',0);
mgrid = getoptions(options,'multigrid',0);
lambda=getoptions(options,'lambda',1e-8);
lambdabis=getoptions(options,'lambdabis',1e-4);
onedim=getoptions(options,'onedim',0);
softthreshold=getoptions(options,'softthreshold',0);
optm = getoptions(options,'optmethod','gradient_desc');
debugging = getoptions(options,'debugging',false);
mom = getoptions(options,'momentum',0.2);
lr = getoptions(options,'lr',0.9);
do_TV = getoptions(options,'do_TV',false);
if ~do_TV 
    lambdabis = 0;
end 

flat = @(x)x(:);
%rng('shuffle');
rand('state',sum(100*clock));

[N1, N2]= size(filters{1}.psi{1}{1}{1});

%using mirror symmetry: cut by half these numbers
if options.mirror | periodinput
    N1=N1/2; N2=N2/2;
end

%Initialization 
out=abs(randn(N1,N2));
%out=ifft2(fft2(out).*filters{1}.phi{1});
out=getoptions(options,'init',out);

% Sbn =norm( reformat(S,options));
% out = out*Sbn / norm(out(:));

state = init_state(options, N1, N2);



Q=filters{1}.Q;
J=size(filters{1}.psi{1},2)/Q;
total_cost=[];
for j=mgrid*(J-1):-1:0
    cond = 1;
    options.minj=j;
    %nS = scnorm(S,options);
    it=1;
    options.momentum = mom;
    options.lr = lr;
    
    if debugging
        fprintf('doing scale %d \n', j)
    end 
    while (it<niters)&(cond==1)
        
       
        [dout,curdist,nS]= bkwscatt_mask(out, S,Mask,  filters, options);

        %complete computation of the gradient:
        %adding ridge regression or Tikhonov regularization
        dout = dout + lambda * out;

        xtmp = out - lr*dout;
        %[xtmp, state_tmp] = update_params(out, state, dout, options);
        
       % curdist = scdist(fwdscatt(xtmp, filters, options),S,options);
        
        if do_TV
            % Take the prox step.
            % NOTE: to deal with the non-periodicity of the TV implementation,
            % we extend xtmp periodically and then keep just a middle
            % portion. I currently take a large margin, but it can
            % probably made only a couple of pixels wide. In the end, this
            % step is so fast that it's inconsequential.
            
            xtmp = [xtmp(:, (N2/2+1):end) xtmp xtmp(:, 1:(N2/2))];
            xtmp = [xtmp((N1/2+1):end, :); xtmp; xtmp(1:(N1/2+1), :)];
         
            xtmp = TV(xtmp, lambdabis, 1, 1);
            xtmp = xtmp((N1/2+1):(N1/2+N1), (N2/2+1):(N2/2+N2));
        end
        
        
       % E = sum(curdist(:).^2) + lambdabis * sum((xtmp(:)<0).*(xtmp(:).^2));
         
         if lambdabis
            E = sum(curdist(:).^2) + lambda * sum((xtmp(:)<0).*(xtmp(:).^2)) + ...
                                   lambdabis * sum(abs(flat(diff(xtmp))));
         else 
             Stmp =fwdscatt(xtmp, filters, options) ;
             E=mean(sc_maskeddist(Stmp, S,Mask, options));

   %         E = sum(curdist(:).^2) + lambda * sum((xtmp(:)<0).*(xtmp(:).^2));
         end 
       
            out = xtmp;
         
%            state = state_tmp;
            total_cost(end+1) = E;

            %clamp to input min and max values (for images)
            out=max(-inmax,min(inmax,out));

%        end 
        
        if debugging
            if mod(it,5-1)==0
                figure(1);
                m = reformat(Mask,options);
                subplot(2,2,1);imagesc(real(out));colormap gray;drawnow;
                subplot(2,2,2);
                hold on;plot(log10(real(total_cost)),'g');
               % subplot(1,3,3);
                
                Svout=fwdscatt(out,filters,options);
                svout = reformat(Svout,options);
                U = size(Stmp{2}{1}.l1);
              %  svout = reshape(svout(2:J*Q+1,:),Q,J,U(1),U(2));
                
                svin = reformat(S,options);
              %  svin = reshape(svin(2:J*Q+1,:),Q,J,U(1),U(2));
               
                subplot(223);plot(svout(:));hold on;plot(svin(:),'r')
                
%                 
%                 auxss=reformat(S,options);
%                 if (size(m,2) ==1)
%                     svout=repmat(m,[1 size(svout,2)]).*svout;
%                     ss = repmat(m,[1 size(svout,2)]).*reformat(S,options);
%                 else
%                     svout = m.*svout;
%                     ss = m.*reformat(S,options);
%                 end 
%                 [~,indx] = max(mean(svout,1));
%                 %indx =1:size(svoutaux,2);
%                 % plot_scattering(cat(2,svout(:,indx),ss(:,indx)),J,options.L1);
%                 subplot(2,2,3);
%                 %plot(svout(2:J*options.L1+1,indx),'b');hold on; plot(ss(2:J*options.L1+1,indx),'r');
%                 plot(svoutaux(2:J*options.L1+1,indx),'b');hold on; plot(auxss(2:J*options.L1+1,indx),'r');
%                 hold off;
%                 drawnow;
%                  subplot(2,2,4);
%                  %plot(svout(J*options.L1+2:end,indx),'b');hold on; plot(ss(J*options.L1+2:end,indx),'r');
%                  plot(svoutaux(J*options.L1+2:end,indx),'b');hold on; plot(auxss(J*options.L1+2:end,indx),'r');
%                  hold off;
               drawnow;
            end
       end 
       cond = max( (curdist(:) > rho*nS(:)) )==1;
       
        it=it+1;
    end
end
if ~cond
    disp(['scatt_synth_mask: Getting out because cond=0 at it=' num2str(it)])
end 

end

function state = init_state(options, N1, N2)

method=getoptions(options,'optmethod','mom');

if isequal(method,'BFGS')
    state.B=eye(N1*N2);
    state.d=zeros(N1*N2,1);
    state.s=state.d;
    state.y=state.d;
    state.kick=0;
else
    state = zeros(N1,N2);
end

end


function [out, stateout] = update_params(in, statein, din, options)

method=getoptions(options,'optmethod','mom');
lr=getoptions(options,'lr',0.1);
rho=getoptions(options,'momentum',0.9);


if strcmp(method,'grad_desc')
  rho = 0;     
end

%momentum/nesterov/gradient descent
stateout = rho * statein - lr * din;
out = in + stateout;


end

function [dout,curdist, refe] = bkwscatt_mask(in, S0,Mask, filters, options)

onedim = getoptions(options,'onedim',0);
destroy = getoptions(options,'destroy',0);
osf = getoptions(options,'osf',0);
coef_stretching = getoptions(options,'coef_stretching',1);

map=@(scale)scale*coef_stretching;

%perform a forward transform first
[Stmp, scratch] = fwdscatt_dicts(in, filters, options);

[curdist,refe] = scdist(Stmp,S0,options);

D = filters.D;

[N1, N2]=size(in);
dout=zeros(N1,N2);

Q=filters{1}.Q;
J=map(size(filters{1}.psi{1},2)/Q);
mgrid=getoptions(options,'multigrid',1);
mmax=size(Stmp,2);

if onedim
    fourier=@fft;
    ifourier=@ifft;
else
    fourier=@fft2;
    ifourier=@ifft2;
end

preco = getoptions(options,'precond',ones(mmax,1));
preco = preco.^2;

for m=mmax:-1:1
    if mgrid
        [codesc_tmp, codeor_tmp] = getcode(Stmp{m});
        [codesc_ref, codeor_ref] = getcode(S0{m});
    end
    for r=1:size(Stmp{m},2)
        %init gradients at their corresponding resolution and backpropagate phi_J output
        if mgrid
            Ir = find( (codesc_tmp(r)== codesc_ref) & ( codeor_tmp(r) == codeor_ref));
            din{m}{r}=blurbp_mask(Stmp{m}{r}, S0{m}{Ir(1)},Mask{m}{Ir(1)}, scratch.X{m}{r}, filters, options, J,preco(m),osf);
        else
            din{m}{r}=blurbp_mask(Stmp{m}{r}, S0{m}{r},Mask{m}{r}, scratch.X{m}{r}, filters, options, J,preco(m),osf);
        end
        if m==mmax & destroy
            din{m}{r} = -2e-1 * din{m}{r};
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
                aux = aux + real(ifourier( fourier(aux0).*filters{m}.dpsi{1+res}{sc}{or}));
               % din{m}{r} = din{m}{r} + real(ifourier( fourier(aux0).*conj(filters{m}.psi{1+res}{sc}{or})));
                %din{m}{r} = din{m}{r} + real(ifft2( fft2(aux0).*(filters{m}.dpsi{1+res}{sc}{or})));
            end
            din{m}{r} = din{m}{r} + aux;
        end
    end
end

dout=din{1}{1};

end


function out=blurbp_mask(Stmp, Sref,Mask, scratch, filters, options,J, prec,osf)

onedim = getoptions(options,'onedim',0);
local = getoptions(options, 'localized',1);
%ninput = getoptions(options, 'ninput', numel(filters{1}.phi{3}));
ninput = getoptions(options, 'ninput', numel(filters{1}.phi{end}));

res = 1+scratch.res;
%ids_f = 2^(max(0,J-res-osf));
ids_f = 2^(max(0,floor(J-res-osf)));

if onedim
    if local
        in = prec*2^((J-1)/2)*(Mask.l1).*(Stmp.l1 - Sref.l1);
        ch = zeros(size(in,1)*ids_f,1);
        ch(1:ids_f:end) = in;
        out = ifft(fft(ch).*conj(filters{1}.phi{res}));
       
    else %delocalized
        out = ninput^2*ones(size(scratch.asignal))*(Mask.l1).*(Stmp.l1 - Sref.l1)/numel(scratch.asignal) ;
    end
else
    if local
        in = prec*2^(J-1)*(Mask.l1).*(Stmp.l1 - Sref.l1);
        in = permute(in, [2 1]);
        ch = zeros(size(in,1)*ids_f,size(in,2));
        ch(1:ids_f:end,:) = in;
       % ch = ifft(fft(ch).*repmat(conj(filters{1}.phi{res}(:,1)), 1 , size(ch,2)));
         ch = ifft(fft(ch).*repmat(conj(filters{1}.phi{res}(1,:)'), 1 , size(ch,2)));
        
        in = permute(ch, [2 1]);
        ch = zeros(size(in,1)*ids_f,size(in,2));
        ch(1:ids_f:end,:)=in;
        out = ifft(fft(ch).*repmat(conj(filters{1}.phi{res}(:,1)), 1 , size(ch,2)));
       
    else %delocalized
        out = ninput^2*ones(size(scratch.asignal))*(Mask.l1).*(Stmp.l1 - Sref.l1)/numel(scratch.asignal) ;
       
    end
    
    
end

end

function [out,rout]=sc_maskeddist(Stmp, Sref,Mask, options)

mmax=size(Stmp,2);
usel2 = getoptions(options,'l2scatt',1);
mgrid=getoptions(options,'multigrid',1);
preco = getoptions(options,'precond',ones(mmax,1));

if usel2
out=zeros(mmax,2);
rout=zeros(mmax,2);
else
out=zeros(mmax,1);
rout=zeros(mmax,1);
end

for m=1:mmax
    if mgrid
        [codesc_tmp, codeor_tmp] = getcode(Stmp{m});
        [codesc_ref, codeor_ref] = getcode(Sref{m});
    end
    for r=1:size(Stmp{m},2)
        if mgrid
            Ir = find( (codesc_tmp(r)== codesc_ref) & ( codeor_tmp(r) == codeor_ref));
            out(m,1)=out(m,1)+preco(m)^2*sum((Mask{m}{r}.l1(:)).*(Stmp{m}{r}.l1(:)-Sref{m}{Ir(1)}.l1(:)).^2);
            rout(m,1)=rout(m,1)+preco(m)^2*sum((Mask{m}{Ir(1)}.l1(:)).*(Sref{m}{Ir(1)}.l1(:)).^2);
            
        else
            out(m,1)=out(m,1)+preco(m)^2*sum((Mask{m}{r}.l1(:)).*(Stmp{m}{r}.l1(:)-Sref{m}{r}.l1(:)).^2);
            rout(m,1)=rout(m,1)+preco(m)^2*sum((Mask{m}{r}.l1(:)).*(Sref{m}{r}.l1(:)).^2);
            
        end
    end
end
out=sqrt(out);
rout=sqrt(rout);

end




