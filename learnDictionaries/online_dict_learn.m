function [D,X,err]=online_dict_learn(Y,k,n,p)
%
% energy: Y-D*X, with target sparsity k, X is of w length
% Y: (d,m) with m patches of dim. d
% k: maximum sparsity
% n: n dims of 
% p: number of vectors that we want to use for the training (p with maximum
%    norm)

verbose = 1;
N = size(Y,2);
Nminibatch = min(256,N);
[d,m]=size(Y); %num. of vectors used for the training

if nargin<4
    p=m;
end 

%Only keep the patches with largest energy.
%[~,I] = sort(sum(Y.^2), 'descend'); %todo:randomly select on each it.
% I = randperm(m);
% Y = Y(:,I(1:p));

eps=1e-5;
ProjC = @(D)D ./ repmat( sqrt(sum(D.^2)+eps), [d 1]);% [w, 1] );

%select some subset of the data to initialize the dictionary
sel = randperm(p); sel = sel(1:n);
D =ProjC(sparse(double(Y(:,sel))));

%init
X = rand(size(D,2),Nminibatch);

At = zeros(size(D,2),size(D,2));
Bt = zeros(size(Y,1),size(D,2));

flat=@(x)x(:);
norm2=@(D)sqrt(sum(flat(abs(D).^2)));

niter=100;
for it=1:niter
    progressbar(it,niter,20);
    
    %obtain randomly a set of elements from the database
    indx = randperm(N);
    Yt = Y(:,indx(1:Nminibatch));
   
    %get coeffs
    X = updateX(Yt,D,k,X,verbose); 
    
    % data to update Dict
    At = At+X*X';
    Bt = Bt+Yt*X';
    D = updateD(At,Bt,D,verbose);

    err(it) = norm2(Yt-D*X);
    if verbose
        subplot(1,3,3);plot(log10(err+eps));drawnow
    end 
    if err(it)<1e-4
        return;
    end 
   
end
if verbose
    subplot(1,3,3);plot(log10(err+eps));
end 
% hold on;

% [~,I] = sort(sum(X.^2,2), 'descend'); %todo:randomly select on each it.
% 
% X = X(I,:);
% D = D(:,I);

end     

function D = updateD(A,B,D,verbose)

[d,~]=size(D);
epsilon = 1e-3;
ProjC = @(D)D ./ repmat( max(sqrt(sum(D.^2)),1), [d, 1] );% d_j = d_j * (1/max(norm(d_j))_2,1)

norm2=@(D)sqrt(sum((abs(D(:)).^2)));

it = 100;
iX = diag(1./diag(A));
t = 1.8/(norm((A*iX)*(A*iX)')+epsilon);% + lambda*k+ epsilon);

for j=1:it
  %  D = ProjC( D-t*(D*X-Y)*X' );
    D = ProjC(D+t*(B-D*A)*iX);
    %for debugging
    Err(j) = mean( 0.5*trace(D'*D*A) - trace(D'*B) );
    if verbose
        if (j>1) 
            if Err(j-1)-Err(j) < 1e-8
                  subplot(1,3,1);plot(log10(Err),'-');drawnow;

                return;
            end 
        end
    end 
end
if verbose
    subplot(1,3,1);plot(log10(Err),'*-');drawnow;
end 
end 


