function Dl2=getDatabase_onY2(Y2,lambda2,J,Ni)
%This function saves from Y2 the values (||x conv Psi_j1| conv Psi_j2|), for all
%j1 (note the two modulus!)
%Since our scattering can be localized, we obtain the set of corresponding patches of size 2^(J-1) (repatch function).
%The size of the output depends on the number of lambda1 for the given lambda2 and
%Ni/2^(J-1)(number of patches). More specifically,
% --- Output: 
% Dl2: (num_lambda1 x size patch, number of patches)
% --- Input: 
% Y2: scratch.X
% lambda2=[j2 l2]
% J : global J 
% Ni: size of the image.

j2 = lambda2(1);
l2 = lambda2(2);

rast = 1;
%num_patches = Ni/(2^(J-1));

for i=1:length(Y2{3})
    if ((Y2{3}{i}.scale(2) == j2) && (Y2{3}{i}.orientation(2) == l2))
       % Dl2(rast,:,:) = repatch(Y2{3}{i}.asignal,num_patches);
        Dl2(rast,:) = Y2{3}{i}.asignal(:);
        rast = rast+1;
  %      [Y2{3}{i}.scale(2) Y2{3}{i}.orientation(2) size(Y2{3}{i}.signal)]
    end 
end
%Dl2 = reshape(Dl2,[size(Dl2,1)*size(Dl2,2),size(Dl2,3)]);

function v = repatch(signal,num_patches)

N = size(signal);
w = N/num_patches;%patch size

flat=@(x)x(:);
X = 1:N/num_patches:(N(1)-w+1);
Y = 1:N/num_patches:(N(2)-w+1);

for i=0:w-1
    for j=0:w-1
        v(i+1,j+1,:)=flat(signal(X(:)+i,Y(:)+j));
    end 
end
v = reshape(v,size(v,1)*size(v,2),size(v,3));
