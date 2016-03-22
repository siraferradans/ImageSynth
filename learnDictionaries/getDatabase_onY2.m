function Dl2=getDatabase_onY2(Y2,lambda2)
%This function saves from Y2 the values (|x conv Psi_j1| conv Psi_j2), for all
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

flat=@(x)x(:);

for i=1:length(Y2{3})
    if ((Y2{3}{i}.scale(2) == j2) && (Y2{3}{i}.orientation(2) == l2))
       % Dl2(rast,:,:) = repatch(Y2{3}{i}.asignal,num_patches);
        Dl2(rast,:) = flat(Y2{3}{i}.signal(1:2:end,1:2:end)); 
        %we know there is overlapping between the windows, we get only 1 every 4.
        
        rast = rast+1;
  %      [Y2{3}{i}.scale(2) Y2{3}{i}.orientation(2) size(Y2{3}{i}.signal)]
    end 
end
