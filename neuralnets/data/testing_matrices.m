% testing the matrices

W = permute(x,[5 4 3 2 1]);
Wr = W(:,:,:,:,1);
Wi = W(:,:,:,:,2);

[Kp,K,w1,w2] = size(WWi);
WWi = reshape(Wi,Kp,K*w1*w2);
WWr = reshape(Wr,Kp,K*w1*w2);

pinvWWr = pinv(WWr);
H = pinvWWr*WWi;
[V,D]=eig(H);
d = diag(D);
d=sort(d,'descend');