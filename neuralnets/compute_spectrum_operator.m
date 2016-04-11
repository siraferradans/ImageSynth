function spectrum=compute_spectrum_operator(k,N)
%assert(numel(size(k))==3);
% number of element should be odd

n_zero=floor((N-size(k,1))/2);
k=padarray(k,[n_zero,n_zero]);

k=fft2(k);



for i=1:size(k,1)
    for j=1:size(k,2)
        K=squeeze(k(i,j,:,:));
        s(i,j,:)=sqrt(svd(K*K'));
    end
end

s=s(:);

spectrum=sort(s,'descend');
end