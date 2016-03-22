function h=visualizing_dict(d)

%show both, real and imaginary part
%sort per peak of maximum frequency 

%real part
[~, peak_lambda1s] = max(real(d).^2, [], 1);

[~, sorting_indices] = sort(peak_lambda1s);

h=figure; 
subplot(121);
imagesc(real(d(:,sorting_indices)));title(['Real dict '])

[~, peak_lambda1s] = max(imag(d).^2, [], 1);

[~, sorting_indices] = sort(peak_lambda1s);

subplot(122);
imagesc(imag(d(:,sorting_indices)));title(['Imag dict ']);