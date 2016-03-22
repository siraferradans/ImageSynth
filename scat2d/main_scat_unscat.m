addpath(genpath('../scatnet-0.2/'));

L=8;
% 
% %get Image info
disp(['Get filters:'])
J = 4;
Ni=64;

%[filters_image, options_image,lpal] = generate_translate_wavelets([Ni Ni], J, L, 'image');
impath='../../data/ILSVRC2012_val_00035280.JPEG' ; %an image in the training database
in=double(imread(impath));
in = imresize(in(:,:,1),[Ni,Ni]);
in = in / norm(in(:));
[Sinput, U] = scattering2d(in, filters_image, options_image);

out=scattering2d_synthesis(Sinput,filters_image, options_image);

out = out-min(out(:))+min(in(:));
subplot(121);imagesc(in);
subplot(122);imagesc(out);