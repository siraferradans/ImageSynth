

%  addpath(genpath('/users/data/ferradan/ImageSynth/code/'));
% % 
% % %Obtain statistics on the videos 
%  L=8;
% % 
% % %get Image info
% disp(['Get filters:'])
% J = 4;
% Ni= 128;%2^J;%sqrt(Na*2^(J-1));
% [filters_image, options_image,lpal] = generate_translate_wavelets([Ni Ni], J, L, 'image');

% %Get representation of some image
% 
impath='../../data/ILSVRC2012_val_00035280.JPEG' ; %an image in the training database
in=double(imread(impath));
in = imresize(in(:,:,1),[Ni,Ni]);
[Sinput, U] = scattering2d(in, filters_image, options_image);

[out,total_cost]=scattering2d_synthesis(Sinput,filters_image, options_image);

%%

% Totalnum_im = 1000;
% path='../../data/' ;
% 
% X=getDabase_forDicts(path,Totalnum_im,J,L,Ni,filters_image, options_image);
% 
% save('./Database1000_J4.mat','X','Totalnum_im','J','Ni','L')
% 
% return;

%%
load('./Database1000_J4.mat');

linearindx = @(lambda)(lambda(1)-1).*L+lambda(2);
disp(['start computing dictionary'])

params.modeD=0; %atoms with norm 1
params.mode=1; %l1 norm on the coefs alpha
params.lambda=0.01;
params.iter = 1000;
params.batchsize=512;
params.numThreads=-1; % number of threads
params.verbose=false;

%X=D*alpha, we want to get the best dictionary and the alphas for our DB
 for j2=2:J
    for l2=1:L
        lambda2=[j2 l2];
        
        XX = X{linearindx(lambda2)};
        IndexP = randperm(size(XX,2));
        %contiguous patches are very similar, thus we need to randomize their
        %position
        params.K=size(XX,1); %squared dictionary

        data = double(XX(:,IndexP));
     
        D_real=mexTrainDL_Memory(real(data),params);
        D_imag=mexTrainDL_Memory(imag(data),params);
   
%         Lambda1 = size(X{linearindx(lambda2)},1);
%         n = ceil(Lambda1*2/3); %dim of the atoms value 2/3 of whole dimension
%         k = ceil(Lambda1*1/2); %sparsity param
%         disp(['k-sparsity parameter: ' num2str(k) ' on ' num2str(n) ' dims' ])
%         [auxD,auxalpha,err]= learn_dict(double(X{linearindx(lambda2)}),k,n);
%         
        h=visualizing_dict(D_real+1i*D_imag);
        figure(h);title(['dicts for j2=' num2str(j2) ' l2=' num2str(l2)]);
        
        D{linearindx(lambda2)} = D_real + 1i*D_imag;
    %    alpha{linearindx(lambda2)}=auxalpha;
    %    relerr= norm(X{linearindx(lambda2)}-auxD*auxalpha)/norm(X{linearindx(lambda2)});
    %    ['For lambda2=[' num2str(j2) ',' num2str(l2) '], err=' num2str(err(end)) 'err rel=' num2str(relerr) ]
    end 
 end 
 save('./Dictionaries.mat','D')
      
%Get representation of some image
filters_image{3} = D;

impath='../../data/ILSVRC2012_val_00035280.JPEG' ; %an image in the training database
in=double(imread(impath));
in = imresize(in(:,:,1),[Ni,Ni]);
[S, U] = scattering2d_dicts(in, filters_image, options_image);

%visualize layers
flat=@(x)x(:);
for i_n = 1:size(S{3}{end}.l1,1)
    visualization_S3(i_n,:) = flat(S{3}{end}.l1(i_n,:,:));
end 
imagesc(visualization_S3)
