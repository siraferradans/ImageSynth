

addpath(genpath('/users/data/ferradan/ImageSynth/code/'));

%Obtain statistics on the videos 
L=8;

%get Image info
disp(['Get filters:'])
J = 3;
Ni= 128;%2^J;%sqrt(Na*2^(J-1));
[filters_image, options_image,lpal] = generate_translate_wavelets([Ni Ni], J, L, 'image');
%for reconstruction 

Totalnum_im = 1;
path='../../data/' ;

X=getDabase_forDicts(path,Totalnum_im,J,L,Ni,filters_image, options_image);

%k = J*L/2; %sparsity value
save('./Database.mat','X','Totalnum_im','J','Ni','L')

disp(['start computing dictionary'])
%X=D,alpha
 for j2=2:J
    for l2=1:L
        [auxD,auxalpha,err]= learn_dict(double(X{linearindx(lambda2)}),k);
        D{linearindx(lambda2)} = auxD;
        alpha{linearindx(lambda2)}=auxalpha;
        ['For lambda2=[' num2str(j2) ',' num2str(l2) '], err=' num2str(err(end))]
    end 
 end 
       
