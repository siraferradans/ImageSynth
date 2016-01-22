

addpath(genpath('/users/data/ferradan/ImageSynth/code/'));

%Obtain statistics on the videos 
L=8;

%get Image info
disp(['Get filters:'])
J = 3;
Ni= 128;%2^J;%sqrt(Na*2^(J-1));
[filters_image, options_image,lpal] = generate_translate_wavelets([Ni Ni], J, L, 'image');
%for reconstruction 
path='../../data/' ;

Totalnum_im = 5;
dirimages = dir([path '*.JPEG']);

linearindx = @(v)(v(1)-1).*L+v(2);
k = J*L/2; %sparsity value

disp('Initialize']);
for j2=2:J
    for l2=1:L
        lambda2=[j2 l2];
         X{linearindx(lambda2)}= [];
    end 
end

disp('and compute DB']);
for i_im = 1:min(size(dirimages,1),Totalnum_im)
    image = double(imread([path dirimages(i_im).name]));
    image = imresize(image(:,:,1),[Ni,Ni]);

    [Si,Y2] = scattering2d(image,filters_image, options_image);

    for j2=2:J
        for l2=1:L
            lambda2=[j2 l2];
            aux=getDatabase_onY2(Y2,lambda2,J,Ni);
            X{linearindx(lambda2)}=cat(2,X{linearindx(lambda2)},aux);
           
         %   filters_image{3}{linearindx(lambda2)} = learn_dict(double(D),k);
        end 
    end 
    
    disp(['[image ' num2str(i_im) ']Size X:'])
    X

end 
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
       
