function X=getDabase_forDicts(path,Totalnum_im,J,L,Ni,filters_image, options_image)

linearindx = @(v)(v(1)-1).*L+v(2);

dirimages = dir([path '*.JPEG']);

disp('Initialize');
for j2=2:J
    for l2=1:L
        lambda2=[j2 l2];
         X{linearindx(lambda2)}= [];
    end 
end

disp('and compute DB');
for i_im = 1:min(size(dirimages,1),Totalnum_im)
    image = double(imread([path dirimages(i_im).name]));
    image = imresize(image(:,:,1),[Ni,Ni]);

    [~,Y2] = scattering2d(image,filters_image, options_image);

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