function out=reformat(S, options)

gpu = getoptions(options, 'gpu', 0);
l2scatt =getoptions(options,'l2scatt',0);
localized = getoptions(options,'localized',0);
options.nscattcoeffs = nsccoeffs(S);

rast=1;
if localized
   d = length(S{1}{1}.l1(:)); 
else 
   d = 1;
end
if gpu
    out=zeros(options.nscattcoeffs,d,'gpuArray');
else
    out=zeros(options.nscattcoeffs,d);
end
 

for m=1:size(S,2)
    for l=1:size(S{m},2)
        if localized
            out(rast,:)=S{m}{l}.l1(:);rast=rast+1;
        else 
            out(rast)=mean(S{m}{l}.l1(:));rast=rast+1;
            if l2scatt
                out(rast)=mean(S{m}{l}.l2(:));rast=rast+1;
            end
        end 
    end
end
if gpu
    out = gather(out);
end
end



