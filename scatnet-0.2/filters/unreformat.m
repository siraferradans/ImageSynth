function S = unreformat(in, Sref, options)

rast=1;
S = Sref;
dims = size(Sref{2}{1}.l1);

for m=1:size(Sref,2)
    for l=1:size(Sref{m},2)
        S{m}{l}.l1 = reshape(in(rast,:),dims); rast=rast+1;
    end
end



