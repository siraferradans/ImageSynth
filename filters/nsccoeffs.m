function rast=nsccoeffs(S)

rast=0;
for m=1:size(S,2)
for l=1:size(S{m},2)
rast=rast+1; 
end
end

end


