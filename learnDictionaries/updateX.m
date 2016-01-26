function X = updateX(Y,D,k,X,verbose)

%Update of the Coefficients X
select = @(A,k)repmat(A(k,:), [size(A,1) 1]);
ProjX = @(X,k)X .* (abs(X) >= select(sort(abs(X), 'descend'),k));

epsilon = 1e-3;
flat=@(x)x(:);
t = 1.5/(norm(flat(D*D')) + epsilon);

norm2=@(D)sqrt(sum((abs(D(:)).^2)));

it = 10000; 
for j=1:it
    X = ProjX(X-t*D'*(D*X-Y),k);%soft thresholding on the grad desc
    
    %for debugging
    Err(j) = norm2(Y-D*X);
   if verbose
    if (j>1) 
        if Err(j-1)-Err(j) < 1e-8
            subplot(1,3,2);plot(log10(Err),'-');drawnow;
            return;
        end 
    end
   end
  
end
if verbose
   subplot(1,3,2);plot(log10(Err),'*-');drawnow;
end   

end 