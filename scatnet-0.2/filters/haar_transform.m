function H = haar_transform(N, options)
%we implement a 1D haar transform of N points

J=floor(log2(N));

rast=1;
for j=1:J
atom=zeros(1,2^j);
atom(1:2^(j-1))=1;
atom(2^(j-1)+1:2^(j))=-1;
for n=1:N-2^j+1
H(rast,n:n+2^j-1)=2^(-j)*atom;rast=rast+1;
end
end

atom=ones(1,2^J);
for n=1:N-2^J+1
H(rast,n:n+2^j-1)=2^(-J)*atom;rast=rast+1;
end



