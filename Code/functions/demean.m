function x=demean(X);

[T,N]=size(X);
m=mean(X);
x=X-repmat(m,T,1);
