function [sim] = HOP(A,L,order,weights)
[n,m]=size(A);
F=zeros(n,m);
for i=1:order
    F=F+A^i.*(weights^(i-1));
end
% sim=F;
[U,S,V]=svd(F);
lambda = diag(S);
lambda(L+1:length(lambda)) = 0;
S=diag(lambda);
sim=U*S*V';

end