function X = randUniformOnSphere(d,n)

%RANDUNIFORMONSPHERE  Draw n sample from the uniform distribution on the unit
%                     sphere S^{d-1} in R^{d}.
%
% The output X is a d-by-n matrix.

if nargin<2, n=1; end
X = randn(d,n);
X = X./repmat(sqrt(sum(X.^2,1)),[d, 1]);
