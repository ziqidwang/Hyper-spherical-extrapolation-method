function [out] = objectf1(x)
global Pf r Dim K wei

for j=1:numel(r)
b=x(1:K)*r(j)+x(K+1:2*K);
% b=(-x(1:K)+sign(x(1:K)).*sqrt(x(1:K).^2+4*x(1:K).*x(K+1:2*K)+4*r(j)^2))/2;
es(j)=0.5*sum(betainc(1-(b/r(j)).^2,(Dim-1)/2,0.5));
end

% out=sum(((log(es)-log(Pf))./log(Pf)).^2);
out=sum(wei.*((log(es)-log(Pf))).^2);
end




