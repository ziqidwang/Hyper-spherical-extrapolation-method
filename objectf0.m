function [out] = objectf0(x)
global Pf r Dim wei

for j=1:numel(r)
b=x;
es(j)=0.5*sum(betainc(1-(b/r(j)).^2,(Dim-1)/2,0.5));
end
% out=sum(((log(es)-log(Pf))./log(Pf)).^2);
 out=sum(wei.*((log(es)-log(Pf))).^2);

% if nargout > 1 % gradient required
% for j=1:numel(r)
% grd(j,:)=(log(sum(betainc(1-(x/r(j)).^2,(Dim-1)/2,0.5)))/log(Pf(j))-1)/log(Pf(j))./sum(betainc(1-(x/r(j)).^2,(Dim-1)/2,0.5)).*(x/r(j)).^(-1).*(1-(x/r(j)).^2).^((Dim-1)/2-1).*(-2.*x/r(j));
% end  
% g=sum(grd);    
% end
end




