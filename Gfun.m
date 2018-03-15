function y = Gfun(A,rs,N_rvs)
    G_xt(:,1) = A - rs;
    G_xt(:,2) = A + rs;  
    y=min(G_xt,[],2);
% y=A-rs';
end