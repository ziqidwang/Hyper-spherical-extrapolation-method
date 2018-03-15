function [out]=extrap(Main,j)
N_rvs = Main.numb;%% number of random variables
N=100;
Ns=200;
Main.R=sqrt(N_rvs); %% sampling radius
count=0; %% record the number of G function evaluations
nt=numel(Main.Amplitude); %% number of thresholds
nr=6;
ratl=0.03;
rath=0.06;
convtol=0.12;
global K Pf r Dim wei
Dim=N_rvs;
%% generate an initial  N samples 
Main.A=Main.Amplitude(end);
if j==1
Main.R=Main.R*1.1;
while 1==1
u=Main.R*randUniformOnSphere(N_rvs,Ns);
rs=Responsef(u,Main);
count=count+Ns;
H_x=Gfun(Main.A,rs,N_rvs);
   H_x=(H_x<=0);
   if sum(H_x)>=Ns*ratl-1
       break
   else
       Main.R=Main.R*1.1;
   end
end
r(1)=Main.R;
Main.R=r(1)*1.1;
while 1==1
u=Main.R*randUniformOnSphere(N_rvs,Ns);
rs=Responsef(u,Main);
count=count+Ns;
H_x=Gfun(Main.A,rs,N_rvs);
   H_x=(H_x<=0);
   if sum(H_x)>=Ns*rath-1
       break
   else
       Main.R=Main.R*1.1;
   end
end
r(nr)=Main.R;
r=r(1):(r(end)-r(1))/(nr-1):r(end);
out.r=r;
else
r=Main.Rs;
end
l=0;
while 1
    l=l+1;
    u0=randUniformOnSphere(N_rvs,N);
    id=1:1:N;
    us=u0;
for j=1:nr   
    u=r(j)*us;
    us=u(:,id)/r(j);
    rs=Responsef(us*r(j),Main);
    count=count+size(us,2);
    H_x=Gfun(Main.A,rs,N_rvs);
    H_x=(H_x<=0);
    id=~H_x;
    pf(l,j)=(sum(H_x)+N-size(us,2))/N;
    Pf(j)=mean(pf(:,j));
end
    if l>=(1-ratl)/ratl/convtol^2/N
      conv=std(pf(:,1))/Pf(1)/sqrt(l);
        if conv<=convtol
            break
        end  
    end
end 

Cov=sqrt((1-Pf)./Pf/l/N);
wei=(log((Pf+1.96*Pf.*Cov)./(Pf-1.96*Pf.*Cov))).^(-2);

for K=1:100
lb=[zeros(1,K)];
ub=[(sqrt(N_rvs)-3.5)*ones(1,K)];
x0=[3.0*ones(1,K)]; 
y1{K}=fmincon(@objectf0,x0,[],[],[],[],lb,ub);
conv(K)=objectf0(y1{K});
if K>2
    if conv(K-2)>conv(K-1)&&conv(K)>conv(K-1)
        id1=K-1;
        out.convid1=conv(id1);
        break
       end 
end
if K==3
    if conv(K-1)>conv(K-2)&&conv(K)>conv(K-2)
    id1=K-2;
    out.convid1=conv(id1);
    break
    end
end
end
id2=id1;
for K=1:id1
lbt=[zeros(1,K)];
ubt=[(sqrt(N_rvs)-3.5)*ones(1,K)];
x0t=[3.0*ones(1,K)]; 
yt{K}=fmincon(@objectf0,x0t,[],[],[],[],lbt,ubt);
lb=[-yt{K}/sqrt(N_rvs)/5,zeros(1,K)];
ub=[yt{K}/sqrt(N_rvs)/5,(sqrt(N_rvs)-3.5)*ones(1,K)];
% lb=[-1-min(yt{K})*ones(1,K)/r(end),zeros(1,K)];
% ub=[1-max(yt{K})*1.0*ones(1,K)/(sqrt(N_rvs)-3.5),(sqrt(N_rvs)-3.5)*ones(1,K)];
x0=[zeros(1,K),yt{K}];
options=optimoptions('fmincon','MaxFunEvals',5000);
y2{K}=fmincon(@objectf1,x0,[],[],[],[],lb,ub,[],options);
conv(K)=objectf1(y2{K});
if K>2
    if conv(K-2)>conv(K-1)&&conv(K)>conv(K-1)
        id2=K-1;
        out.convid2=conv(id2);
        break
       end 
end
if K==3
    if conv(K-1)>conv(K-2)&&conv(K)>conv(K-2)
    id2=K-2;
    out.convid2=conv(id2);
    break
    end
end
if K==id1
    out.convid2=conv(K);
end
end

x1=y1{id1};
x2=y2{id2};
rt=sqrt(N_rvs)-3.5:(r(1)-sqrt(N_rvs)+3.5)/100:r(1);
rd=[rt(1:end-1),r];
Pfd=[zeros(1,100),Pf];
for j=1:numel(rd)
b1=x1;
b2=x2(1:id2)*rd(j)+x2(id2+1:2*id2);
% b=(-x(1:id)+sign(x(1:id)).*sqrt(x(1:id).^2+4*x(1:id).*x(id+1:2*id)+4*r(j)^2))/2;
es1(j)=0.5*sum(betainc(1-(b1/rd(j)).^2,(Dim-1)/2,0.5));
es2(j)=0.5*sum(betainc(1-(b2/rd(j)).^2,(Dim-1)/2,0.5));
end

Pff1=sum(normcdf(-x1));
beta1=-norminv(Pff1);
pff=0;
upr=sqrt(N_rvs)+3.5;
lowr=sqrt(N_rvs)-3.5;
dr=(upr-lowr)/10000;
rs=lowr:dr:upr;
for j=1:numel(rs)
    b=x2(1:id2)*rs(j)+x2(id2+1:2*id2);
    if max(abs(b/rs(j)))>1
    lpdff=-Inf;    
    else
    lpdff=log(0.5)+log(sum(betainc(1-(b/rs(j)).^2,(N_rvs-1)/2,0.5)))+log(dr)+chipdf(rs(j),N_rvs);
    end
    pff=exp(lpdff)+pff;
end
Pff2=mean(pff);
beta2=-norminv(Pff2);
out.count=count;
out.beta1=beta1;
out.id1=id1;
out.Pf1=Pff1;
out.beta2=beta2;
out.id2=id2;
out.Pf2=Pff2;
end
