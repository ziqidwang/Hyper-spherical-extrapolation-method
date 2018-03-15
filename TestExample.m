clear all
clc
tic 
Main.numb=100;%% number of random variables
Main.Amplitude=3;%% Threshold value for the G function 
Main.Run=1; %% number of runs

for j=1:Main.Run
out=extrap(Main,j);
beta1(j)=out.beta1
Pf1(j)=out.Pf1;
id1(j)=out.id1;
beta2(j)=out.beta2
Pf2(j)=out.Pf2;
id2(j)=out.id2;
count(j)=out.count;
conv1(j)=out.convid1;
conv2(j)=out.convid2;
if j==1
Main.Rs=out.r;
end
end

toc
