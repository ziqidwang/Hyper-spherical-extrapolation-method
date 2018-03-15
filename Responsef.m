function rs=Responsef(u,Main)
kap=-0.0;
rs=sum(u,1)/sqrt(Main.numb)+kap/4*(u(1,:)-u(2,:)).^2;
end