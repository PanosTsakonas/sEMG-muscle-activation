function dadt=activation(t,a,u)

dadt=0;
Tact=15*10^-3;
Tdeact=50*10^-3;
dadt=(ppval(u,t)./Tact+(1-ppval(u,t))./Tdeact).*(ppval(u,t)-a);

end