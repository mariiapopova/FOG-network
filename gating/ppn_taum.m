function tau=ppn_taum(V)
alpha=0.32.*(V+55)./(1-exp(-(V+55)./4));
beta=-0.28.*(V+28)./(1-exp((V+28)./5));
tau=1./(alpha+beta);
return