function tau=ppn_taumk(V)
alpha=0.032*(V+63.8)./(1-exp(-(V+63.8)./5));
beta=0.5*(exp(-(V+68.8)./40));
tau=1./(alpha+beta);
return