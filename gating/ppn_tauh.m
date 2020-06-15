function tau=ppn_tauh(V)
alpha=0.12.*exp(-(V+51)./18);
beta=4./(1+exp(-(V+28)./5));
tau=1./(alpha+beta);
return