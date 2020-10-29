function ck = findck(cij)

R = length(cij);
t0 = (cij*cij').*cij;
maxto = max(t0(~eye(R)));
C = clustering_coef_bu(cij)';
d = degrees_und(cij);

ck = zeros(R,R,maxto);
for i=1:R
    for j=1:R
        ntmp = double(sum(cij([i j],:))==2); ntmp([i j]) = -1; id = find(ntmp==1);
        ck(i,j,1:length(id)) = C(id).*d(id).*(d(id)-1)/2 - cij(i,j);
    end
end
