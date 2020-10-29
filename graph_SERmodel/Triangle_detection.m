function [t0,t00,t01,t02,t11,t12,t13] = findtriangles(cij)

R = length(cij);
t0 = cij*cij';
maxto = max(t0(~eye(R)));

t0 = t0.*cij;
t00 = zeros(R,R,2); t01 = zeros(R,R,maxto,2); t02 = zeros(R);
t11 = zeros(R,R,maxto,4); t12 = zeros(R,R,maxto*(maxto-1)/2,3); t13 = zeros(R);

for i=1:R
    for j=i+1:R
        ntmp = double(sum(cij([i j],:))==2); ntmp([i j]) = -1;
        idn = find(ntmp==1); idnn = find(ntmp==0); sidn = length(idn);

        t00(i,j,:) = diag(cij([i j],idnn)*cij(idnn,idnn)*cij(idnn,[i j]))/2;
        for k=1:sidn, t01(i,j,k,:) = diag(cij([i j],idn(k))*cij(idn(k),idnn)*cij(idnn,[i j])); end
        t02(i,j) = trace(cij([i j],idn)*cij(idn,idn)*cij(idn,[i j]))/2;

        if sidn>0
            s00 = find(ntmp==0&cij(i,:)==0&cij(j,:)==0);
            s10 = find(ntmp==0&cij(i,:)==1&cij(j,:)==0);
            s01 = find(ntmp==0&cij(i,:)==0&cij(j,:)==1);

            t11(i,j,1:sidn,1) = diag(cij(idn,s00)*cij(s00,s00)*cij(s00,idn))/2;
            t11(i,j,1:sidn,2) = diag(cij(idn,s10)*cij(s10,s00)*cij(s00,idn));
            t11(i,j,1:sidn,3) = diag(cij(idn,s01)*cij(s01,s00)*cij(s00,idn));
            t11(i,j,1:sidn,4) = diag(cij(idn,s10)*cij(s10,s01)*cij(s01,idn));

            if sidn>1                
                t12(i,j,1:sidn*(sidn-1)/2,1) = squareform((cij(idn,s00)*cij(s00,idn)).*cij(idn,idn));
                t12(i,j,1:sidn*(sidn-1)/2,2) = squareform((cij(idn,s10)*cij(s10,idn)).*cij(idn,idn));
                t12(i,j,1:sidn*(sidn-1)/2,3) = squareform((cij(idn,s01)*cij(s01,idn)).*cij(idn,idn));
            end
            t13(i,j) = trace(cij(idn,idn)^3)/6;
        end
        t00(j,i,:) = t00(i,j,:); t01(j,i,:,:) = t01(i,j,:,:);
        t11(j,i,:,:) = t11(i,j,:,:); t12(j,i,:,:) = t12(i,j,:,:);
    end
end
t02 = t02+t02';
t13 = t13+t13';
