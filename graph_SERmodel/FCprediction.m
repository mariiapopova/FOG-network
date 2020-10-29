function [fc2, fc1] = analyticalSERd(y,cij,t0,t00,t01,t02,t11,t12)
if nargin<3, [t0,t00,t01,t02,t11,t12] = findtriangles(cij); end

R = length(cij);
rto = cij*cij';
ydiff = zeros(1,3); zdiff = zeros(3,2);
for z=1:3, ydiff(z) = prod(y(setdiff(1:3,z))); zdiff(z,:) = setdiff(1:3,z); end
zid = [2 3 1];

fc2 = zeros(R);
if nargout>1
    fc1 = fc2;
    ck = findck(cij);
end

for z=1:3
%% level 0            
    npt0 = (1-y(z)).^t0;
    npt00 = ones(R,R,2); npt00b = ones(R,R,2,2); npt00c = npt00b;
    npt01d0 = ones(R,R,2); npt01d1 = npt01d0; npt01d2 = ones(R,R,2,2);
    npt01d0b = npt01d2; npt01d0c = npt01d2; npt01d1b = npt01d2; npt01d2b = npt01d2;
    npt02d2 = (1-2*ydiff(z)).^(t02/2); npt02d2b = ones(R,R,2); npt02d2c = ones(R,R,2); npt02d2d = npt02d2c;

    for zz=1:2
        czz = setdiff(1:2,zz);
        npt00(:,:,zz) = (1-2*ydiff(z)).^t00(:,:,zz);
        npt01d0(:,:,zz) = prod(1-y(setdiff(1:3,[z zid(z)])).*(1-(1-y(zid(z))).^t01(:,:,:,zz)),3);
        npt01d1(:,:,zz) = prod(1-y(zid(z)).*(1-(1-y(setdiff(1:3,[z zid(z)]))).^t01(:,:,:,zz)),3);
        if isempty(setdiff(zid(zdiff(z,czz)),zdiff(zdiff(z,zz),:))), npt02d2b(:,:,zz) = (1-2*(1-cij).*ydiff(zdiff(z,zz))).^(t02/2);
        else                                                         npt02d2c(:,:,zz) = (1-2*(1-cij).*ydiff(zdiff(z,zz))).^(t02/2);
        end
        if isempty(find(zdiff(z,:)==zid(zdiff(z,zz)),1)), npt02d2d(:,:,zz) = npt02d2c(:,:,zz); end
        
        cijtmp = repmat(cij,[1 1 size(t01,3)]);
        cijmask0 = 1-double(length(unique([zdiff(z,:) setdiff(zdiff(zdiff(z,zz),:),zid(zdiff(z,zz)))]))==3).*cijtmp;
        cijmask1 = 1-double(length(unique([zdiff(z,:) zid(zdiff(z,zz))]))==3).*cijtmp;
        cijmask2 = 1-double(length(unique([zdiff(z,:) zdiff(zdiff(z,zz),1)]))==3).*cijtmp;
        cijmask3 = 1-double(length(unique([zdiff(z,:) zdiff(zdiff(z,zz),2)]))==3).*cijtmp;
        cijmask2b = zid(zdiff(z,czz))==zdiff(zdiff(z,zz),1);
        cijmask3b = zid(zdiff(z,czz))==zdiff(zdiff(z,zz),2);
        for zzz=1:2
            npt00b(:,:,zz,zzz) = (1-2*ydiff(zdiff(z,zz))).^t00(:,:,zzz);
            
            npt01d0b(:,:,zz,zzz) = prod(1-cijmask0.*y(setdiff(zdiff(zdiff(z,zz),:),zid(zdiff(z,zz)))).*(1-(1-y(zid(zdiff(z,zz)))).^t01(:,:,:,zzz)),3);
            if isempty(find(zdiff(z,:)==zid(zdiff(z,zz)),1))
                npt00c(:,:,zz,zzz) = npt00b(:,:,zz,zzz);
                npt01d0c(:,:,zz,zzz) = npt01d0b(:,:,zz,zzz);
            end
            
            npt01d1b(:,:,zz,zzz) = prod(1-cijmask1.*y(zid(zdiff(z,zz))).*(1-(1-y(setdiff(zdiff(zdiff(z,zz),:),zid(zdiff(z,zz))))).^t01(:,:,:,zzz)),3);
            
            npt01d2(:,:,zz,zzz) = prod(1-cijmask2*(1-cijmask2b).*y(zdiff(zdiff(z,zz),1)).*(1-(1-y(zdiff(zdiff(z,zz),2))).^t01(:,:,:,zzz))...
                                        -cijmask3*(1-cijmask3b).*y(zdiff(zdiff(z,zz),2)).*(1-(1-y(zdiff(zdiff(z,zz),1))).^t01(:,:,:,zzz)),3);
            npt01d2b(:,:,zz,zzz) = prod(1-cijmask2.*cijmask2b.*y(zdiff(zdiff(z,zz),1)).*(1-(1-y(zdiff(zdiff(z,zz),2))).^t01(:,:,:,zzz))...
                                         -cijmask3.*cijmask3b.*y(zdiff(zdiff(z,zz),2)).*(1-(1-y(zdiff(zdiff(z,zz),1))).^t01(:,:,:,zzz)),3);
        end
    end

%% level 1
    npt11 = prod(1-y(zid(z)).*(1-(1-2*ydiff(zid(z))).^t11(:,:,:,1)),3);
    npt11b = prod(1-y(z).*(1-(1-2*ydiff(z)).^t11(:,:,:,1))-y(setdiff(1:3,[z zid(z)])).*(1-(1-2*ydiff(setdiff(1:3,[z zid(z)]))).^t11(:,:,:,1)),3);
    npt11c = ones(R,R,2); npt11e = npt11c;
    for zz=1:2
        cijtmp = repmat(cij,[1 1 size(t11,3)]);
        cijmask = 1-double(length(unique([zdiff(z,:) zid(zdiff(z,zz))]))==3).*cijtmp;
        npt11c(:,:,zz) = prod(1-cijmask.*y(zid(zdiff(z,zz))).*(1-(1-2*ydiff(zid(zdiff(z,zz)))).^t11(:,:,:,1)),3);

        a = hist([zdiff(z,:) zid(zdiff(z,zz))],1:3);
        if ~isempty(find(a==0,1)), npt11e(:,:,zz) = (1-(1-cij).*y(a==0)).^rto;
        else                       npt11e(:,:,zz) = (1-(1-cij).*sum(y(zdiff(z,:)))).^rto;
        end
    end
    npt11d = prod(1-y(setdiff(1:3,zid(zdiff(z,:)))).*(1-(1-2*ydiff(setdiff(1:3,zid(zdiff(z,:))))).^t11(:,:,:,1)),3);
    npt11f = (1-(1-cij).*y(z)).^rto;

    npt12 = prod(1-2*y(zid(z)).*y(z).*(1-(1-y(setdiff(1:3,[z zid(z)]))).^t12(:,:,:,1)),3);
    npt12b = prod(1-2*y(setdiff(1:3,[z zid(z)])).*y(z).*(1-(1-y(zid(z))).^t12(:,:,:,1)),3);
    npt12c = prod(1-2*ydiff(z).*(1-(1-y(z)).^t12(:,:,:,1)),3);

%% FC2
    fc2 = fc2 + y(z)^2.*...
            (1-prod(npt01d1,3).*npt02d2.*npt11.*npt12.*...
               (1-(1-npt00(:,:,1).*npt01d0(:,:,1)).*(1-npt00(:,:,2).*npt01d0(:,:,2))).*...
               (1-prod(npt00,3).*prod(npt01d0,3).*(1-npt11b.*npt12b)));

    fc2 = fc2 + ydiff(z).*...
            npt0.*npt01d2b(:,:,1,1).*npt01d2b(:,:,2,2).*prod(npt02d2b,3).*...
            (1-(1-npt00b(:,:,1,1).*npt01d2(:,:,1,1).*npt02d2c(:,:,1)).*...
               (1-npt00b(:,:,2,2).*npt01d2(:,:,2,2).*npt02d2c(:,:,2))).*...
            (1-(1-npt11e(:,:,1).*(1-npt01d1b(:,:,1,1).*npt11c(:,:,1))).*...
               (1-npt11e(:,:,2).*(1-npt01d1b(:,:,2,2).*npt11c(:,:,2))).*...
               (1-npt00b(:,:,1,1).*npt01d0b(:,:,1,1).*npt00b(:,:,2,2).*npt01d0b(:,:,2,2).*prod(npt02d2c,3).*npt11f.*(1-npt11d)).*...
               (1-npt00c(:,:,1,1).*npt01d0c(:,:,1,1).*npt00c(:,:,2,2).*npt01d0c(:,:,2,2).*prod(npt02d2d,3).*npt11f.*(1-npt12c)));

    fc2 = fc2 + ydiff(z).*...
            npt0.*npt01d2b(:,:,2,1).*npt01d2b(:,:,1,2).*prod(npt02d2b,3).*...
            (1-(1-npt00b(:,:,2,1).*npt01d2(:,:,2,1).*npt02d2c(:,:,2)).*...
               (1-npt00b(:,:,1,2).*npt01d2(:,:,1,2).*npt02d2c(:,:,1))).*...
            (1-(1-npt11e(:,:,1).*(1-npt01d1b(:,:,2,1).*npt11c(:,:,1))).*...
               (1-npt11e(:,:,2).*(1-npt01d1b(:,:,1,2).*npt11c(:,:,2))).*...
               (1-npt00b(:,:,2,1).*npt01d0b(:,:,2,1).*npt00b(:,:,1,2).*npt01d0b(:,:,1,2).*prod(npt02d2c,3).*npt11f.*(1-npt11d)).*...
               (1-npt00c(:,:,2,1).*npt01d0c(:,:,2,1).*npt00c(:,:,1,2).*npt01d0c(:,:,1,2).*prod(npt02d2d,3).*npt11f.*(1-npt12c)));
end
fc2 = fc2/3;

%% FC1
if nargout>1
    fc1 = (1 - 2*cij.*(ydiff(3).*(1-(1-y(3)).^t0) + ydiff(2).*(1-(1-y(2)).^t0) + ydiff(1).*(1-(1-y(1)).^t0))).*...
            (1 - prod(y(3)*(1-2*ydiff(3)).^ck + y(2).*(1-2*ydiff(2)).^ck + y(1).*(1-2*ydiff(1)).^ck,3))/3;
end
