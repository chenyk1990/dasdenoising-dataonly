function [dn,ds] = das_pwsmooth_lop_mf(dip,w1,n1,n2,ns,order,eps,ndn,nds,type_mf,ifsmooth,dn,ds)
% str_pwsmooth_lop: plane-wave smoothing
%
% INPUT:
% dn: model   noisy data
% dip: slope (2D array)
% ns:       spray radius
% order:    PWD order
% eps: regularization (default:0.01);

% OUTPUT:
% ds:  smoothed data
%
% Copyright (c) 2021 by the Society of Exploration Geophysicists.
% You must read and accept usage terms at:
%    https://software.seg.org/disclaimer.txt before use.

if ndn~=nds
    error('Wrong size %d != %d',ndn,nds);
end

if ns~=1
    fprintf('nds=%d, n1=%d, n2=%d\n',nds,n1,n2);
    ds=zeros(nds,1);
%     [ dn,ds ] = yc_adjnull( adj,add,ndn,nds,dn,ds );
    
    ns2=2*ns+1;%spray diameter
    n12=n1*n2;
    
    u=zeros(n1,ns2,n2);
    utmp=zeros(n12*ns2,1);
    w=zeros(ns2,1);
    % w1=zeros(n1,n2);
    
    for is=0:ns2-1
        w(is+1)=ns+1-abs(is-ns);
    end
    
    
    % for Normalization
    % t=zeros(n12,1);

%     [dn,utmp]=pwspray_lop(0,0,n12,n12*ns2,dn,utmp,dip,ns,n1,n2,order,eps);
    [utmp] = str_pwspray_lop2d(dn,dip,ns,order,eps);
    size(utmp)
    u=reshape(utmp,n1,ns2,n2);
    
    for i2=1:n2
        if type_mf==0   %MF
            u(:,:,i2)= das_mf(u(:,:,i2),ns2,1,2);
        else            %SVMF
            u(:,:,i2)= das_svmf(u(:,:,i2),ns2,1,2);
        end
    end
    
    fprintf('size of u(:,:,i2)\n');
    [n1,n2,ns]
    if ifsmooth %with smoothing or only with median filtering
        for i2=0:n2-1
            for i1=0:n1-1
                ws=w1(i1+1,i2+1);
                for is=0:ns2-1
                    ds(i2*n1+i1+1)=ds(i2*n1+i1+1)+u(i1+1,is+1,i2+1)*w(is+1)*ws;
                end
            end
        end
    else
        for i2=0:n2-1
            for i1=0:n1-1
                ds(i2*n1+i1+1)=ds(i2*n1+i1+1)+u(i1+1,ns+1,i2+1);
            end
        end
        
    end
    
    
else
    ds=dn;
end


return




