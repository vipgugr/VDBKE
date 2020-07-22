function [y] = deconvL2_w_undeterminedBC(x,h,we,max_it,weighr_x,weighr_y,weighr_xx,weighr_yy,weighr_xy )
%Undetermined BC
if (~exist('max_it','var'))
   max_it=200;
end
[n,m]=size(x);
[M1,M2]=size(h);
m1=floor((M1-1)/2);
m2=floor((M2-1)/2);
hr=rot90(h,2);
K=35;
b=conv2(x,hr);
y=x([ones(1,m1),1:end,end*ones(1,m1)],[ones(1,m2),1:end,end*ones(1,m2)]);

n=n+M1-1;
m=m+M2-1;
if M1>K||M2>K
    Hr=fft2(hr,n,m);
    Hrc=conj(Hr);
end

if (~exist('weighr_x','var'))
  weighr_x=ones(n,m-1);
  weighr_y=ones(n-1,m);
  weighr_xx=zeros(n,m-2);
  weighr_yy=zeros(n-2,m);
  weighr_xy=zeros(n-1,m-1);
end
dxf=[1 -1];
dyf=[1;-1];
dyyf=[-1; 2; -1];
dxxf=[-1, 2, -1];
dxyf=[-1 1;1 -1];

if M1>K||M2>K
    temp=real(ifft2(fft2(y).*Hrc));
    temp(n-M1+2:n,:)=0;
    temp(1:n-M1+1,m-M2+2:m)=0;
    Ax=real(ifft2(fft2(temp).*Hr));
else
    Ax=conv2(conv2(y,h,'valid'),hr);
end

if we>0
Ax=Ax+we*conv2(weighr_x.*conv2(y,rot90(dxf,2),'valid'),dxf);
Ax=Ax+we*conv2(weighr_y.*conv2(y,rot90(dyf,2),'valid'),dyf);
Ax=Ax+we*conv2(weighr_xx.*conv2(y,rot90(dxxf,2),'valid'),dxxf);
Ax=Ax+we*conv2(weighr_yy.*conv2(y,rot90(dyyf,2),'valid'),dyyf);
Ax=Ax+we*conv2(weighr_xy.*conv2(y,rot90(dxyf,2),'valid'),dxyf);
end

r = b - Ax;

for iter = 1:max_it  
     rho = (r(:)'*r(:));

     if ( iter > 1 )                       % direction vector
        beta = rho / rho_1;
        p = r + beta*p;
     else
        p = r;
     end
    if M1>K||M2>K
        temp=real(ifft2(fft2(p).*Hrc));
        temp(n-M1+2:n,:)=0;
        temp(1:n-M1+1,m-M2+2:m)=0;
        Ap=real(ifft2(fft2(temp).*Hr));
    else
        Ap=conv2(conv2(p,h,'valid'),hr);
    end
    if we>0
     Ap=Ap+we*conv2(weighr_x.*conv2(p,rot90(dxf,2),'valid'),dxf);
     Ap=Ap+we*conv2(weighr_y.*conv2(p,rot90(dyf,2),'valid'),dyf);
     Ap=Ap+we*(conv2(weighr_xx.*conv2(p,rot90(dxxf,2),'valid'),dxxf));
     Ap=Ap+we*(conv2(weighr_yy.*conv2(p,rot90(dyyf,2),'valid'),dyyf));
     Ap=Ap+we*(conv2(weighr_xy.*conv2(p,rot90(dxyf,2),'valid'),dxyf));
    end

     q = Ap;
     alpha = rho / (p(:)'*q(:) );
     y = y + alpha * p;                    % update approximation vector

     r = r - alpha*q;                      % compute residual

     rho_1 = rho;
end
disp(rho);
end
