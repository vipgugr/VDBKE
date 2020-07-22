function [y]=deconvSps_undeterminedBC(x,h,we,max_it)
%note: size(h) is expected to be odd in both dimensions 
%undetermined BC
if (~exist('max_it','var'))
   max_it  =200;
end

[n,m]=size(x);
[M1,M2]=size(h);
n=n+M1-1;m=m+M2-1;
dxf=[1 -1];
dyf=[1;-1];
dyyf=[-1; 2; -1];
dxxf=[-1, 2, -1];
dxyf=[-1 1;1 -1];

weight_x=ones(n,m-1);
weight_y=ones(n-1,m);
weight_xx=ones(n,m-2);
weight_yy=ones(n-2,m);
weight_xy=ones(n-1,m-1);

[y]=deconvL2_w_undeterminedBC(x,h,we,max_it,weight_x,weight_y,weight_xx,weight_yy,weight_xy );
if we>0
    w0=0.1;
    exp_a=0.8;
    thr_e=0.01; 

    for t=1:2

    dy=conv2(y,rot90(dyf,2),'valid');
    dx=conv2(y,rot90(dxf,2),'valid');
    dyy=conv2(y,rot90(dyyf,2),'valid');
    dxx=conv2(y,rot90(dxxf,2),'valid');
    dxy=conv2(y,rot90(dxyf,2),'valid');


    weight_x=w0*max(abs(dx),thr_e).^(exp_a-2); 
    weight_y=w0*max(abs(dy),thr_e).^(exp_a-2);
    weight_xx=0.25*w0*max(abs(dxx),thr_e).^(exp_a-2); 
    weight_yy=0.25*w0*max(abs(dyy),thr_e).^(exp_a-2);
    weight_xy=0.25*w0*max(abs(dxy),thr_e).^(exp_a-2);

    [y]=deconvL2_w_undeterminedBC(x,h,we,max_it,weight_x,weight_y,weight_xx,weight_yy,weight_xy );
    end
end
m1=floor((M1-1)/2);
m2=floor((M2-1)/2);
y=y(m1+1:n-m1,m2+1:m-m2);
end


