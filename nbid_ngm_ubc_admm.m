function [x_fov, x, pcost, fcost,hatf] = nbid_ngm_ubc_admm(B,k,pars)
%Non-blind image deconvolution using NGM and undetermined boundary
%condition. This is the implementation of Alg.1 
%   Input:Blurred image B, Kernel k, NGM weight lambda, Number
%   of iterations N1(outer loop) and N2 (inner loop) 
%   Output:Deblurred image x

lambda1=pars.lambda1; % weight on the NGM term

if isfield(pars,'lambda_min')
    lambda_min=pars.lambda_min; % weight the constraint term
else
    lambda_min=lambda1*5;
end
if isfield(pars,'lambda_max')
    lambda_max=pars.lambda_max; % weight the constraint term
else
    lambda_max=1;
end
if isfield(pars,'cost_display')
    cost_display=pars.cost_display; 
else
    cost_display=0;
end
if isfield(pars,'IF')
    IF=pars.IF; % Increase Factor
else
    IF=sqrt(2);
end
if isfield(pars,'N2')
    N2=pars.N2; % inner loop iterations
else
    N2=2;
end
if isfield(pars,'N1')
    N1=pars.N1; % outer loop iterations
else
    N1=20;
end
if isfield(pars,'lambda_u')
    lambda_u=pars.lambda_u;%penalty weight on the boundary constraint
else
    lambda_u = 0.1;
end
if isfield(pars,'xv_iter')
    xv_iter =pars.xv_iter;%penalty weight on the boundary constraint
else
    xv_iter = 1;
end
lambda_v=lambda_min;


dx = [1 -1];
dy = dx';
[m,n]=size(B);
hks1=floor(size(k,1)/2);
hks2=floor(size(k,2)/2);
M=m+2*hks1;
N=n+2*hks2;

K=psf2otf(k,[M,N]);
Kt=conj(K);

Dx=psf2otf(dx,[M,N]);
Dy=psf2otf(dy,[M,N]);
Dxt=conj(Dx);
Dyt=conj(Dy);
DtD=Dxt.*Dx+Dyt.*Dy;
KtK=Kt.*K;

x=pars.x0; % initial latent image
if size(x)==size(B)
    x=padarray(x,[hks1 hks2],'replicate','both'); % pads x with replicate boundary
end

MtB=padarray(B,[hks1 hks2],0,'both'); % pads B with zero boundary
MtM=padarray(ones(m,n),[hks1 hks2],0,'both');

x1=[diff(x, 1, 2), x(:,1) - x(:,end)]; 
x2=[diff(x, 1, 1); x(1,:) - x(end,:)];

u=padarray(B,[hks1 hks2],'replicate','both');
du=zeros(size(u));% dual variable for u;

fcost=[];
pcost=fcost;
clear Dx Dy Dxt Dyt;
i=1;
X = fft2(x);

while i<=N1%lambda_v<=lambda_max

    %update u
    Ax = real(ifft2(X.*K));

%     Ax=conv2(padarray(x,[hks1 hks2],'circular','both'),k,'valid');
    u=(MtB+lambda_u*(Ax+du))./(MtM+lambda_u);
    % update dual variable
    du=du+Ax-u;
%     Ktu=conv2(padarray(u-du,[hks1 hks2],'circular','both'),kt,'valid');
    Ktu = fft2(u-du).*Kt;
    
    invA=1./(KtK+lambda_v/lambda_u*DtD);
    lambda=lambda1/lambda_v;

    %update v
    for j=1:xv_iter
        temp1=abs(x1);
        temp2=abs(x2);
        v1=temp1;
        v2=temp2;
        s1=sign(x1);
        s2=sign(x2);   
        for t=1:N2
            G1=mean(v1(:));
            G2=mean(v2(:));
            beta1=lambda*G1./(v1+G1).^2;
            beta2=lambda*G2./(v2+G2).^2;
            v1=max(temp1-beta1,0);
            v2=max(temp2-beta2,0);
        end
        v1=s1.*v1;
        v2=s2.*v2;
    end
%   Update x
    temp1=-[v1(:,1)-v1(:,end), diff(v1,1,2)];
    temp2=-[v2(1,:)-v2(end,:); diff(v2,1,1)];
    X= Ktu + lambda_v/lambda_u*fft2((temp1+temp2));
    X=invA.*X;
    x=real(ifft2(X));
    x1=[diff(x, 1, 2), x(:,1) - x(:,end)]; 
    x2=[diff(x, 1, 1); x(1,:) - x(end,:)];
    

   %pcost at iteration i
   if cost_display
    temp1=abs(x1(:));
    temp2=abs(x2(:));
    G1=mean(temp1(:));
    G2=mean(temp2(:));
    pcost(i)=lambda1*(sum(temp1./(temp1+G1))+sum(temp2./(temp2+G2)));
    fcost(i)=0.5*norm(Ax(hks1+1:m+hks1,hks2+1:n+hks2)-B,'fro')^2+pcost(i);
    hatf(i)=fcost(i)+lambda_v*(norm(v1-x1,'fro')^2+norm(v2-x2,'fro')^2);
    fprintf('Outer iteration %d: fcost=%.6f hatf=%.6f pcost=%.6f \n',i,fcost(i),hatf(i),pcost(i)); 
   end
   lambda_v=min(lambda_v*IF,lambda_max);
   i=i+1;
end
x_fov=x(hks1+1:m+hks1,hks2+1:n+hks2);
if nargout == 3
    temp1=abs(x1(:));
    temp2=abs(x2(:));
    G1=mean(temp1(:));
    G2=mean(temp2(:));
    pcost(1)=lambda1*(sum(temp1./(temp1+G1))+sum(temp2./(temp2+G2)));
    pcost(2)=lambda_v*(norm(x1-v1,'fro')^2 +norm(x2-v2,'fro')^2);
end
%     fcost=0.5*norm(Ax(hks1+1:m+hks1,hks2+1:n+hks2)-B,'fro')^2+pcost;
%     fprintf('NGM, lambda_v=%f: fcost=%.4f pcost=%.4f\n',lambda_v,fcost(i),pcost(i)); 
%     fprintf('Outer iteration pcost=%.4f\n',pcost(i)); 

end

