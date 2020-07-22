function [x_fov, x, opt] = firls_deb_ubc(y,h,opt)
%fast IRLS for image deblurring with undetermined boundary conditions
% input kernel size must be odd
% This method is a faster version of the CG based method Zhou et al [1].
% It incorpates the techniques of mask decoupling [3] and FIRLS [2].
% References: 
% [1] X. Zhou, et al, A boundary condition based deconvolution framework for 
% image deblurring, J. Comput. Appl. Math. 261 (2014) 14¨C29.
% [2] X. Zhou, et al, Fast iteratively reweighted least squares for 
% lp regularized image deconvolution and reconstruction, IEEE ICIP, 2014,

% [3]  M. S. C. Almeida and M. A. T. Figueiredo; "Deconvolving images 
% with unknown boundaries using the alternating direction method of multipliers", 
% IEEE Transactions on Image Processing, vol. 22, No. 8, pp. 3074-3086, 2013.

% If you use this code to do NBID, please cite the above papers.
% Copyright (c) Xu Zhou <xuzhou@buaa.edu.cn>, Beihang University, 2015. 
%-------------------------------------------------------------------------
% Input: blurred image y, blur kernel h, and options opt
% Output: 
% x_fov : the deblurred image in the field of view, of the same size as y
% x: the full deblurred image.
% Please refer to [2] for more details on FIRLS and opt.

%------------------------------------------------------------------------- 

[M1, M2]=size(y);
[m1, m2]=size(h);
hks1=floor(m1/2);
hks2=floor(m2/2);
n1=M1+m1-1;
n2=M2+m2-1;
x=padarray(y,[hks1,hks2],'replicate','both');


% first and second order dirivative filters
dxf=[0 0 0;0 1 -1;0 0 0];
dyf=[0 0 0;0 1 0;0 -1 0];
dyyf=[0 -1 0;0 2 0;0 -1 0];
dxxf=[0 0 0;-1 2 -1;0 0 0];
dxyf=[0 0 0;0 1 -1;0 -1 1];

dxfr=rot90(dxf,2);
dyfr=rot90(dyf,2);
dxxfr=rot90(dxxf,2);
dyyfr=rot90(dyyf,2);
dxyfr=rot90(dxyf,2);

H=psf2otf(h,[n1 n2]);
Ht=conj(H);
Hx=psf2otf(dxf,[n1 n2]);
Hy=psf2otf(dyf,[n1 n2]);
Hxx=psf2otf(dxxf,[n1 n2]);
Hyy=psf2otf(dyyf,[n1 n2]);
Hxy=psf2otf(dxyf,[n1 n2]);

HH=H.*Ht;
HHx=Hx.*conj(Hx);
HHy=Hy.*conj(Hy);
HHxx=Hxx.*conj(Hxx);
HHyy=Hyy.*conj(Hyy);
HHxy=Hxy.*conj(Hxy);

% RR=(HHx+HHy);
RR=(HHx+HHy+HHxx+HHyy+HHxy);
%import parameters
lambda=opt.lambda;

w0=0.25;

% n_tol=opt.n_tol;
if isfield(opt,'alpha')
    alpha = opt.alpha ;
else
    alpha = 2/3;
end
if isfield(opt,'beta_a')
    beta_a = opt.beta_a ;
else
    beta_a = lambda*alpha*(20/255)^(alpha-2);
end
if isfield(opt,'lambda_u')
    lambda_u = opt.lambda_u ;
else
    lambda_u=min(0.1,5000*lambda);
end

if isfield(opt,'isnr_display')
else
    opt.isnr_display = 0;
end
if opt.isnr_display==1
    I=opt.groundtruth;
end
if isfield(opt,'cost_display')
else
    opt.cost_display = 0;
end
if isfield(opt,'inner_iter')
    N2=opt.inner_iter;
else
    N2=4;
end
if isfield(opt,'out_iter')
    N1=opt.out_iter;
else
    N1=5;
end
if isfield(opt,'epsilon')
    epsilon=opt.epsilon;
else
    epsilon=0.01;
end
c=alpha*lambda;
beta=alpha*lambda/epsilon^(2-alpha);

iter=0;
xpad=padarray(x,[1 1],'circular','both');
dx=conv2(xpad,dxf,'valid');
dy=conv2(xpad,dyf,'valid');
dxx=conv2(xpad,dxxf,'valid');
dyy=conv2(xpad,dyyf,'valid');
dxy=conv2(xpad,dxyf,'valid');

adx=abs(dx);
ady=abs(dy); 
adxx=abs(dxx);
adyy=abs(dyy);
adxy=abs(dxy);
% initialization for dual variables 
du=zeros(n1,n2);
dvx=du;
dvy=dvx;
dvxx=dvx;
dvyy=dvx;
dvxy=dvx;

X=fft2(x);
Ax=real(ifft2(H.*X));
invA=(HH+beta_a/lambda_u*RR);
totiter=0;

% Wx=ones(size(adx))*lambda;
% Wy=ones(size(ady))*lambda;
% Wxx=ones(size(adxx))*lambda*w0;
% Wyy=ones(size(adyy))*lambda*w0;
% Wxy=ones(size(adxy))*lambda*w0;


while iter<N1
    iter=iter+1;
%     x_old=x;
%     epsilon=(c/beta)^(1/(2-alpha));
    Wx=min(beta,c*adx.^(alpha-2));
    Wy=min(beta,c*ady.^(alpha-2));
    Wxx=min(beta,c*adxx.^(alpha-2))*w0;
    Wyy=min(beta,c*adyy.^(alpha-2))*w0;
    Wxy=min(beta,c*adxy.^(alpha-2))*w0;
    %admm
    i=0;
    while i<N2
        i=i+1;
        totiter=totiter+1;
        %u subproblem       
        u=Ax+du;
        u(hks1+1:end-hks1,hks2+1:end-hks2)=(y+lambda_u*u(hks1+1:end-hks1,hks2+1:end-hks2))/(1+lambda_u);
        
        %v subproblem
        vx=beta_a*(dx+dvx)./(Wx+beta_a);
        vy=beta_a*(dy+dvy)./(Wy+beta_a);
        vxx=beta_a*(dxx+dvxx)./(Wxx+beta_a);
        vyy=beta_a*(dyy+dvyy)./(Wyy+beta_a);
        vxy=beta_a*(dxy+dvxy)./(Wxy+beta_a);
        %update dual varible
        du=du-u+Ax;
        dvx=dvx-vx+dx;
        dvy=dvy-vy+dy;
        dvxx=dvxx-vxx+dxx;
        dvyy=dvyy-vyy+dyy;
        dvxy=dvxy-vxy+dxy;
        %x subproblem
        Y=fft2(u-du).*Ht;
        
        tempx=vx-dvx;
        tempy=vy-dvy;
        tempxx=vxx-dvxx;
        tempyy=vyy-dvyy;
        tempxy=vxy-dvxy;
        
%         tempx=-[tempx(:,1)-tempx(:,n2) diff(tempx,1,2)];
%         tempy=-[tempy(1,:)-tempy(n1,:); diff(tempy,1,1)];
        tempx=conv2(padarray(tempx,[1 1],'circular','both'),dxfr,'valid');
        tempy=conv2(padarray(tempy,[1 1],'circular','both'),dyfr,'valid');
        tempxx=conv2(padarray(tempxx,[1 1],'circular','both'),dxxfr,'valid');
        tempyy=conv2(padarray(tempyy,[1 1],'circular','both'),dyyfr,'valid');
        tempxy=conv2(padarray(tempxy,[1 1],'circular','both'),dxyfr,'valid');
        
        X=Y+beta_a/lambda_u*fft2(tempx+tempy+tempxx+tempyy+tempxy);
        X=X./invA;
        Ax=real(ifft2(H.*X));
        x=real(ifft2(X));
        
        %calculate current gradient image
%         dx=[diff(x,1,2),x(:,1)-x(:,n2)];
%         dy=[diff(x,1,1);x(1,:)-x(n1,:)];
        xpad=padarray(x,[1 1],'circular','both');
        dx=conv2(xpad,dxf,'valid');
        dy=conv2(xpad,dyf,'valid');
        dxx=conv2(xpad,dxxf,'valid');
        dyy=conv2(xpad,dyyf,'valid');
        dxy=conv2(xpad,dxyf,'valid');
        adx=abs(dx);
        ady=abs(dy); 
        adxx=abs(dxx);
        adyy=abs(dyy);
        adxy=abs(dxy);
    end
    
    if opt.cost_display==1
        r=Ax(hks1+1:end-hks1,hks2+1:end-hks2)-y;
        opt.cost1(iter)=0.5*norm(r,'fro')^2;
        opt.cost2(iter)=lambda*sum(sum(adx.^(alpha)+ady.^(alpha)+...
            adxx.^(alpha)+adyy.^(alpha)+adxy.^(alpha)));
        opt.cost3(iter)=opt.cost1(iter)+opt.cost2(iter);
        fprintf('Outiter=%d,costf=%f,',iter,opt.cost3(iter));
    end
    if opt.isnr_display==1
        opt.isnr(iter)=20*log10(norm(y-I,'fro')/norm(x(hks1+1:end-hks1,hks2+1:end-hks2)-I,'fro'));
        fprintf('isnr=%f,beta=%f\n',opt.isnr(iter),beta);
    else
        fprintf('beta=%f\n',beta);
    end    
%     e(iter)=max(abs(x(:)-x_old(:)));
end
x_fov=x(hks1+1:end-hks1,hks2+1:end-hks2);
end

