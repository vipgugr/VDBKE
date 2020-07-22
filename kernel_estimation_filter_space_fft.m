function [alpha,fx,stepsize] = kernel_estimation_filter_space_fft(k,x,y,opt)
%Gradient projection method for dirichlet minimization
%   fft version
%   Detailed explanation goes here
lambda=opt.lambda;
lambda_C=opt.lambda_C;
C=opt.Laplacian_filter;
max_iter=opt.max_iter;
ba=opt.back_alpha;
bb=opt.back_beta;
lb=opt.lower_bound;
ng_min=opt.ng_min;
cost_display=opt.cost_display;
[ks1, ks2]=size(k);
alpha=opt.alpha0;
% compute A, Ad and b
[Ad,b,xtAx,Ax] = dirichlet_Adbc_fft(x,y,ks1,ks2,lambda_C,C);
iter=0;
% compute cost
[fx(iter+1),atAa,Adta,bta,Sa,den1,Xalpha] = dirichlet_cost_by_fft(alpha,lambda,xtAx,Ad,b);
costcalls=1;
if cost_display
fprintf('iteration=%d,cost=%f\n',iter,fx(iter+1));
end
stepsize=Sa;%min(Sa,2000);
while iter<max_iter
%     if ng<ng_min
%         break;
%     end
    % compute gradient
    den2=den1^2/(4*Sa+2);
    g=lambda*(alpha-1).*(psi(1,alpha)-psi(1,Sa));
    Aa=Ax(Xalpha);
    g=g+(2*Aa+Ad)/den1+b/(2*Sa)-(atAa+Adta)/den2-bta/(2*Sa^2);    
    %backtracking for step size
    d=-g;
%     tmax=Sa;
%     tmax=(Sa+stepsize)/2;
    tmax=min(Sa,stepsize*1.2);% maximum step size
%     fprintf('tmax=%f,',tmax);
    t=tmax;
    temp=alpha+t*d;
    temp=max(temp,lb);
    [ftemp,atAa,Adta,bta,Sa,den1,Xalpha] = dirichlet_cost_by_fft(temp,lambda,xtAx,Ad,b);
    costcalls=costcalls+1;   
    dg=sum(sum((temp-alpha).*g));
    fx1=fx(iter+1)+ba*dg;
    
    while ftemp>fx1
        t=bb*t;
        temp=alpha+t*d;
        temp=max(temp,lb);
        %call cost
        [ftemp,atAa,Adta,bta,Sa,den1,Xalpha] = dirichlet_cost_by_fft(temp,lambda,xtAx,Ad,b);
        costcalls=costcalls+1;
        dg=sum(sum((temp-alpha).*g));
        fx1=fx(iter+1)+ba*dg; 
%         if t<10^(-8)
%             break;
%         end
    end
    iter=iter+1;
    stepsize=t;
    alpha=temp;
    rf=abs((ftemp-fx(iter))/ftemp);
    fx(iter+1)=ftemp;

    if cost_display        
    fprintf('iteration=%d,costcalls=%d,cost=%f,rf=%f,step_size=%f\n',iter,costcalls,fx(iter+1),rf,t);
    end
    if t<10^(-2)||rf<ng_min%max(abs(g(:)))<ng_min
        break;
    end  
end
fprintf('DIter=%d,costcalls=%d,cost=%.4f,cost=%.4f,rf=%f,Sa=%.0f\n',iter,costcalls,fx(1),fx(iter+1),rf,Sa);
% fprintf('Iter=%d,costcalls=%d,cost=%.4f,cost=%.4f,rf=%f,|g|=%.6f,t=%.1f\n',...
%     iter,costcalls,fx(1),fx(iter+1),rf,max(abs(g(:))),t);
end
function [f,atAa,Adta,bta,Sa,den1,Xalpha] = dirichlet_cost_by_fft(alpha,lambda,xtAx,Ad,b)
%compute cost by convolution   
Sa=sum(alpha(:));
den1=2*Sa*(Sa+1);
[atAa, Xalpha]=xtAx(alpha);
Adta=sum(sum(Ad.*alpha));
bta=sum(sum(b.*alpha));
temp=(alpha-1).*(psi(alpha)-psi(Sa));
f=lambda*(sum(temp(:))+gammaln(Sa)-sum(gammaln(alpha(:))))+(atAa+Adta)/den1+bta/(2*Sa);
end
