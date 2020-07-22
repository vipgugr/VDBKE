clear all;
close all;
lambda = 0.0002;
cont= 1;
lambda2 = 1;
img_pars.lambda_u=.1;
img_pars.lambda_max = lambda2;
img_pars.lambda1 = lambda;
img_pars.xv_iter=1;
img_pars.N2=2;
img_pars.N1=100;
img_pars.IF=sqrt(2);
img_pars.cost_display=1;
N = img_pars.N1;
fontsize = 14;
datadir='deblur_iccp2013_testset_640';
dir = 'plots';
load PSF
%%
if cont ==1
    lambda1 = 5*lambda;
else
    lambda1 = lambda2;
end
img_pars.lambda_min = lambda1;    
imgs = 1:80;
kers = 1:8;
imgs = 5;
kers = 1;
PCOST = cell(length(imgs),length(kers));
FCOST = cell(length(imgs),length(kers));
HATF = cell(length(imgs),length(kers));
for i = imgs
    for j=kers    
        B = imread(sprintf('./%s/%d_%d_blurred.png',datadir,i,j));
        B = double(B)/255;
        img_pars.x0 = B;
        k = PSFs{j};
        k = rot90(k,2);
        tic
        [x_fov, x, pcost, fcost,hatf] = nbid_ngm_ubc_admm(B,k,img_pars);
        toc
        FCOST{i,j}=fcost;
        PCOST{i,j}=pcost;
        HATF{i,j}=hatf;
        imwrite(x,sprintf('./%s/deblurred_%d_%d_%d_%d.png',dir,cont,lambda*10000,i,j));
%         figure(1),plot(1:N,fcost,'.-r',1:N,pcost,'-b',1:N,hatf,'*-k');
        figure(1),plot(1:N,fcost,'.-r',1:N,pcost,'-b');
        xlabel('Iteration','fontsize',fontsize);ylabel('Cost','fontsize',fontsize);
        title(sprintf('image %d kernel %d',i,j),'fontsize',fontsize);
        % set(gca,'Fontsize',fontsize)
%         h_legend = legend('$F(x)$','$R_x(x)$','$\hat F(x)$');
        h_legend = legend('$F(x)$','$R_x(x)$');
        set(h_legend,'Interpreter','latex','FontSize',fontsize,'Location','Best');
        % subplot(2,1,1),plot(pcost);xlabel('Iteration');ylabel('NGM(x)');
        % subplot(2,1,2),plot(fcost);xlabel('Iteration');ylabel('$\frac{1}{2}\|Hx-y\|_2^2+0.0002*NGM(x)$','Interpreter','latex');
%         eval(sprintf('print -dpng ./%s/costngm_%d_%d_%d_%d_%d_%d.png',dir,cont,i,j,lambda*10000,lambda2));
        eval(sprintf('print -dpng ./%s/cost_%d_%d_%d_%d.png',dir,cont,lambda*10000,i,j));
    end
end
save(sprintf('./%s/PCOST_%d_%d.mat',dir,cont,lambda*10000),'PCOST');
save(sprintf('./%s/FCOST_%d_%d.mat',dir,cont,lambda*10000),'FCOST');
save(sprintf('./%s/HATF_%d_%d.mat',dir,cont,lambda*10000),'HATF');
%%


