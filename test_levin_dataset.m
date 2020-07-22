clear all
addpath ./test_data;
load  sizeL 
load errl.mat
dirname='results';
opts.kernel_est_win =[];% [100 100 355 355]; 
opts.gamma_correct = 1.0;
opts.kernel_init = 3; 
% non-blind settings
opts.use_ycbcr = 1;
fn=[];
NGMtime=0;

%% ngm+dirichlet+ubc
kernel_pars.lambda_C=.01; % regularization weight on the kernel prior term
% kernel_pars.Laplacian_filter=[0 -1 0;-1 4 -1;0 -1 0]/4;
kernel_pars.Laplacian_filter=[0 0 0;0 1 0;0 0 0];
opts.kernel_pars=kernel_pars;
L=[0 0 0 8 0 4 8 8];
lambda=0.00015; 
opts.img_pars.lambda1=lambda;  % regularization weight on the image prior term
opts.img_pars.lambda_min=5*lambda;
opts.img_pars.lambda_max=1;
ssdes=zeros(4,8);
ers=ssdes;
kernels=cell(4,8);
deblurs=cell(4,8);
etimes = ssdes;
for i=1:4
    for j=1:8
        aks=L(j);
        eval(sprintf('load test_data/im%02d_flit%02d',i+4,j))
        opts.blur=uint8(y*255); 
        ks = max(size(f))+aks;
        opts.kernel_size = [ks ks];
        tic
        [blur, kernel] = ms_ngm_dirichlet_ubc_img(fn, opts);
        etimes(i,j) = toc;
        deblur=deconvSps_undeterminedBC(y,kernel,0.002,70);
        deblur=double(uint8(deblur*255))/255;
%         figure(2),imagesc(rot90(kernel,4)),colormap gray;
%         figure(3),imshow(deblur);   
        [ssde]=comp_upto_shift(deblur,x);
        er=ssde/errL(i,j,1);
        fprintf('ssde=%f,er=%f\n',ssde,er); 
        ssdes(i,j)=ssde;
        ers(i,j)=er;
        deblurs{i,j}=deblur;
        kernels{i,j}=kernel;
        sclk=floor(min(size(deblur(:,:,1))./size(kernel))/3);
        k=kernel/max(kernel(:));
        k_sz=size(k);
        kex=deblur;
        for s=1:size(deblur,3)
         kex(1:k_sz(1)*sclk,1:k_sz(2)*sclk,s)=imresize(k,sclk,'nearest');
        end
        imwrite(kex,sprintf('%s/NGM_%d_%d_%.4f_%.4f.png',dirname,i,j,ssde,er));
    end
end
 save(sprintf('./%s/results.mat',dirname), 'ssdes', 'ers', 'kernels', 'deblurs','etimes');