clear all;
addpath ./images
dir='./images';

opts.kernel_est_win =[]; 
% set initial downsampling size for really large images
opts.prescale = 1;
% set this to 1 for no gamma correction - default 1.0
opts.gamma_correct = 1.0;
opts.kernel_init = 3; 


% non-blind settings
firls = [];
firls.lambda=.0002;
firls.alpha=2/3;
opts.firls = firls;
opts.use_ycbcr = 0;

fn=[];
% fn='4_input.jpg'; ks = [51 51];
% fn='Blurry3_8.tif';  ks = [37 125];
% fn='y1.tif';

%% motion blur + noise
fn='4_input.jpg'; ks = [51 51];

% non-blind settings
firls = [];
firls.lambda=.0004;
firls.alpha=2/3;
opts.firls = firls;
opts.use_ycbcr = 0;

lambda=0.0005; % use a large lambda to handle the noise
kernel_pars.lambda_C=100;
opts.kernel_pars=kernel_pars;
opts.kernel_size = ks;
opts.img_pars.lambda1=lambda;
opts.img_pars.lambda_min=5*lambda;%;
opts.img_pars.lambda_max=1;%
tic
[blur, kernel,deblur] = ms_ngm_dirichlet_ubc_img(fn, opts);
toc
etime=toc;
figure,imagesc(kernel)%,colormap gray;
figure,imshow(deblur);

% Save results
fname=fn(1:end-4);
kernelname=sprintf('%s/%s_k_%f_%.2f.png',dir,fname,lambda,kernel_pars.lambda_C);
filename=sprintf('%s/%s_%f_%.2f.png',dir,fname,lambda,kernel_pars.lambda_C);
k=kernel/max(kernel(:));
imwrite(uint8(deblur*255),filename);
imwrite(uint8(k*255),kernelname);

k_sz=size(k);
sclk=min(size(deblur(:,:,1))./opts.kernel_size/3);
% sclk=100/k_sz;
sk=round(k_sz*sclk);
kex=deblur;
for i=1:size(deblur,3)
 kex(1:sk(1),1:sk(2),i)=imresize(k,sk,'nearest');
end
imwrite(uint8(255*kex),sprintf('%s/%s_kdeb_%f_%.2f.png',dir,fname,lambda,kernel_pars.lambda_C));


%% motion blur (large)
fn='Blurry3_8.png';  ks = [125 37];
% non-blind settings
firls = [];
firls.lambda=.0002;
firls.alpha=2/3;
opts.firls = firls;
opts.use_ycbcr = 0;


kernel_pars.lambda_C=100;
kernel_par.lower_bound = 0.1; % use a small lb for large blur kernel
opts.kernel_pars=kernel_pars;
opts.kernel_size = ks;
lambda=0.0002;
opts.img_pars.lambda1=lambda;
opts.img_pars.lambda_min=5*lambda;%;
opts.img_pars.lambda_max=1;%
tic
[blur, kernel,deblur] = ms_ngm_dirichlet_ubc_img(fn, opts);
toc
etime=toc;
figure,imagesc(kernel)%,colormap gray;
figure,imshow(deblur);

% Save results
fname=fn(1:end-4);
kernelname=sprintf('%s/%s_k_%f_%.2f.png',dir,fname,lambda,kernel_pars.lambda_C);
filename=sprintf('%s/%s_%f_%.2f.png',dir,fname,lambda,kernel_pars.lambda_C);
k=kernel/max(kernel(:));
imwrite(uint8(deblur*255),filename);
imwrite(uint8(k*255),kernelname);

k_sz=size(k);
sclk=min(size(deblur(:,:,1))./opts.kernel_size/3);
% sclk=100/k_sz;
sk=round(k_sz*sclk);
kex=deblur;
for i=1:size(deblur,3)
 kex(1:sk(1),1:sk(2),i)=imresize(k,sk,'nearest');
end
imwrite(uint8(255*kex),sprintf('%s/%s_kdeb_%f_%.2f.png',dir,fname,lambda,kernel_pars.lambda_C));

%% atmostpheric blur
fn='y1.tif'; 
% fn='y7.tif';
% fn='y9.tif';
ks = [15 15];
% non-blind settings
firls = [];
firls.lambda=.0002;
firls.alpha=2/3;
opts.firls = firls;
opts.use_ycbcr = 0;


kernel_pars.lambda_C=10;
opts.kernel_pars=kernel_pars;
opts.kernel_size = ks;
lambda=0.0002;
etimes=zeros(length(lambda));
opts.img_pars.lambda1=lambda;
opts.img_pars.lambda_min=5*lambda;%;
opts.img_pars.lambda_max=1;%
tic
[blur, kernel,deblur,option] = ms_ngm_dirichlet_ubc_img(fn, opts);
toc
etime=toc;
figure,imagesc(kernel)%,colormap gray;
figure,imshow(deblur);

% %% NBID using zoran's method, EPLL
% % you need to download the zoran's code from  http://people.csail.mit.edu/danielzoran/epllcode.zip
% noiseSD = 0.005;
% patchSize = 8;
% addpath ./epllcode
% ks = floor(size(kernel, 1)/2);
% B = im2double(imread(fn));
% y = padarray(B, [1 1]*ks, 'replicate', 'both');
% for a=1:4
%   y = edgetaper(y, kernel);
% end
% load GSModel_8x8_200_2M_noDC_zeromean.mat
% excludeList = [];
% prior = @(Z,patchSize,noiseSD,imsize) aprxMAPGMM(Z,patchSize,noiseSD,imsize,GS,excludeList);
% 
% % comment this line if you want the total cost calculated
% LogLFunc = [];
% 
% tic
% [cleanI] = EPLLhalfQuadraticSplitDeblur(y,64/noiseSD^2,kernel,patchSize,50*[1 2 4 8 16 32 64],1,prior,y,LogLFunc);
% toc
% 
% deblur = cleanI(ks+1:end-ks,ks+1:end-ks);

% Save results 
fname=fn(1:end-4);
kernelname=sprintf('%s/%s_k_%f_%.2f.png',dir,fname,lambda,kernel_pars.lambda_C);
filename=sprintf('%s/%s_%f_%.2f.png',dir,fname,lambda,kernel_pars.lambda_C);
k=kernel/max(kernel(:));
imwrite(uint8(deblur*255),filename);
imwrite(uint8(k*255),kernelname);

k_sz=size(k);
sclk=min(size(deblur(:,:,1))./opts.kernel_size/3);
% sclk=100/k_sz;
sk=round(k_sz*sclk);
kex=deblur;
for i=1:size(deblur,3)
 kex(1:sk(1),1:sk(2),i)=imresize(k,sk,'nearest');
end
imwrite(uint8(255*kex),sprintf('%s/%s_kdeb_%f_%.2f.png',dir,fname,lambda,kernel_pars.lambda_C));