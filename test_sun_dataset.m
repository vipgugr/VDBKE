% test on Sun's data set
% Please download Sun's dataset via this link http://cs.brown.edu/~lbsun/deblur2013/deblur_iccp2013_testset_640.zip
% and unzip the dataset to the file deblur_iccp2013_testset_640.
clear all;
addpath ./deblur_iccp2013_testset_640
dir='./results';
opts.kernel_est_win =[]; 
opts.prescale = 1;
opts.gamma_correct = 1.0;
% non-blind settings
firls = [];

opts.use_ycbcr = 1;
%% Blind deconvolution

kernel_pars.lambda_C = 100;
kernel_pars.mode = 0;
opts.kernel_pars=kernel_pars;
% opts.ss_bid=1;
L=[0 0 0 8 0 6 4 4] + [19 17 15 27 13 21 23 23];

opts.firls.lambda=.0002;
lambda=0.0002;
opts.img_pars.lambda1=lambda;
opts.img_pars.lambda_min=5*lambda;%;0.01*lambda/0.001;
opts.img_pars.lambda_max=1;%*lambda/0.001;
% imgs = [5 27 41 55 56];
% kers = 1:8;
imgs = 5;kers=1;
for i=imgs
    for j=kers
    img=i;
    ker=j;
    opts.kernel_size=[L(ker),L(ker)];
    fn=sprintf('%d_%d_blurred.png',img,ker);
    B=im2double(imread(fn));
    opts.blur=B;
    tic
    [blur, kernel,deblur] = ms_ngm_dirichlet_ubc_img([], opts);
    toc
    kernel_name = sprintf('%s/kernel_%d_%d.mat',dir,i,j);
    save(kernel_name, 'kernel');
    imwrite(imresize(kernel/max(kernel(:)),[256 256],'nearest'),sprintf('%s/kernel_%d_%d.png',dir,i,j));
    imwrite(deblur,sprintf('%s/deblurred_%d_%d.png',dir,i,j));
    end
end


