function [x, k, alpha0] = ss_ngm_dirichlet_ubc_img(y, x, k,alpha0, pars)
% Do single-scale blind deconvolution using the input initializations
% x, k and alpha. 

[m, n] = size(y);
[k1, k2] = size(k);
khs1 = floor(k1/2);
khs2 = floor(k2/2);

xk_iter=pars.xk_iter;
img_pars=pars.img_pars;
img_pars.x0=x;

%using gradients for kernel estimation 
y2{1}=diff(y(khs1+1:end-khs1,khs2+1:end-khs2),1,1);
y2{2}=diff(y(khs1+1:end-khs1,khs2+1:end-khs2),1,2);

alpha=alpha0;
ker_opts=pars.kernel_pars;
lambda1 = img_pars.lambda1;
lambda_min =  img_pars.lambda_min;
if lambda1 <0.0005
    delta_lambda = 0.00005;
else
   delta_lambda = 0;
end
% lambda_min =img_pars.lambda_min;
% lambda_max = img_pars.lambda_max;
for i=1:xk_iter
  % fast x estimation with ngm and ubc
    img_pars.lambda1 = lambda1+delta_lambda*max((6-i),0);
    img_pars.lambda_min = lambda_min/lambda1*img_pars.lambda1;
    [x, x_full] = nbid_ngm_ubc_admm(y,k,img_pars);
    img_pars.x0=x_full;%using x_full as initial point for the next iteration
    x1{1}=diff(x,1,1);
    x1{2}=diff(x,1,2);
  % set up options for the kernel estimation
    ker_opts.alpha0=alpha;
    [alpha,fcost] = kernel_estimation_filter_space_fft(k,x1,y2,ker_opts);
    alpha0=reshape(alpha,k1,k2);
    if ker_opts.mode
        k=(alpha0-1)/(sum(alpha(:))-k1*k2);%update kernel using mode
    else
        k=alpha0/sum(alpha(:));%update kernel using expectation (default)
    end
%     Ccost = ker_opts.lambda_C*norm(conv2(k,C,'same'),'fro')^2;
%     Fcost = f+p+Ccost;
%     fprintf('Iteration=%d,F=%.4f,f=%.4f,p=%.4f,Ccost=%.4f\n',i,Fcost,f,p,Ccost);
     fprintf('Iteration=%d\n',i);
    % evolutions of kernel and latent sharp image
%     figure(1);
%     subplot(1,2,1);imagesc(k),colormap gray;   colorbar, drawnow;
%     subplot(1,2,2);imshow(x), drawnow;
%    figure(1);imagesc(k),colormap gray;   colorbar, drawnow;
%    warning off;figure(2), imshow(x), drawnow;
%   figure(1),imagesc(alpha),colorbar,drawnow;
    if i>=5
        max_k = 1;%max(abs(k(:)));
        r_k = max(abs(k(:)-k_old(:)))/max_k;
        %r_f = abs(fcost(end)-fcost(1))/abs(fcost(1));
        if length(fcost)<=2||r_k<=pars.k_tol%||r_f<=1e-4
            break;
        end
    end
    k_old = k;
end;
end