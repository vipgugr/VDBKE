This code implements the method in 
X. Zhou, J. Mateos, F. Zhou, R. Molina, and A. K. Katsaggelos, 
"Variational Dirichlet Blur Kernel Estimation," 
to appear in IEEE Transactions on Image Processing, 2015.

If you're using this code, please cite this paper.

Author: Xu Zhou (xuzhou@buaa.edu.cn), Beihang University.
Please send any questions/comments/suggestions to <xuzhou@buaa.edu.cn>. 
------------------------------------------------------------------------------

The flow of the code is

Kernel Estimation: 
ms_ngm_dirichlet_ubc_img.m
    -> ss_ngm_dirichlet_ubc_img.m 
        -> nbid_ngm_ubc_admm.m (latent image estimation )
            -> kernel_estimation_filter_space_fft.m (VD kernel estimation)
NBID:
-> firls_deb_ubc.m: Final image reconstruction

ms_ngm_dirichlet_ubc_img performs multicale BID.
ss_ngm_dirichlet_ubc_img does single scale BID, which alternatively 
estimates the image and kernel. 
After the kernel estimation, we then run firls_deb_ubc.m (NBID) to estimate
the final image.

Non-blind deconvolution: firls_deb_ubc.m 
This is a faster version of the method in Zhou et al. JCAM 2014. 

------------------------------------------------------------------------------

There are 4 scripts for how to run these functions to reproduce the results 
in the paper. 

test_levin_dataset.m: 
This script reproduces the results of our method on Levin's dataset.

test_real_data.m:
This script reproduces the results on real data in our paper. 

test_sun_dataset.m
This script reproduces the results on Sun's dataset 640.
Please download the dataset from his project page.

test_on_convergence_sundataset.m
This script tests the convergence of Alg.1 presented in the paper.

------------------------------------------------------------------------------
The code is reimplemented for readability and it might not give the exact same
results as the paper above.
------------------------------------------------------------------------------

-------------------------------------------
Copyright and disclaimer
-------------------------------------------

(c) 2015 X. Zhou, J. Mateos, F. Zhou, R. Molina, and A.K. Katsaggelos.

The programs are granted free of charge for research and education 
purposes only. Scientific results produced using the software provided 
shall acknowledge the use of the implementation provided by us. If you 
plan to use it for non-scientific purposes, don't hesitate to contact us.

Because the programs are licensed free of charge, there is no warranty 
for the program, to the extent permitted by applicable law. Except when 
otherwise stated in writing the copyright holders and/or other parties 
provide the program "as is" without warranty of any kind, either 
expressed or implied, including, but not limited to, the implied 
warranties of merchantability and fitness for a particular purpose. 
The entire risk as to the quality and performance of the program is with 
you. Should the program prove defective, you assume the cost of all 
necessary servicing, repair or correction.

In no event unless required by applicable law or agreed to in writing 
will any copyright holder, or any other party who may modify and/or 
redistribute the program, be liable to you for damages, including any 
general, special, incidental or consequential damages arising out of 
the use or inability to use the program (including but not limited to 
loss of data or data being rendered inaccurate or losses sustained by 
you or third parties or a failure of the program to operate with any 
other programs), even if such holder or other party has been advised 
of the possibility of such damages. 
