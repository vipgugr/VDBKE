# VDBKE

# Variational Dirichlet Blur Kernel Estimation

A Matlab program that implements the variational Dirichlet blur kernel estimation method in X. Zhou, J. Mateos, F. Zhou, R. Molina, and A.K. Katsaggelos, “Variational Dirichlet Blur Kernel Estimation”, IEEE Transactions on Image Processing, vol. 24, no. 12, 5127-5139, December 2015. doi: [10.1109/TIP.2015.2478407](http://doi.org/10.1109/TIP.2015.2478407)
[[![DownloadPDF](http://decsai.ugr.es/vip/images/pdficon.gif)(2.86MB)](http://decsai.ugr.es/vip/files/journals/Variation%20Dirichlet%20Blur%20Kernel%20Estimation.pdf)]

## Abstract

Blind image deconvolution involves two key objectives, latent image and blur estimation. For latent image estimation, we propose a fast deconvolution algorithm, which uses an image prior of nondimensional Gaussianity measure to enforce sparsity and an undetermined boundary condition methodology to reduce boundary artifacts. For blur estimation, a linear inverse problem with normalization and nonnegative constraints must be solved. However, the normalization constraint is ignored in many blind image deblurring methods, mainly because it makes the problem less tractable. In this paper, we show that the normalization constraint can be very naturally incorporated into the estimation process by using a Dirichlet distribution to approximate the posterior distribution of the blur. Making use of variational Dirichlet approximation, we provide a blur posterior approximation that takes into account the uncertainty of the estimate and removes noise in the estimated kernel. Experiments with synthetic and real data demonstrate that the proposed method is very competitive to state-of-the-art blind image restoration methods.

## Examples

### Motion blur
		
| ![Input](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/Blurry3_8.png)  | ![Xu](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/xu0_3_8.png)![XuPSF](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/kernel_xu0_3_8.png) | ![Ours](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/ngm_3_8.png)![OurPSF](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/kernel_ngm_3_8.png)  |
|-------|---------------|------|
| Input |	Xu et al. [1] |	Ours |


 
### Atmospheric blur
| ![Input](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/y7.png) | ![Babacan](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/kdeb_derin_7.png) | ![Xu](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/kdeb_xu_7.png) | ![Ors](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/kdeb_zhou_7.png) |
|-------|--------------------|----------------|------|
| Input |	Babacan et al. [2] |	Xu et al. [1] |	Ours |

 
### Numerical results of Algorithm 1:

[Zip file containing the numerical results of Algorithm 1](http://decsai.ugr.es/vip/resources/VDBlurKernelEstimation/plots.zip)

## References

[1] L. Xu, et al., “Unnatural l0 sparse representation for natural image deblurring,” in CVPR, 2013.

[2] S. Babacan, et al., “Bayesian blind deconvolution with general sparse image priors,” in ECCV, 2012.

## Disclaimer

The programs are granted free of charge for research and education purposes only. Scientific results produced using the software provided shall acknowledge the use of the implementation provided by us. If you plan to use it for non-scientific purposes, don't hesitate to contact us.

Because the programs are licensed free of charge, there is no warranty for the program, to the extent permitted by applicable law. except when otherwise stated in writing the copyright holders and/or other parties provide the program "as is" without warranty of any kind, either expressed or implied, including, but not limited to, the implied warranties of merchantability and fitness for a particular purpose. The entire risk as to the quality and performance of the program is with you. Should the program prove defective, you assume the cost of all necessary servicing, repair or correction.

In no event unless required by applicable law or agreed to in writing will any copyright holder, or any other party who may modify and/or redistribute the program, be liable to you for damages, including any general, special, incidental or consequential damages arising out of the use or inability to use the program (including but not limited to loss of data or data being rendered inaccurate or losses sustained by you or third parties or a failure of the program to operate with any other programs), even if such holder or other party has been advised of the possibility of such damages. 
