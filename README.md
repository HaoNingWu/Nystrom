# Nystrom
This is a collection of MATLAB codes to reproduce all the figures in our paper "Spherical configurations and quadrature methods for integral equations of the second kind," which is going to appear in SIAM Journal on Numerical Analysis.

# Reproducing figures
* Please download the sphere_approx_toolbox_v3.0 (available [here](https://1drv.ms/u/s!AmzdJkQhNBOrhlhZ7TNzdUYOb7X1?e=rfGNGn)) and add it onto path before running the codes.
* Please go to /sphere_approx_toolbox_v3.0/utilities/ and change the bold part of the path '**/Users/haoningwu/Documents/MATLAB/BypaHyper**/sphere_approx_toolbox_v3.0/data/xx' in **loadMD.m**, **loadME.m**, and **loadStd.m** to your own path storing the sphere_approx_toolbox_v3.0. Otherwise, MATLAB would report error:
  >The file '/Users/haoningwu/Documents/MATLAB/BypaHyper/sphere_approx_toolbox_v3.0/data/xx/xxxx' could not be opened because: No such file or
directory

* Please install [Chebfun](https://www.chebfun.org/). We need `legpoly` for the evaluation of modified moments.

# Correspondence
* Please feel free to contact Hao-Ning (email available at [https://haoningwu.github.io/](https://haoningwu.github.io/)) for any questions.


