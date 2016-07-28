# Minding the Gaps for Block Frank-Wolfe Optimization of Structured SVMs

This is a Matlab implementation of the Structured SVM (SSVM) solvers proposed in the ICML-2016 paper [Minding the Gaps for Block Frank-Wolfe Optimization of Structured SVMs](http://jmlr.org/proceedings/papers/v48/osokin16-supp.pdf).

This code is based on the [BCFWstruct](https://github.com/ppletscher/BCFWstruct) library and is organized in a similar way:

* `solvers` contains the optimization methods
* `applications` contains the application-dependent code, such as MAP decoding or the feature map computation. The source code includes the following application:
	- sequence prediction for Optical Character Recognition (OCR), [OCR dataset](http://ai.stanford.edu/~btaskar/ocr/)
	- sequence prediction for text chunking, [CONLL dataset](http://www.cnts.ua.ac.be/conll2000/chunking/)
	- binary image segmentation, [HorseSeg dataset](https://pub.ist.ac.at/~akolesnikov/HDSeg/)
	- human pose estimation, [LSP dataset](http://www.comp.leeds.ac.uk/mat4saj/lsp.html)
* `experiments` contains scripts to reproduce experiments of the [ICML 2016 paper](http://jmlr.org/proceedings/papers/v48/osokin16-supp.pdf).
* `data` is initally empty and is used to store the data files required for the demos.

The code is released under Apache v2 License allowing to use the code in any way you want.
All the datasets are available on the [project webpage](http://www.di.ens.fr/sierra/research/gapBCFW/).

### Citation

If you are using this software please cite the following paper in any resulting publication:

>@InProceedings{osokin2016gapBCFW, <br>
    author      = {Osokin, Anton and Alayrac, Jean-Baptiste and Lukasewitz, Isabella and Dokania, Puneet K. and Lacoste-Julien, Simon},<br>
    title       = {Minding the Gaps for Block {F}rank-{W}olfe Optimization of Structured {SVM}s},<br>
    booktitle   = {Proceedings of The 33rd International Conference on Machine Learning (ICML)},<br>
    year        = {2016} }

### Authors

* [Anton Osokin](http://www.di.ens.fr/~osokin/)*
* [Jean-Baptiste Alayrac](http://www.di.ens.fr/~alayrac/)*
* [Isabella Lukasewitz](https://www.linkedin.com/in/isabella-lukasewitz-19159359)
* [Puneet K. Dokania](http://puneetkdokania.github.io/)
* [Simon Lacoste-Julien](http://www.di.ens.fr/~slacoste/)

\* Both authors contributed equally.

### Installation

1. You need a working installation of MATLAB (with [C++ compiler configured](http://www.mathworks.com/help/matlab/matlab_external/what-you-need-to-build-mex-files.html)).

2. Clone the git repo (`git clone https://github.com/aosokin/gapBCFW.git`) or download the [zip archive](https://github.com/aosokin/gapBCFW/archive/master.zip).

3. Obtain the data files required to run the demos. You can use our scripts to download the data: `application/OCR/download_ocr.m`, `application/text_chunking/download_conll.m`, `application/binary_segmentation/download_horseSeg.m`, `applications/pose_estimation/download_LSP.m`

4. Compile the binaries using `compile_BCFW.m`

5. Run `application/OCR/demo_ocr.m`, `application/text_chunking/demo_conll.m`, `application/binary_segmentation/demo_horseSeg.m`, `applications/pose_estimation/demo_LSP.m`

The code was tested on Ubuntu 12.04, gcc 4.6.3, Matlab-R2014b, but should run on other systems as well.

### Usage

If you want to use our solvers for your own structured prediction problem you will need to implement several functions:

* The feature map.
* The maximization oracle.
* The loss function.
* (required for some regimes) The label hashing function.

You can find example implementations in the `applications` folder. For an overview of the exact usage and supported options, please check the Matlab documentation of the solvers: `solvers/solver_BCFW_hybrid.m` and `solvers/solver_multiLambda_BCFW_hybrid.m`.

Note that the interface of our code is very similar to the interfaces of [BCFWstruct](https://github.com/ppletscher/BCFWstruct) and the [Matlab wrapper to SVM^struct](http://www.robots.ox.ac.uk/~vedaldi//svmstruct.html).
Users of these packages can easily use our package as well.

### Reproducing our experiments

Folder `experiments` contains code to reproduce Figures 3, 5, 6 of the [ICML 2016 paper](http://jmlr.org/proceedings/papers/v48/osokin16-supp.pdf).

To just reproduce the plots, go to subfolder `experiments/plots_icml2016`, download our result files by running `download_results.m`, run plotting scripts
* `plots_icml2016_BCWH_ocrLarge.m`
* `plots_icml2016_BCWH_connl.m`
* `plots_icml2016_BCWH_horseSegSmall.m`
* `plots_icml2016_BCWH_horseSegMedium.m`
* `plots_icml2016_BCWH_horseSegLarge.m`
* `plots_icml2016_BCWH_lspSmall.m`.

To rerun the whole experiment, use scripts 
* `experiments/ocr_dataset/ocr_large_BCFW_hybrid_run_experiments.m`
* `experiments/conll_dataset/conll_BCFW_hybrid_run_experiments.m`
* `experiments/horseSeg_dataset/horse_small_BCFW_hybrid_run_experiments.m`
* `experiments/horseSeg_dataset/horse_medium_BCFW_hybrid_run_experiments.m`
* `experiments/horseSeg_dataset/horse_large_BCFW_hybrid_run_experiments.m`
* `experiments/LSP_dataset/lsp_small_BCFW_hybrid_run_experiments.m`

Note, that to run all experiments you will need significant computational resources. We recommend using a cluster to run the experiments in reasonable time.