This folder contains the code that denoise data and approximate their partial derivatives with various methods.

Run **DEMO.py** specifying the desired user inputs to perform the data denoising part of the study. Groundtruth files are in the Data/ folder and noise can get added using appropriate 'generate_' file.

**DEMO.py** - contains the code to train a new surface fitting artificial neural network (ANN) or make predictions using the ANN, finite difference, spline, local bi-spline, local NCV bi-spline, or global NCV bi-spline methods. All ANN parameters are stored in the checkpoints folder. All predictions are automatically stored in the data folder. 

**make_rmse_tables.py** - Run to generate RMSE tables comparing the different denoising methods. This code computes relative mean square errors (RMSEs) between the true function/derivative values and the approximations.

**make_animations.py** - Run to make animations comparing the true solution, predicted solution, and the data. This code saves a number of plots in the animations folder and combines them into a .gif file. The plots are deleted after the .gif file is complete.

**prediction_functions.py** - contains the code to forward evaluate the ANN, finite difference, spline, local bi-spline, local NCV bi-spline, and global NCV bi-spline methods

**custom_functions.py** - contains functions used for training new ANN models. 

**surface_fitter.py** - contains the SurfNN class which handles loading data, training new ANNs, and making predictions. This code can also save figures in the plots folder during training to illustrate ANN convergence.

**spline_custom_functions.py** - contains the code used for implementing the global NCV bi-spline method

**generate_advection_diffusion_data.py, generate_fisher_data.py, generate_fisher_nonlin_data.py**  - code to add proportional error to precomputed noiseless datasets and save the noisy versions in the data folder.

**generate_new_general_data.py** - code to add proportional error to precomputed noiseless datasets given by prefix + '\_groundtruth.mat' and saves the noisy versions in the data file. Requires user input for the name of the groundtruth file

**generate Fisher groundtruth/**- this folder contains the files used to create groundtruth data for general Fisher KPP
