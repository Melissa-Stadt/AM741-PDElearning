AM741_data.m contains the RMSE and PDEFIND data in the correct form, and calls the relevent plotting functions.

plot_PDEFIND_bar.m plots the learned PDE scores in a bar chart

plot_RMSE.m plots the RMSE values for u, ut, ux, and uxx as a function of noise level for each denoising method

plot_normRMSE_heatmap.m plots a heatmap of RMSE values normalized for each noise level. Depending on the input, this function can compare different denoising techniques, ANN results with a different number of neurons (extension 1), or ANN results with different amounts of data (extension 3).
