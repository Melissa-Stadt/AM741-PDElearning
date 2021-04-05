%%% Inputs 
    % u, ut, ux, uxx are 6x4 (noise x X) matrixes
        % X = denoising method, number of neurons, % data reduction
    % x_var = string specifying the x variable
        % Options:
            % 'denoise'
            % 'neurons'
            % 'reduction'
    % title_name is a string containing the title of your plot
%%% Output: heatmap of RMSE values normalized for each noise level figure

function plot_normRMSE_heatmap(u,ut,ux,uxx,x_var,title_name)

if strcmp(x_var,'denoise')
    x_axis = 'Denoising Method';
    x_labels = {'FD','LCVSP','LNCVSP','ANN'};
elseif strcmp(x_var,'neurons')
    x_axis = 'Number of Neurons';
    x_labels = {'500','1000','2000','3000'};
elseif strcmp(x_var,'reduction')
    x_axis = '% Reduction in Data';
    x_labels = {'0','10','50','80'};
end


% Normalize RMSE for each noise level between 0 and 1
norm_u = (u - min(u,[],2))./(max(u,[],2)-min(u,[],2));
norm_ut = (ut - min(ut,[],2))./(max(ut,[],2)-min(ut,[],2));
norm_ux = (ux - min(ux,[],2))./(max(ux,[],2)-min(ux,[],2));
norm_uxx = (uxx - min(uxx,[],2))./(max(uxx,[],2)-min(uxx,[],2));

% Noise level
y_labels = categorical({'0.00','0.01','0.05','0.10','0.25','0.50'});
y_axis = 'Noise Level';

figure
sgtitle(title_name,'fontsize',14,'fontweight','bold');
subplot(2,2,1)
heatmap(x_labels,y_labels,norm_u);
title('u');
xlabel(x_axis);
ylabel(y_axis);
subplot(2,2,2)
heatmap(x_labels,y_labels,norm_ut);
title('u_t');
xlabel(x_axis);
ylabel(y_axis);
subplot(2,2,3)
heatmap(x_labels,y_labels,norm_ux);
title('u_x');
xlabel(x_axis);
ylabel(y_axis);
subplot(2,2,4)
heatmap(x_labels,y_labels,norm_uxx);
title('u_{xx}');
xlabel(x_axis);
ylabel(y_axis);

end