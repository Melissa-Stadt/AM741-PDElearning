%%% Inputs
    % score = 6x4 (noisexdenoising method) matrix storing the learned PDEs
        % column 1 = FD, 2 = LCVSP, 3= LNCVSP, 4 = ANN
    % title_name = string containing the title of your plot
%%% Output: bar plot of PDEFIND scores 
function plot_PDEFIND_bar(score,title_name)

xlabs = categorical({'0.00','0.01','0.05','0.10','0.25','0.50'});

figure
b = bar(xlabs,score);
b(1).FaceColor = [0.9290 0.6940 0.1250];
b(2).FaceColor = [0.4660 0.6740 0.1880];
b(3).FaceColor = [0.3010 0.7450 0.9330];
b(4).FaceColor = [0.6350 0.0780 0.1840];
legend('FD','LCVSP','LNCVSP','ANN');
set(legend,'interpreter','latex','fontsize',12);
set(gca,'fontsize',12);
xlabel('Noise Level','interpreter','latex','fontsize',14);
ylabel('PDE Score','interpreter','latex','fontsize',14);
title(title_name,'interpreter','latex','fontsize',14);
set(gca,'yscale','log');
ylim([0,1.1]);

end