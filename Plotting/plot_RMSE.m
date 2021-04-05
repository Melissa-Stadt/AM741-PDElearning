%%% Inputs 
    % u, ut, ux, uxx are 6x4 (noise x denoising method) matrixes
        % column 1 = FD, 2 = LCVSP, 3= LNCVSP, 4 = ANN
    % title_name is a string containing the title of your plot
        % example: 'advection_diffusion', 'fishers', etc.
%%% Output: RMSE figure

function plot_RMSE(u,ut,ux,uxx,title_name)

n = [0,0.01,0.05,0.1,0.25,0.5];

figure
sgtitle(title_name,'interpreter','latex','fontsize',16) 
subplot(2,2,1)
plot(n,u(:,1),'-o','color',[0.9290 0.6940 0.1250],'linewidth',1,'markerfacecolor',[0.9290 0.6940 0.1250]);
hold on
plot(n,u(:,2),'-o','color',[0.4660 0.6740 0.1880],'linewidth',1,'markerfacecolor',[0.4660 0.6740 0.1880]);
plot(n,u(:,3),'-o','color',[0.3010 0.7450 0.9330],'linewidth',1,'markerfacecolor',[0.3010 0.7450 0.9330]);
plot(n,u(:,4),'-o','color',[0.6350 0.0780 0.1840],'linewidth',1,'markerfacecolor',[0.6350 0.0780 0.1840]);
hold off
legend('FD','LCVSP','LNCVSP','ANN','interpreter','latex','fontsize',12);
title('$\mathbf{u}$','interpreter','latex','fontsize',14);
xlabel('Noise level','interpreter','latex','fontsize',14);
ylabel('RMSE','interpreter','latex','fontsize',14)
subplot(2,2,2)
plot(n,ut(:,1),'-o','color',[0.9290 0.6940 0.1250],'linewidth',1,'markerfacecolor',[0.9290 0.6940 0.1250]);
hold on
plot(n,ut(:,2),'-o','color',[0.4660 0.6740 0.1880],'linewidth',1,'markerfacecolor',[0.4660 0.6740 0.1880]);
plot(n,ut(:,3),'-o','color',[0.3010 0.7450 0.9330],'linewidth',1,'markerfacecolor',[0.3010 0.7450 0.9330]);
plot(n,ut(:,4),'-o','color',[0.6350 0.0780 0.1840],'linewidth',1,'markerfacecolor',[0.6350 0.0780 0.1840]);
hold off
title('$\mathbf{u_t}$','interpreter','latex','fontsize',14);
xlabel('Noise level','interpreter','latex','fontsize',14);
ylabel('log(RMSE)','interpreter','latex','fontsize',14)
set(gca,'yscale','log');
subplot(2,2,3)
plot(n,ux(:,1),'-o','color',[0.9290 0.6940 0.1250],'linewidth',1,'markerfacecolor',[0.9290 0.6940 0.1250]);
hold on
plot(n,ux(:,2),'-o','color',[0.4660 0.6740 0.1880],'linewidth',1,'markerfacecolor',[0.4660 0.6740 0.1880]);
plot(n,ux(:,3),'-o','color',[0.3010 0.7450 0.9330],'linewidth',1,'markerfacecolor',[0.3010 0.7450 0.9330]);
plot(n,ux(:,4),'-o','color',[0.6350 0.0780 0.1840],'linewidth',1,'markerfacecolor',[0.6350 0.0780 0.1840]);
hold off
title('$\mathbf{u_x}$','interpreter','latex','fontsize',14);
xlabel('Noise level','interpreter','latex','fontsize',14);
ylabel('log(RMSE)','interpreter','latex','fontsize',14)
set(gca,'yscale','log');
subplot(2,2,4)
plot(n,uxx(:,1),'-o','color',[0.9290 0.6940 0.1250],'linewidth',1,'markerfacecolor',[0.9290 0.6940 0.1250]);
hold on
plot(n,uxx(:,2),'-o','color',[0.4660 0.6740 0.1880],'linewidth',1,'markerfacecolor',[0.4660 0.6740 0.1880]);
plot(n,uxx(:,3),'-o','color',[0.3010 0.7450 0.9330],'linewidth',1,'markerfacecolor',[0.3010 0.7450 0.9330]);
plot(n,uxx(:,4),'-o','color',[0.6350 0.0780 0.1840],'linewidth',1,'markerfacecolor',[0.6350 0.0780 0.1840]);
hold off
title('$\mathbf{u_{xx}}$','interpreter','latex','fontsize',14);
xlabel('Noise level','interpreter','latex','fontsize',14);
ylabel('log(RMSE)','interpreter','latex','fontsize',14)
set(gca,'yscale','log');

end