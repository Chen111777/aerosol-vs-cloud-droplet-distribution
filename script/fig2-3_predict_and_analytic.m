clear; clc; close all
aer_tick= [50,100,200,500,1000,2000,5000,10000,20000,50000];

mnt_interval = 3;
mnt_tick = 75:mnt_interval:135;
nmnt = length(mnt_tick);
len_aer = length(aer_tick);
nregion=3;
InterFileName='rbar_epsilon_sigma';
InputFileName = 'vars_vs_na_4p2_exp1.mat';
load(InputFileName)
% [rbar,epsilon,sigma]=deal(zeros(nregion,nmnt*len_aer));
% for i=1:nregion
%     temp=ydat2(:,:,i);
%     rbar(i,:)=temp(:);
%     temp=ydat3(:,:,i);
%     epsilon(i,:)=temp(:);
%     temp=ydat4(:,:,i);
%     sigma(i,:)=temp(:);
% end
% save(InterFileName,'rbar','epsilon','sigma')
%% solve Equation 13
% Initial guess
x0 = [.5; .5; .5; .5; 2]; 
% Optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
% Run optimization
[x_opt, fval] = fmincon(@objective_function, x0, [], [], [], [], [0; 0; 0; 0; 0], [2; 2; 2; 1; 25], []);
x_opt
%% fig.2
C = figure('position',[597.8,228.2,850.2,301.8]);
for i_r=1:nregion
    temp=ydat2(:,:,i_r);
    rbar_m(i_r,:)=nanmean(temp);
    rbar_e(i_r,:)=nanstd(temp);
    temp=ydat3(:,:,i_r);
    eps_m(i_r,:)=nanmean(temp);
    eps_e(i_r,:)=nanstd(temp);
    temp=ydat4(:,:,i_r);
    sig_m(i_r,:)=nanmean(temp);
    sig_e(i_r,:)=nanstd(temp);
end

load('rbar_epsilon_sigma')
rbar_fst=min(rbar,[],2)+0.2;
rbar_lst=max(rbar,[],2);
for i=1:nregion
    % to draw fitting curve, abscissa is needed
    rbar_predict_pic(i,:)=rbar_fst(i):(rbar_lst(i)-rbar_fst(i))/100:rbar_lst(i);
end
eps_predict_pic=epsilon_expression(x_opt,rbar_predict_pic);

func_pic_errorbar_and_plot(rbar_m,eps_m,rbar_e,eps_e,...
    rbar_predict_pic,eps_predict_pic,1,'¦Å','linear')
eps_predict=epsilon_expression(x_opt,rbar);
eps_predict=eps_predict(:);
eps_predict(isnan(eps_predict))=[];
eps_actual=epsilon(:);
eps_actual(isnan(eps_actual))=[];
rms=sqrt(sum((eps_predict-eps_actual).^2)/length(eps_actual))
r2=corrcoef(eps_predict,eps_actual);
r2=r2(1,2)
text(16,0.305,{['R^2=',num2str(r2,'%.2f')];['RMSE=',num2str(rms,'%.2f')]},...
    'FontName','Times New Roman','fontsize',13)
ylim([0.13,0.42])
xlim([x_opt(5),25])

func_pic_errorbar_and_plot(rbar_m,sig_m,rbar_e,sig_e,...
    rbar_predict_pic,eps_predict_pic.*rbar_predict_pic,2,'$\sigma (\mu$m)','log')
sig_predict=epsilon_expression(x_opt,rbar).*rbar;
sig_predict=sig_predict(:);
sig_predict(isnan(sig_predict))=[];
sig_actual=sigma(:);
sig_actual(isnan(sig_actual))=[];
rms=sqrt(sum((sig_predict-sig_actual).^2)/length(sig_actual))
r2=corrcoef(sig_predict,sig_actual);
r2=r2(1,2)
text(3.329,3.045,0,{['R^2=',num2str(r2,'%.2f')];['RMSE=',num2str(rms,'%.2f')]},...
    'FontName','Times New Roman','fontsize',13)
xlim([3,25])
legend({'0.0<AF¡Ü0.4, Predicted';'0.0<AF¡Ü0.4, LES data';'0.4<AF¡Ü0.7, Predicted';'0.4<AF¡Ü0.7, LES data';...
    '0.7<AF¡Ü1.0, Predicted';'0.7<AF¡Ü1.0, LES data'},...
    'box','off','position',[0.755,0.328,0.225,0.454],'FontSize',10)
OutFigName = 'epsilon_sigma_predict';
% print('-dpng',C,OutFigName,'-r450')
%% fig.3
close all;
B=figure('position',[488,83,560,690]);
X = repmat(x_opt(5):0.3:25.3,31,1);
[~,b] = size(X);
Y = repmat([0.2:0.02:0.8]',1,b);
E1 = x_opt(5)^4./X.^4*x_opt(4)^2;
E2 = Y.^2./4.*(1-x_opt(5)^2./X.^2).^2;
Z=sqrt(E1+E2);
func_pic_contour(X,Y,Z,1,'¦Å from Eq. 11')
func_pic_contour(X,Y,E1,2,'E_1 from Eq. 11')
func_pic_contour(X,Y,E2,3,'E_2 from Eq. 11')
func_pic_contour(X,Y,E1./Z.^2,4,'E_1/¦Å^2 from Eq. 11')
func_pic_contour(X,Y,Z.*X,5,'\sigma from Eq. 12')
OutFigName = 'theoretical_joint_distribution';
% print('-dpng',B,OutFigName,'-r450')

function f = objective_function(x)
load('rbar_epsilon_sigma','rbar','epsilon')
    f=0;
    for i=1:3
        rbar0orbar1_sq=x(5)^2./rbar(i,:).^2;
        epsilon_pre=sqrt(x(4)^2*rbar0orbar1_sq.^2+...
                        x(i).^2/4.*(1-rbar0orbar1_sq).^2);
        f = f+nansum((epsilon(i,:)-epsilon_pre).^2) ;
    end
end

function epsilon_pre = epsilon_expression(x,rbar)
    epsilon_pre=zeros(size(rbar));
    for i=1:3
        rbar0orbar1_sq=x(5)^2./rbar(i,:).^2;
        epsilon_pre(i,:)=sqrt(x(4)^2*rbar0orbar1_sq.^2+...
                        x(i).^2/4.*(1-rbar0orbar1_sq).^2);
    end
end

function func_pic_errorbar_and_plot(x_act,y_act,x_err,y_err,x_pre,y_pre,i_p,ylb,yscl)
pstn_list=[0.09,0.17,0.295,0.733;
    0.455,0.17,0.295,0.733];
para_xylbl = {'FontSize',14,'interpreter','latex'};
para_axis = {'linewidth',1,'FontName','Times New Roman','FontSize',11};
nbwd = 'abcdefghijklmn';
clr_list = [223,122,94; 227,216,183; 130,178,154; 60,64,91]/255;
clr_list2=[0.900,0.125,0.098; 0.929,0.694,0.125; 0.206,0.774,0.158];
subplot('position',pstn_list(i_p,:))
[nregion,~]=size(x_act);
for i_r=1:nregion
    hold on; 
    plot(x_pre(i_r,:),y_pre(i_r,:),'LineWidth',2,'color',clr_list2(i_r,:))
    hold on; 
    errorbar(x_act(i_r,:),y_act(i_r,:),y_err(i_r,:),y_err(i_r,:),x_err(i_r,:),x_err(i_r,:),...
        'color',clr_list(i_r,:)*3/5,'LineStyle','none','LineWidth',1,'Marker','o',...
        'MarkerFaceColor','w','MarkerEdgeColor','k','MarkerSize',3)
end

box('on')
grid('on')
set(gca,para_axis{:},'yscale',yscl,'xscale',yscl)
if strcmp(yscl,'log')
    set(gca,'ytick',1:20,'xtick',0:5:30)
end
title(['(',nbwd(i_p),') '],'FontSize',15,'FontWeight','bold',...
       'units','normalized','position',[0.106,0.861,0])
xlabel('$\overline{r} (\mu$m)',para_xylbl{:})
ylabel(ylb,para_xylbl{:})
end

function func_pic_contour(X,Y,Z1,i_p,txtwd)
pstn_list=[0.11,0.73,0.34,0.22;
    0.56,0.73,0.34,0.22;
    0.11,0.405,0.34,0.22;
    0.56,0.405,0.34,0.22;
    0.11,0.08,0.34,0.22;
    0.56,0.08,0.34,0.22];
nbwd = 'abcdefghijklmn';
subplot('position',pstn_list(i_p,:));
para_xylbl = {'FontSize',12,'FontWeight','bold'};
para_axis = {'linewidth',1,'FontName','Times New Roman','xtick',0:5:100,...
    'FontSize',10};

if i_p==2 || i_p==4
    hold on
    contour(X,Y,Z1,[0.0001,0.0002,0.0005,0.001,0.002,0.005,...
        0.01,0.02,0.05,0.1,0.2,0.5],'showtext','on','LineWidth',1.5,'labelspacing',100);
elseif i_p==5
    contour(X,Y,Z1,1:1:10,'showtext','on','LineWidth',1.5,'labelspacing',100);
else
    contour(X,Y,Z1,'showtext','on','LineWidth',1.5); 
end
set(gca,para_axis{:})
text(0.05,1.08,['(',nbwd(i_p),') ',txtwd],'units','normalized',...
        'FontSize',13,'FontWeight','bold','FontName','Times New Roman')
box on
grid on
xlabel('$\overline{r} (\mu m)$','interpreter','latex',para_xylbl{:})
ylabel('$\sigma_{S_{m}}/\overline{S}_{m}$','interpreter','latex',para_xylbl{:})
xlim([2.7,25])
end