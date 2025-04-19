clear; clc; close all
global x_opt
aer_tick= [50,100,200,500,1000,2000,5000,10000,20000,50000];
len_aer = length(aer_tick);
mnt_interval = 4;
mnt_tick = 50:mnt_interval:135;
nmnt = length(mnt_tick);
nregion=3;
expname = 'tke100';
InterFileName=['rbar_epsilon_sigma_',expname];
InputFileName = ['vars_vs_na_4p2_',expname,'.mat'];
load(InputFileName)
[rbar,epsilon,sigma]=deal(zeros(nregion,nmnt*len_aer));
for i=1:nregion
    temp=ydat2(:,:,i);
    rbar(i,:)=temp(:);
    temp=ydat4(:,:,i);
    epsilon(i,:)=temp(:);
    temp=ydat3(:,:,i);
    sigma(i,:)=temp(:);
end
save(InterFileName,'rbar','epsilon','sigma')
%% solve Equation 13
% Initial guess
x0 = [.5; .5; .5; .5; 2]; 
% Optimization options
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
% Run optimization
[x_opt, fval] = fmincon(@objective_function, x0, [], [], [], [], [0; 0; 0; 0; 0], [2; 2; 2; 1; 25], []);
x_opt
%% fig.3
C = figure('position',[597.8,228.2,380,416.8]);
for i_r=1:nregion
    temp=ydat2(:,:,i_r);
    rbar_m(i_r,:)=nanmean(temp);
    rbar_e(i_r,:)=nanstd(temp);
    temp=ydat4(:,:,i_r);
    eps_m(i_r,:)=nanmean(temp);
    eps_e(i_r,:)=nanstd(temp);
end

load(InterFileName)
rbar_fst=min(rbar,[],2);
rbar_lst=max(rbar,[],2);
for i=1:nregion
    % to draw fitting curve, abscissa is needed
    rbar_predict_pic(i,:)=rbar_fst(i):(rbar_lst(i)-rbar_fst(i))/100:rbar_lst(i);
end
eps_predict_pic=epsilon_expression(x_opt,rbar_predict_pic);

func_pic_errorbar_and_plot(rbar_m,eps_m,rbar_e,eps_e,...
    rbar_predict_pic,eps_predict_pic,'  ','linear')
set(gca,'position',[0.168,0.16,0.73,0.631397849462366])
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
ylim([0.12,0.4])
xlim([x_opt(5),25])

% OutFigName = 'epsilon_predict';
% print('-dpng',C,OutFigName,'-r450')
%% fig.4
% close all;
B=figure('position',[488,83,600,450]);
X = repmat(x_opt(5):0.1:20,31,1);
[~,b] = size(X);
sig_s_over_mean_s = repmat([0.2:0.02:0.8]',1,b);
E1 = x_opt(5)^4./X.^4*x_opt(4)^2;
E2 = sig_s_over_mean_s.^2./4.*(1-x_opt(5)^2./X.^2).^2;
% xx = x_opt(5)^2 ./ X.^2;
% yy = x_opt(4)^2;
% zz = sig_s_over_mean_s.^2./4;
% E1 = (xx.^2 + 2 .* zz .* xx .* (1-xx)) .* yy ./ (1 + 2 .* zz .* (1 - xx));
% E2 = zz .* (1-xx) .^ 2 ./ (1 + 2 .* zz .* (1 - xx));
dispersion=sqrt(E1+E2);
func_pic_contour(X,sig_s_over_mean_s,dispersion,1,'  ')
func_pic_contour(X,sig_s_over_mean_s,E1,2,'E_1')
func_pic_contour(X,sig_s_over_mean_s,E2,3,'E_2')
func_pic_contour(X,sig_s_over_mean_s,dispersion.*X,4,'\sigma (\mu m)')
OutFigName = 'theoretical_joint_distribution';
print('-dpng',B,OutFigName,'-r450')

function f = objective_function(x)
expname = 'tke100';
load(['rbar_epsilon_sigma_',expname],'rbar','epsilon')
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

function func_pic_errorbar_and_plot(x_act,y_act,x_err,y_err,x_pre,y_pre,ylb,yscl)
para_xylbl = {'FontSize',16,'interpreter','latex'};
para_axis = {'linewidth',1,'FontName','Times New Roman','FontSize',11};
clr_list = [223,122,94; 227,216,183; 130,178,154; 60,64,91]/255;
clr_list2=[0.900,0.125,0.098; 0.929,0.694,0.125; 0.206,0.774,0.158];
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
xlabel('$\overline{r} (\mu$m)',para_xylbl{:})
ylabel(ylb,para_xylbl{:})%     ͼ  
lgwd = {'0.0$<$AF$\le$0.4, Fit';'0.0$<$AF$\le$0.4, LES';'0.4$<$AF$\le$0.7, Fit';'0.4$<$AF$\le$0.7, LES';...
    '0.7$<$AF$\le$1.0, Fit';'0.7$<$AF$\le$1.0, LES'};
legend(lgwd, ...
    'Location', 'NorthOutside', 'Orientation', 'horizontal', 'NumColumns', 2, ...
    'FontSize', 10, 'FontName', 'Times New Roman','box','off',...
    'orientation','horizontal','interpreter','latex','position',[0.089,0.809,0.894,0.136]);
end

function func_pic_contour(X,Y,Z1,i_p,txtwd)
global x_opt
pstn_list=[0.11,0.60,0.263,0.32;
    0.58,0.60,0.263,0.32;
    0.11,0.11,0.263,0.32;
    0.58,0.11,0.263,0.32];
% rainbow=textread('BlAqGrYeOrReVi200.rgb');
% rainbow=rainbow(:,1:3)/255;
xmin = 2.7;
xmax = 20;
ymin = 0.2;
ymax = 0.8;
nbwd = 'abcdefghijklmn';
ax = subplot(2,2,i_p);
para_xylbl = {'FontSize',13,'FontWeight','bold'};
para_axis = {'linewidth',1,'FontName','Times New Roman','xtick',0:5:100,...
    'FontSize',9};

Z1(X<xmin | X>xmax | Y<ymin | Y>ymax) = nan;
pcolor(X,Y,Z1); shading flat 
colormap('jet')
if i_p==1
caxis([0.09,0.47])
end
cb = colorbar;
cb.FontSize = 7;
v=axis;
xlmt = v(1:2);
dis_sm = [0.717,0.438,0.324];
hold on
plot(xlmt,[dis_sm(1),dis_sm(1)],'--w','LineWidth',1.5)
hold on
plot(xlmt,[dis_sm(2),dis_sm(2)],'--w','LineWidth',1.5)
hold on
plot(xlmt,[dis_sm(3),dis_sm(3)],'--w','LineWidth',1.5)
axis(v)
hold on
min_y = Y(:,1);
min_x = x_opt(5) * sqrt((2*x_opt(4)./min_y).^2+1);
% if i_p==1 || i_p == 4
    plot(min_x,min_y,'w','LineWidth',2.)
% end
ytk = sort([ax.YTick,dis_sm]);
set(gca,para_axis{:},'ytick',ytk,...
    'position',pstn_list(i_p,:),...
    'YTickLabel',num2str(ytk','%.2f'))
text(0.0,1.10,['(',nbwd(i_p),')  ',txtwd],'units','normalized',...
        'FontSize',14,'FontWeight','bold','FontName','Times New Roman')
box on
xlabel('$\overline{r} (\mu m)$','interpreter','latex',para_xylbl{:})
ylabel('$\sigma_{S_{m}}/\overline{S}_{m}$','interpreter','latex',para_xylbl{:})
xlim([2.7,20])
end
