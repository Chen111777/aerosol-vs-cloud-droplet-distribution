clear; clc
close all;
G = 10^6;
t = 60 * 20;

%-------initialize droplet spectrum---------------%
nbins_r = 200;
n_c = 100;
B=figure('position',[315.4,502.6,650,259.4]);
%%
mu_ln_r0 = log(1.87); % 对数正态分布的对数均值，调整以满足需求
sigma_ln_r0 = 0.57;  % 对数正态分布的对数标准差，调整以满足需求

step_r = 5 / (nbins_r - 1);
r0 = 2.^(0:step_r:5); % radius before growth, 1~64, unit: um
col = step_r; % bin width for r

% 对数正态分布的 PDF
ln_r0 = log(r0); % 计算对数半径
nr0_over_col = exp(-((ln_r0 - mu_ln_r0).^2) / (2 * sigma_ln_r0^2)); % log-normal distribution

% 归一化
n_r0 = n_c * nr0_over_col .* col / sum(nr0_over_col .* col);

% 计算统计量
r0_bar = sum(n_r0 .* r0) / n_c % 平均半径
sigma0 = sqrt(sum(n_r0 .* (r0 - r0_bar).^2) / n_c); % 标准差
epsilon0 = sigma0 / r0_bar % 相对离散度

%-------initialize PDF of S ---------------%
nbins_s = 200;
sigma_ln_s_tick = 0.1:0.025:0.75; % lower limit determine the upper limit of epsilon_s
mu_ln_s_tick = -8:0.03:0; % lower limit determine the upper limit of r_bar
len_mu = length(sigma_ln_s_tick);
len_lmd = length(mu_ln_s_tick);
step_s = 5/(nbins_s-1);
s = 10.^(-4:step_s:1);
col_s = step_s; % bin width for r
[record_eps_s,record_rbar,record_E1,record_E2] = deal(zeros(len_lmd,len_mu));
imu = 0;
for sigma_ln_s = sigma_ln_s_tick
    imu = imu+1;
    ilmd = 0;
    for mu_ln_s = mu_ln_s_tick
        ilmd = ilmd+1;
        ln_s = log(s); % 计算对数s
        fs_over_col = exp(-((ln_s - mu_ln_s).^2) / (2 * sigma_ln_s^2)); % gamma distribution
        f_s = fs_over_col.*col_s/sum(fs_over_col.*col_s);
        s_bar = sum(f_s.*s);
        record_eps_s(ilmd,imu) = sqrt(sum(f_s.*(s-s_bar).^2))/s_bar;
    %--------------------solve-----------------------%
        r0_mat=repmat(r0*10^-4,nbins_s,1); % unit: cm
        s_mat = repmat(s'/100,1,nbins_r);
        r=sqrt(r0_mat.^2+2.*s_mat.*t./G);
        r_bar_s = sum(n_r0.*r,2)/n_c; % Eq. 10
        sigma2_s = sum(n_r0.*(r-r_bar_s).^2,2)/n_c; % Eq. 11
        r_bar = sum(f_s'.*r_bar_s); % Eq.12
        record_rbar(ilmd,imu) = r_bar;
        record_E1(ilmd,imu) = sum(f_s'.*sigma2_s)./r_bar^2;
        record_E2(ilmd,imu) = sum(f_s'.*(r_bar_s-r_bar).^2)./r_bar^2;
    end
end

X = record_rbar*10^4;
Y = record_eps_s;
E1 = record_E1;
E2 = record_E2;
Z=sqrt(E1+E2);
func_pic_contour(X,Y,Z,1,'Lognormal')
%
lmd_r = 0.90; % larger lmd_r, smaller r0_bar, smaller epslon
mu_r = 0.59; % larger mu_r, larger r0_bar, smaller epsilon
step_r = 5/(nbins_r-1);
r0 = 2.^(0:step_r:5); % radius before growth, 1~64, unit: um
col = 2.^(0+step_r/2:step_r:5+step_r/2)-...
    2.^(0-step_r/2:step_r:5-step_r/2); % bin width for r
nr0_over_col = r0.^mu_r.*exp(-lmd_r*r0); % gamma distribution
n_r0 = n_c*nr0_over_col.*col/sum(nr0_over_col.*col);
r0_bar = sum(n_r0.*r0)/n_c
sigma0 = sqrt(sum(n_r0.*(r0-r0_bar).^2)/n_c);
epsilon0 = sigma0/r0_bar

%-------initialize PDF of S ---------------%
nbins_s = 200;
mu_tick = 2.^(-3:0.2:8); % lower limit determine the upper limit of epsilon_s
lmd_tick = 2.^(3:0.04:15.5); % lower limit determine the upper limit of r_bar
len_mu = length(mu_tick);
len_lmd = length(lmd_tick);
step_s = 5/(nbins_s-1);
s = 10.^(-4:step_s:1);
col_s = 10.^(-4+step_s/2:step_s:1+step_s/2)-...
    10.^(-4-step_s/2:step_s:1-step_s/2); % bin width for r
[record_eps_s,record_rbar,record_E1,record_E2] = deal(zeros(len_lmd,len_mu));
imu = 0;
for mu = mu_tick
    imu = imu+1;
    ilmd = 0;
    for lmd = lmd_tick
        ilmd = ilmd+1;
        fs_over_col = s.^mu.*exp(-lmd*s); % gamma distribution
        f_s = fs_over_col.*col_s/sum(fs_over_col.*col_s);
        s_bar = sum(f_s.*s);
        record_eps_s(ilmd,imu) = sqrt(sum(f_s.*(s-s_bar).^2))/s_bar;
    %--------------------solve-----------------------%
        r0_mat=repmat(r0*10^-4,nbins_s,1); % unit: cm
        s_mat = repmat(s'/100,1,nbins_r);
        r=sqrt(r0_mat.^2+2.*s_mat.*t./G);
        r_bar_s = sum(n_r0.*r,2)/n_c; % Eq. 10
        sigma2_s = sum(n_r0.*(r-r_bar_s).^2,2)/n_c; % Eq. 11
        r_bar = sum(f_s'.*r_bar_s); % Eq.12
        record_rbar(ilmd,imu) = r_bar;
        record_E1(ilmd,imu) = sum(f_s'.*sigma2_s)./r_bar^2;
        record_E2(ilmd,imu) = sum(f_s'.*(r_bar_s-r_bar).^2)./r_bar^2;
    end
end

X = record_rbar*10^4;
Y = record_eps_s;
E1 = record_E1;
E2 = record_E2;
Z=sqrt(E1+E2);
func_pic_contour(X,Y,Z,2,'Gamma')


OutFigName = 'support_2distributions';
% print('-dpng',gcf,OutFigName,'-r450')
%%
function func_pic_contour(X,Y,Z,i_p,txtwd)
pstn_list=[0.14,0.1677,0.2926,0.6804;
    0.5536,0.1677,0.2926,0.6804];
xmin = 2.7;
xmax = 20;
ymin = 0.2;
ymax = 0.8;
nbwd = 'abcdefghijklmn';
ax = subplot(1,2,i_p);
para_xylbl = {'FontSize',13,'FontWeight','bold'};
para_axis = {'linewidth',1,'FontName','Times New Roman','xtick',0:5:100,...
    'FontSize',10};

Z1 = Z;
Z1(X<xmin | X>xmax | Y<ymin | Y>ymax+0.01) = nan;
pcolor(X,Y,Z1); shading flat 
colormap('jet')
caxis([0.09,0.47])
if i_p==2
colorbar
end

hold on

y_range = ymin:0.05:ymax;
x_range = xmin:0.5:xmax;
[xq, yq] = meshgrid(x_range, y_range); % 生成新的横纵坐标网格
values_interp = griddata(X(~isnan(Z)), Y(~isnan(Z)), Z(~isnan(Z)), xq, yq, 'cubic'); % 插值
[~, idx_min] = min(values_interp, [], 2); % 找到每一行的最小值索引
x_min = x_range(idx_min); 
plot(x_min,y_range,'w','LineWidth',1.5)

ax = gca;
set(gca,para_axis{:},...
    'position',pstn_list(i_p,:))
text(0.0,1.10,['(',nbwd(i_p),')  ',txtwd],'units','normalized',...
        'FontSize',14,'FontWeight','bold','FontName','Times New Roman')
box on
xlabel('$\overline{r} (\mu m)$','interpreter','latex',para_xylbl{:})
ylabel('$\sigma_{S_{m}}/\overline{S}_{m}$','interpreter','latex',para_xylbl{:})
xlim([xmin,xmax])
ylim([ymin,ymax])
end