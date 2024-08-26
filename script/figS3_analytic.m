clear
G=10^6;
t = 60*20;
%-------initialize droplet spectrum---------------%
nbins_r = 100;
n_c=100;
lmd_r = 0.90; % larger lmd_r, smaller r0_bar, smaller epslon
mu_r = 0.60; % larger mu_r, larger r0_bar, smaller epsilon
step_r = 5/(nbins_r-1);
r0 = 2.^(0:step_r:5); % radius before growth, 1~64, unit: um
col = 2.^(0+step_r/2:step_r:5+step_r/2)-...
    2.^(0-step_r/2:step_r:5-step_r/2); % bin width for r
nr0_over_col = r0.^mu_r.*exp(-lmd_r*r0); % gamma distribution
n_r0 = n_c*nr0_over_col.*col/sum(nr0_over_col.*col);
r0_bar = sum(n_r0.*r0)/n_c
sigma0 = sqrt(sum(n_r0.*(r0-r0_bar).^2)/n_c);
epsilon0 = sigma0/r0_bar
%%
%-------initialize PDF of S ---------------%
nbins_s = 200;
mu_tick = 2.^(-1:0.3:8); % lower limit determine the upper limit of epsilon_s
lmd_tick = 2.^(2:0.3:15.5); % lower limit determine the upper limit of r_bar
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
%%
close all;
figure('position',[488,83,560,690])
x_opt=[0.731536370882398;0.449332201293307;...
    0.339500258052524;0.546134418674883;...
    2.42153922153402];
X = record_rbar*10^4;
Y = record_eps_s;
E1 = record_E1;
E2 = record_E2;
Z=sqrt(E1+E2);
pic_contour(X,Y,Z,1,'¦Å from Eq. 10')
pic_contour(X,Y,E1,2,'E_1 from Eq. 10')
pic_contour(X,Y,E2,3,'E_2 from Eq. 10')
pic_contour(X,Y,E1./Z.^2,4,'E_1/¦Å^2 from Eq. 10')
pic_contour(X,Y,Z.*X,5,'\sigma from Eq. 9')
OutFigName = 'support_theoretical_joint_distribution';
% print('-dpng',gcf,OutFigName,'-r450')


%%
function pic_contour(X,Y,Z1,i_p,txtwd)
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
ylim([0.2,0.8])
end