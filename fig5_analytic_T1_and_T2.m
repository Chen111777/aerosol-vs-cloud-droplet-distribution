clear
close all
nbins_r = 50;
nbins_s = 100;
G=10^6;
k_tick = 10.^(log10(0.0005):0.015:log10(3));
t_tick = 10.^(log10(60)-0.02:0.02:4);
s_bar100 = 0.1; % s_bar100 = s_bar (%). If s_bar100 is doubled and t_tick is halved, the result is unchanged.
len_t_tick=length(t_tick);
len_k_tick=length(k_tick);
[record_rbar,record_T1,record_T2] = deal(zeros(len_k_tick,len_t_tick));
%-------initialize droplet spectrum---------------%
n_c=100;
lmd_r = 0.8; % 1.47982500000000;
mu_r = 2.5; % 5.17399999999998;
step_r = 6/(nbins_r-1);
r0 = 2.^(0:step_r:6); % radius before growth, 1~64, unit: um
col = 2.^(0+step_r/2:step_r:6+step_r/2)-...
    2.^(0-step_r/2:step_r:6-step_r/2); % bin width for r
nr0_over_col = r0.^mu_r.*exp(-lmd_r*r0); % gamma distribution
n_r0 = n_c*nr0_over_col.*col/sum(nr0_over_col.*col);
r0_bar = sum(n_r0.*r0)/n_c
sigma0 = sqrt(sum(n_r0.*(r0-r0_bar).^2)/n_c);
epsilon0 = sigma0/r0_bar
i_k = 0;
%-------initialize supersaturation distribution function-----------%
for k=k_tick
    i_k = i_k+1;
    s_bar = s_bar100/100;
    step_s = s_bar*2/(nbins_s-1);
    s = 0:step_s:s_bar*2;
    s_std = s_bar*k;
    f_s = exp(-(s-s_bar).^2./2./(s_std)^2);
    f_s = f_s/sum(f_s); % Eq. 5, truncated normal distribution
    record_epss(i_k) = sqrt(sum(f_s.*(s-s_bar).^2))/s_bar;
%--------------------solve-----------------------%
    r0_mat=repmat(r0*10^-4,nbins_s,1); % unit: cm
    s_mat = repmat(s',1,nbins_r);
    i_t = 0;
    for t=t_tick
        i_t = i_t+1;
        r=sqrt(r0_mat.^2+2.*s_mat.*t./G);
        r_bar_s = sum(n_r0.*r,2)/n_c; % Eq. 10
        sigma2_s = sum(n_r0.*(r-r_bar_s).^2,2)/n_c; % Eq. 11
        r_bar = sum(f_s'.*r_bar_s); % Eq.12
        record_rbar(i_k,i_t) = r_bar;
        record_T1(i_k,i_t) = sum(f_s'.*sigma2_s)./r_bar^2;
        record_T2(i_k,i_t) = sum(f_s'.*(r_bar_s-r_bar).^2)./r_bar^2;
    end
end
%%
%--------epsilon vs. r_bar and epsilon_s contour--------%
figure('position',[488,49.4,600,450]);
X = record_rbar*10^4;
Y = repmat(record_epss',1,len_t_tick);
eps=sqrt(record_T1+record_T2); % Eq. 14
pic_pcolor_epsilon(X,Y,record_T1,1)
pic_pcolor_epsilon(X,Y,record_T2,3)

simp_T1 = r0_bar^4./X.^4*epsilon0^2; % Eq. 21
simp_T2 = Y.^2./4.*(1-r0_bar^2./X.^2).^2; % Eq. 22
simp_eps=sqrt(simp_T1+simp_T2); % Eq. 23
pic_pcolor_epsilon(X,Y,simp_T1,2)
pic_pcolor_epsilon(X,Y,simp_T2,4)
print('-dpng',gcf,'fig5','-r450')

%------------epsilon_s vs. k plot-----------------%
% figure('position',[507.4,301,352.0000000000001,273.6])
% plot(k_tick,record_epss,'k','LineWidth',1.5)
% xlabel('k','FontSize',16,'FontWeight','bold')
% ylabel('\epsilon_S','FontSize',16,'FontWeight','bold')
% set(gca,'linewidth',1,'FontName','Times New Roman')
% grid('on')

function []=pic_pcolor_epsilon(X,Y,Z1,i_p)
pst = [0.12,0.580,0.371,0.352;
    0.57,0.580,0.371,0.352;
    0.12,0.123,0.371,0.352;
    0.57,0.123,0.371,0.352];
titlewd = {'(a1) T_1';'(a2) Simplified T_1';...
    '(b1) T_2';'(b2) Simplified T_2'} ;
ax = subplot('position',pst(i_p,:));
para_xylbl = {'FontSize',16,'FontWeight','bold'};
para_axis = {'linewidth',1,'FontName','Times New Roman'};
if i_p==1 || i_p==2
    contour(X,Y,Z1,[0.0001,0.0002,0.0005,0.001,0.002,0.005,0.01,0.02,0.05,0.1],'showtext','on','LineWidth',1.5); 
elseif i_p==3 || i_p==4
    contour(X,Y,Z1,[0.002,0.01,0.02,0.03:0.02:0.11],'showtext','on','LineWidth',1.5); 
end

% show minimal points
% if i_p<5
%     min_idx = Z1==min(Z1,[],2);
%     hold on; plot(X(min_idx),Y(min_idx),'k','LineWidth',2)
% end

set(gca,para_axis{:},...
    'FontSize',12)
if i_p>2
xlabel('$\mathbf{\overline{r} (\mu m)}$','interpreter','latex',para_xylbl{:})
end
if i_p==1 || i_p==3
ylabel('\epsilon_S',para_xylbl{:})
end
v=axis;
xlim([v(1),40])
caxis([0,0.09])
title(titlewd(i_p),'FontSize',14,'FontWeight','bold')
end