clear; clc; close all
paths = {'~/WRFV4.5.1/a50';'~/WRFV4.5.1/a100';'~/WRFV4.5.1/a200';'~/WRFV4.5.1/a500';...
'~/WRFV4.5.1/a1000';'~/WRFV4.5.1/a2000';'~/WRFV4.5.1/a5000';...
'~/WRFV4.5.1/a10000';'~/WRFV4.5.1/a20000';'~/WRFV4.5.1/a50000'};
aer_tick= [50,100,200,500,1000,2000,5000,10000,20000,50000];
expname = 'exp1';

global pstn_list
mnt_interval = 3;
mnt_tick = 75:mnt_interval:135;
InputFileName = ['lwcad_',expname,'_75to135.mat'];
InterFileName = ['dat_DispersionS',expname,'.mat'];
OutFigName = ['Supersaturation_',expname];
len_aer = length(aer_tick);

pstn_list=[0.08,0.21,0.25,0.75;
    0.395,0.21,0.25,0.75;
    0.71,0.21,0.25,0.75];
%% prepare data for figure
load(InputFileName)
[dat_epss,dat_sbar,dat_sigs] = func_cloudy_region_supertsaturation(paths,mnt_tick,lwcad);
save(InterFileName,'dat_epss','dat_sigs','dat_sbar')
%% Figure 2
load(InterFileName)
B = figure('position',[488,258.6,850,270]);
[Y1,Y2,Y3] = deal(zeros(len_aer,3));
for ia = 1:len_aer
    Y1(ia,:) = quantile(dat_epss,[0.25,0.5,0.75]);
    Y2(ia,:) = quantile(dat_sbar,[0.25,0.5,0.75]);
    Y3(ia,:) = quantile(dat_sigs,[0.25,0.5,0.75]);
end
func_fig_shade(aer_tick,Y2,1,[50,50000],'$N_{a}$ (cm$^{-3}$)','auto','$\overline{S}$ (\%)')
func_fig_shade(aer_tick,Y3,2,[50,50000],'$N_{a}$ (cm$^{-3}$)','auto','$\sigma_{S}$ (\%)')
func_fig_shade(aer_tick,Y1,3,[50,50000],'$N_{a}$ (cm$^{-3}$)','auto','$\sigma_{S}/\overline{S}$')

% print('-dpng',B,OutFigName,'-r450')
%%
function func_fig_shade(xdat,ydat,i_p,xlmt,xlbl,ylmt,ylbl)
global pstn_list
para_xylbl = {'FontSize',14,'FontWeight','bold'};
para_axis = {'linewidth',1,'FontName','Times New Roman','FontSize',12,...
    'xtick',[10,100,1000,10^4,10^5,10^6]};
clr_list = [223,122,94;227 216 183;130 178 154;60 64 91]/255;
ttwd = {'(a)';'(b)';'(c)'};
subplot('position',pstn_list(i_p,:))

plot(xdat,ydat(:,2),'-o','Color',clr_list(4,:),'LineWidth',2,'markersize',4);
y_shade=[ydat(:,3);flipud(ydat(:,1))];
x_shade = [xdat';flipud(xdat')];
hold on; fill(x_shade,y_shade,...
   'm','FaceColor',clr_list(4,:),'FaceAlpha',0.15,...
   'EdgeColor','none','handlevisibility','off');

if i_p==3
    plot(xlmt,[0.73,0.73],'--','Color',clr_list(1,:),'LineWidth',2)
    hold on
    plot(xlmt,[0.45,0.45],'--','Color',clr_list(2,:),'LineWidth',2)
    hold on
    plot(xlmt,[0.34,0.34],'--','Color',clr_list(3,:),'LineWidth',2)
    lgdwd = {'0.0$<$AF$\le$0.4';'0.4$<$AF$\le$0.7';'0.7$<$AF$\le$1.0'};
    legend(lgdwd,'interpreter','Latex','location','northeast')
    legend('box','off','AutoUpdate','off')
end
box('on')
xlim(xlmt)
ylim(ylmt)
v = axis;
if v(2)-v(1)>1000
    xscl = 'log';
else
    xscl = 'linear';
end
if v(4)-v(3)>1000
    yscl = 'log';
else
    yscl = 'linear';
end
set(gca,para_axis{:},'xscale',xscl,'yscale',yscl,...
   'XMinorGrid','off','YMinorGrid','off')
xlabel(xlbl,'interpreter','latex',para_xylbl{:})
ylabel(ylbl,'interpreter','latex',para_xylbl{:})
grid('on')
title(cell2mat(ttwd(i_p)),'unit','normalized','position',[0.102,0.876,0],'FontSize',16,'FontWeight','bold')
end

function [dat_epss,dat_sbar,dat_sigs] = func_cloudy_region_supertsaturation(paths,mnt_tick)
a=2.53E12; % coefficients for S
b=5.42E3;
r = 2.^(1:1/3:35/3); % bin radius
m = 4/3*pi.*r.*r.*r/10^9; % bin mass
nmnt = length(mnt_tick);
len_aer=length(paths);
[nx,~,nz] = size(double(ncread([cell2mat(paths(1)),'/wrfbin_d01_0001-01-01_00:00:00'],'ALT')));
[dat_epss,dat_sbar,dat_sigs]=deal(zeros(nmnt,len_aer));
for ia = 1:len_aer
    ia
    im = 0;
    for mnt = mnt_tick
        im = im+1;
        ncpath = [cell2mat(paths(ia)),'/wrfbin_d01_0001-01-01_0',...
              num2str(floor(mnt/60),'%01d'),':',num2str(mod(mnt,60),'%02d'),':00'];
        [nd,qc] = deal(zeros(nx,nx,nz));
        for ibin=1:33
            qc_bin=double(ncread(ncpath,['ff1i',num2str(ibin,'%02d')]));
            qc = qc+qc_bin;
            nd=nd+qc_bin./m(ibin);
        end
        qv = double(ncread(ncpath,'QVAPOR'));
        P = double(ncread(ncpath,'P'))+double(ncread(ncpath,'PB'));
        TH=double(ncread(ncpath,'T'))+300;
        T=TH.*(P./10^5).^0.286;
        es=a/10*exp(-b/T);
        ew=qv.*P./(0.622+0.378*qv);
        s = ew./es*100-100;
        rho=1/double(ncread(ncpath,'ALT'));% m3/kg
        nd = nd.*rho;
        icloud = qc>10^-5; 
        s = s(icloud);
        nd = nd(icloud);
        sbar = sum(s.*nd)/sum(nd);
        sigs = sqrt(sum((s-sbar).^2.*nd)/sum(nd));
        epss = sigs/sbar;
        dat_epss(im,ia)=epss;
        dat_sbar(im,ia)=sbar;
        dat_sigs(im,ia)=sigs;
    end
end
end