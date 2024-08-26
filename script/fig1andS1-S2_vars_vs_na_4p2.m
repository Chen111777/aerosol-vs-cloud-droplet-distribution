clear; clc; close all
paths = {'~/WRFV4.5.1/a50_exp2';'~/WRFV4.5.1/a500_exp2';...
'~/WRFV4.5.1/a2000_exp2';'~/WRFV4.5.1/a5000_exp2';...
'~/WRFV4.5.1/a10000_exp2';'~/WRFV4.5.1/a50000_exp2'}; % for fig. S1-S2
% paths = {'~/WRFV4.5.1/a50';'~/WRFV4.5.1/a100';'~/WRFV4.5.1/a200';'~/WRFV4.5.1/a500';...
% '~/WRFV4.5.1/a1000';'~/WRFV4.5.1/a2000';'~/WRFV4.5.1/a5000';...
% '~/WRFV4.5.1/a10000';'~/WRFV4.5.1/a20000';'~/WRFV4.5.1/a50000'}; % for fig. 1
aer_tick= [50,500,2000,5000,10000,50000]; % for fig. S1-S2
% aer_tick= [50,100,200,500,1000,2000,5000,10000,20000,50000]; % for fig. 1
aer_id = [1,4,6]; % for fig. S1-S2
% aer_id = [1,7,10]; % for fig. 1
expname = 'exp2';

global af_lmt pstn_list
mnt_interval = 3;
mnt_tick = 75:mnt_interval:135;
nmnt = length(mnt_tick);
InputFileName = ['lwcad_',expname,'_75to135.mat'];
InterFileName = ['vars_vs_na_4p2_',expname,'.mat'];
OutFigName = ['vars_vs_na_4p2_',expname];
len_aer = length(aer_tick);
af_lmt=[0,0.4,0.7,1];
nregion = 4;
%% prepare data for figure
g=9.81;
path_example = [cell2mat(paths(1)),'/wrfbin_d01_0001-01-01_01:54:00'];
ph = double(ncread(path_example,'PHB'))+double(ncread(path_example,'PH')); 
zz = squeeze(mean(mean(ph(:,:,2:end)+ph(:,:,1:end-1))/2)/g/1000); % unit: km
load(InputFileName)
[xdat1,xdat2,xdat3,xdat4]=deal(zeros(nmnt,length(zz),length(aer_id)));
[ydat1,ydat2,ydat3,ydat4]=deal(zeros(nmnt,len_aer,nregion));
for ia = 1:len_aer
    ia
    [dat_eps,dat_rbar,dat_nd,dat_sig] ...
        = func_af_range_average(paths(ia),mnt_tick,ia,lwcad);
    for iline=1:nregion
        ydat1(:,ia,iline)=squeeze(dat_nd(:,:,iline))';
        ydat2(:,ia,iline)=squeeze(dat_rbar(:,:,iline))';
        ydat3(:,ia,iline)=squeeze(dat_eps(:,:,iline))';
        ydat4(:,ia,iline)=squeeze(dat_sig(:,:,iline))';
    end
end

for i_id=1:length(aer_id)
    i_id
    [dat_eps,dat_rbar,dat_nd,dat_sig] ...
        = func_cloudy_region_profile(paths(aer_id(i_id)),mnt_tick);
    xdat1(:,:,i_id)=dat_nd';
    xdat2(:,:,i_id)=dat_rbar';
    xdat3(:,:,i_id)=dat_eps';
    xdat4(:,:,i_id)=dat_sig';
end

save(InterFileName,'ydat1','ydat2','ydat3','ydat4','xdat1','xdat2','xdat3','xdat4','zz')
%% Figure 1
load(InterFileName)
pstn_list = [0.105,0.805,0.345,0.175;
    0.575,0.805,0.345,0.175;
    0.105,0.56,0.345,0.175;
    0.575,0.56,0.345,0.175;
    0.105,0.315,0.345,0.175;
    0.575,0.315,0.345,0.175;
    0.105,0.07,0.345,0.175;
    0.575,0.07,0.345,0.175];
B=figure('position',[571.4,49.8,490,733.6]);

%-------------right column------------------%
func_fig_shade(aer_tick,ydat1,2,[50,50000],'$N_{a}$ (cm$^{-3}$)','$N_{c}$ (cm$^{-3}$)')
ylim([1,10^4])
lgwd = {'0.0$<$AF$\le$0.4';'0.4$<$AF$\le$0.7';'0.7$<$AF$\le$1.0';'0.0$<$AF$\le$1.0'};
legend(lgwd,'fontsize',8.5,'position',[0.676,0.801,0.241,0.081],'interpreter','latex')
legend('boxoff')
func_fig_shade(aer_tick,ydat2,4,[50,50000],'$N_{a}$ (cm$^{-3}$)','$\overline{r} (\mu$m)')
func_fig_shade(aer_tick,ydat3,6,[50,50000],'$N_{a}$ (cm$^{-3}$)','¦Å')
ylim([0.14,0.4])%[0.12,0.35])
func_fig_shade(aer_tick,ydat4,8,[50,50000],'$N_{a}$ (cm$^{-3}$)','$\sigma (\mu$m)')


%-------------left column------------------%
para_xylbl = {'FontSize',12,'interpreter','latex'};
para_axis = {'linewidth',1,'FontName','Times New Roman',...
    'FontSize',10};
xlb = {'$N_{c}$ (cm$^{-3})$';'$\overline{r} (\mu$m)';'¦Å';'$\sigma (\mu$m)'};
nbwd = 'abcdefghijklmn';
icol = 1;
for irow = 1:4
    i_p=icol+2*irow-2;
   ax = subplot('position',pstn_list(i_p,:));
   eval(['xdat = xdat',num2str(irow),';'])
   for i_line = 1:length(aer_id)
       dati = squeeze(xdat(:,:,i_line));
       hold on;
       dati(:,sum(~isnan(dati))<=5)=nan;
       h = plot(nanmean(dati),zz);
       set(h,'Linewidth',1.5)
       clr = get(h,'color');
       hold on;
       x_shade=[quantile(dati,0.75)';flipud(quantile(dati,0.25)')];
       y_shade=[zz;flipud(zz)];
       sh = fill(x_shade(~isnan(x_shade)),y_shade(~isnan(x_shade)),...
           'm','FaceColor',clr,'FaceAlpha',0.15,...
           'EdgeColor','none','handlevisibility','off');
   end
   set(ax,para_axis{:})
   ylabel('Altitude (km)',para_xylbl{:})
   xlabel(cell2mat(xlb(irow)),para_xylbl{:})
   box('on')
   grid('on')
   ylim([0.7,3.1])%2.3
   title(['(',nbwd(i_p),') '],'fontsize',15,...
       'units','normalized','position',[0.102,0.833,0],'FontWeight','bold');
   if i_p==1
        set(gca, 'XMinorGrid','off', 'xscale','log')
   end
   if i_p==5
        lgwd = {'50 cm$^{-3}$';'5,000 cm$^{-3}$';'50,000 cm$^{-3}$'};
        legend(lgwd,'position',[0.208,0.335,0.244,0.065],'interpreter','latex',...
            'units','normalized','FontSize',9.0)
        legend(ax,'boxoff')
   end
end
% print('-dpng',B,OutFigName,'-r450')
%%
function func_fig_shade(xdat,ydat,i_p,xlmt,xlbl,ylbl)
global pstn_list
para_xylbl = {'FontSize',12,'interpreter','latex'};
para_axis = {'linewidth',1,'FontName','Times New Roman','FontSize',10,'xtick',[10,100,1000,10^4,10^5,10^6]};
nbwd = 'abcdefghijklmn';
clr_list = [223,122,94;227 216 183;130 178 154;60 64 91]/255;
subplot('position',pstn_list(i_p,:))
[~,~,n_line]=size(ydat);
[n_xdat,~]=size(xdat);
for i_line = 1:n_line
    dati = squeeze(ydat(:,:,i_line));
    if n_xdat==1
        x=xdat;
    elseif n_xdat==n_line
        x=xdat(i_line,:);
    end
    y=nanmean(dati);
    hold on; plot(x,y,'-o','Color',clr_list(i_line,:),'LineWidth',2,'markersize',4);
    y_shade=[quantile(dati,0.75)';flipud(quantile(dati,0.25)')];
    x_shade = [x';flipud(x')];
    hold on; fill(x_shade,y_shade,...
       'm','FaceColor',clr_list(i_line,:),'FaceAlpha',0.15,...
       'EdgeColor','none','handlevisibility','off');
end
box('on')
xlim(xlmt)
ax_range = axis;
[xscl,yscl] = deal('linear');
if ax_range(2)-ax_range(1)>1000
    xscl = 'log';
end
if ax_range(4)-ax_range(3)>1000
    yscl = 'log';
end
set(gca,para_axis{:},'xscale',xscl,'yscale',yscl,...
   'XMinorGrid','off','YMinorGrid','off')
xlabel(xlbl,para_xylbl{:})
ylabel(ylbl,para_xylbl{:})
grid('on')
title(['(',nbwd(i_p),')'],'fontsize',15,'FontWeight','bold','unit','normalized','position',[0.102,0.833,0])
end

function [dat_eps,dat_rbar,dat_nd,dat_sig] = func_af_range_average(path1,mnt_tick,ia,lwcad)
global af_lmt
g=9.81;
nregion = 4;
nmnt = length(mnt_tick);
mnt_interval = mnt_tick(2)-mnt_tick(1);
[nx,~,nz] = size(double(ncread([cell2mat(path1),'/wrfbin_d01_0001-01-01_00:00:00'],'ALT')));
[dat_nd,dat_rbar,dat_eps]=deal(nan(nmnt,1,nregion));
im = 0;
for mnt = mnt_tick
    im = im+1;
    ncfile = [cell2mat(path1),'/wrfbin_d01_0001-01-01_0',...
        num2str(floor(mnt/60),'%01d'),':',num2str(mod(mnt,60),'%02d'),':00'];
    [qc,eps0,rbar0,sig0,nd0] = func_get_grid_properties(ncfile,nx,nz);

    ph = double(ncread(ncfile,'PHB'))+double(ncread(ncfile,'PH'));
    z0 = (ph(:,:,2:end)+ph(:,:,1:end-1))/2000/g; % unit: km
    rho=1/double(ncread(ncfile,'ALT'));% kg/m3, or 10^-6 kg/cm3
    nd0 = nd0.*rho;
    lwc = qc.*rho;
    af = lwc./repmat(permute(lwcad(:,(mnt-mnt_tick(1)+mnt_interval)/mnt_interval,ia),[3,2,1]),nx,nx);
    izb = find(any(any(qc>10^-5,1),2),1)-1; % cloud base index
    hacb=z0-repmat(z0(:,:,izb),1,1,nz); % height above cloud base
    for ir = 1:nregion
        if ir==nregion
            icloud = qc>10^-5 & hacb>0.2;
        else
            icloud = qc>10^-5 & af>af_lmt(ir) & af<=af_lmt(ir+1) & hacb>0.2;
        end
        dat_eps(im,1,ir) = mean(eps0(icloud));
        dat_rbar(im,1,ir) = mean(rbar0(icloud));
        dat_nd(im,1,ir) = mean(nd0(icloud));
        dat_sig(im,1,ir) = mean(sig0(icloud));
    end
end
end

function [qc,eps,rbar,sig,nd] = func_get_grid_properties(ncfile,nx,nz)
r = 2.^(1:1/3:35/3); % bin radius
m = 4/3*pi.*r.*r.*r/10^9; % bin mass
nbins = 33;
[qc,nd,rbar,sig] = deal(zeros(nx,nx,nz));
for ibin=1:nbins
    qc_bin=double(ncread(ncfile,['ff1i',num2str(ibin,'%02d')]));
    qc = qc+qc_bin;
    nd_bin=qc_bin./m(ibin); % 10^6 kg-1
    rbar = rbar+nd_bin.*r(ibin);
    nd = nd+nd_bin;
end
rbar = rbar./nd;
icloud = qc>10^-5;
for ibin=1:nbins
    qc_bin=double(ncread(ncfile,['ff1i',num2str(ibin,'%02d')]));
    qc_bin(~icloud) = nan;
    nd_bin=qc_bin./m(ibin);
    sig = sig+nd_bin.*(rbar-r(ibin)).^2;
end
sig = sqrt(sig./nd);
eps = sig./rbar;
end

function [dat_eps,dat_rbar,dat_nd,dat_sig] =func_cloudy_region_profile(path1,mnt_tick)
r = 2.^(1:1/3:35/3); % bin radius
m = 4/3*pi.*r.*r.*r/10^9; % bin mass
nbins = 33;
nmnt = length(mnt_tick);
[nx,~,nz] = size(double(ncread([cell2mat(path1),'/wrfbin_d01_0001-01-01_00:00:00'],'ALT')));
[dat_nd,dat_rbar,dat_eps,dat_sig]=deal(nan(nz,nmnt));
im = 0;
for mnt = mnt_tick
    im = im+1;
    ncfile = [cell2mat(path1),'/wrfbin_d01_0001-01-01_0',...
        num2str(floor(mnt/60),'%01d'),':',num2str(mod(mnt,60),'%02d'),':00'];
    qc = double(ncread(ncfile,'QCLOUD'));
    icloud = qc>10^-5;
    
    [nd_grid,rsum_grid] = deal(zeros(size(qc)));
    rho=1/double(ncread(ncfile,'ALT'));
    for ibin=1:nbins
        qc_bin=double(ncread(ncfile,['ff1i',num2str(ibin,'%02d')]));
        nd_bin=qc_bin./m(ibin).*rho; % cm-3
        rsum_grid = rsum_grid+nd_bin.*r(ibin);
        nd_grid = nd_grid+nd_bin;
    end
    rsum_grid(~icloud)=0;
    nd_grid(~icloud)=0;
    ndsum_z=squeeze(sum(sum(nd_grid)));
    ndsum_z(ndsum_z==0)=nan;
    rbar_z = squeeze(sum(sum(rsum_grid)))./ndsum_z;
    sig_z=zeros(nz,1);
    for ibin=1:nbins
        qc_bin=double(ncread(ncfile,['ff1i',num2str(ibin,'%02d')]));
        qc_bin(~icloud) = 0;
        nd_bin=qc_bin./m(ibin).*rho;
        ndsum_bin=squeeze(sum(sum(nd_bin)));
        sig_z = sig_z+ndsum_bin.*(rbar_z-r(ibin)).^2;
    end
    sig_z = sqrt(sig_z./ndsum_z);
    eps_z = sig_z./rbar_z;
    nd_z=ndsum_z./squeeze(sum(sum(icloud)));

    dat_eps(:,im) = eps_z;
    dat_rbar(:,im) = rbar_z;
    dat_nd(:,im) = nd_z;
    dat_sig(:,im) = sig_z;
    end
end
end
