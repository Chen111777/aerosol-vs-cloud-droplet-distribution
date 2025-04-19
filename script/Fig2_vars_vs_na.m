clear; clc; close all
% paths = {'~/submit2_wrf4.5.1/a50_';...
% '~/submit2_wrf4.5.1/a100_';...
% '~/submit2_wrf4.5.1/a200_';...
% '~/submit2_wrf4.5.1/a500_';...
% '~/submit2_wrf4.5.1/a1000_';...
% '~/submit2_wrf4.5.1/a2000_';...
% '~/submit2_wrf4.5.1/a5000_';...
% '~/submit2_wrf4.5.1/a10000_';...
% '~/submit2_wrf4.5.1/a20000_';...
% '~/submit2_wrf4.5.1/a50000_'};
paths = {'~/submit2_wrf4.5.1/a50_';...
'~/submit2_wrf4.5.1/a1000_';...
'~/submit2_wrf4.5.1/a5000_';...
'~/submit2_wrf4.5.1/a50000_'};
% aer_tick= [50,100,200,500,1000,2000,5000,10000,20000,50000]; % for fig. 1
aer_tick= [50,1000,5000,50000]; % for fig. S1 and S2
% aer_id = [1,7,10]; % for fig. 1
aer_id = [1,3,4]; % for fig. S1 and S2


global af_lmt pstn_list expname
expname = 'sma100';
mnt_interval = 4;
mnt_tick = 50:mnt_interval:134;
nmnt = length(mnt_tick);
lwcad_file = ['lwcad_',expname,'.mat'];
InterFileName = ['vars_vs_na_4p2_',expname,'.mat'];
OutFigName = ['vars_vs_na_4p2_',expname];
len_aer = length(aer_tick);
af_lmt=[0,0.4,0.7,1];
nregion = 4;
%% prepare data for figure
g=9.81;
path_example = [cell2mat(paths(1)),expname,'/wrfbin_d01_0001-01-01_01:54:00'];
ph = double(ncread(path_example,'PHB'))+double(ncread(path_example,'PH')); 
zz = squeeze(mean(mean(ph(:,:,2:end)+ph(:,:,1:end-1))/2)/g/1000); % unit: km
load(lwcad_file)
[xdat1,xdat2,xdat3,xdat4]=deal(zeros(nmnt,length(zz),length(aer_id)));
[ydat1,ydat2,ydat3,ydat4]=deal(zeros(nmnt,len_aer,nregion));
for ia = 1:len_aer
    ia
    [dat_eps,dat_rbar,dat_nd,dat_sig] ...
        = func_af_range_average(paths(ia),mnt_tick,ia,lwcad);
    for iline=1:nregion
        ydat1(:,ia,iline)=squeeze(dat_nd(:,:,iline))';
        ydat2(:,ia,iline)=squeeze(dat_rbar(:,:,iline))';
        ydat3(:,ia,iline)=squeeze(dat_sig(:,:,iline))';
        ydat4(:,ia,iline)=squeeze(dat_eps(:,:,iline))';
    end
end

for i_id=1:length(aer_id)
    i_id
    [dat_eps,dat_rbar,dat_nd,dat_sig] ...
        = func_cloudy_region_profile(paths(aer_id(i_id)),mnt_tick);
    xdat1(:,:,i_id)=dat_nd';
    xdat2(:,:,i_id)=dat_rbar';
    xdat3(:,:,i_id)=dat_sig';
    xdat4(:,:,i_id)=dat_eps';
end

save(InterFileName,'ydat1','ydat2','ydat3','ydat4','xdat1','xdat2','xdat3','xdat4','zz','mnt_tick')
%% Figure 2
load(InterFileName)
pstn_list = [0.120,0.805,0.310,0.160;
    0.580,0.805,0.310,0.160;
    0.120,0.56,0.310,0.160;
    0.580,0.56,0.310,0.160;
    0.120,0.315,0.310,0.160;
    0.580,0.315,0.310,0.160;
    0.120,0.07,0.310,0.160;
    0.580,0.07,0.310,0.160];
B=figure('position',[571.4,49.8,500,733.6]);

%-------------right column------------------%
x_draw1 = aer_tick;
x_draw2 = ydat1;
x_draw2=nanmean(x_draw2(:,:,end));

ax1 = func_fig_shade(x_draw1,[],ydat1,2,'$N_{a}$(cm$^{-3}$)','$N_{c}$(cm$^{-3}$)');
set(ax1,'YLim',[1,10^4],'ytick',10.^[0:5])
lgwd = {'0.0$<$AF$\le$0.4';'0.4$<$AF$\le$0.7';'0.7$<$AF$\le$1.0';'0.0$<$AF$\le$1.0'};
legend(lgwd,'fontsize',6.5,...
    'position',[0.690,0.805,0.201,0.071],'interpreter','latex')
legend('boxoff')

func_fig_shade(x_draw1, x_draw2, ydat2,4,'$N_{a}$(cm$^{-3}$)','$\overline{r} (\mu$m)');

func_fig_shade(x_draw1, x_draw2, ydat3,6,'$N_{a}$(cm$^{-3}$)','$\sigma (\mu$m)');

ax1 = func_fig_shade(x_draw1, x_draw2, ydat4,8,'$N_{a}$(cm$^{-3}$)','  ');
set(ax1,'YLim',[0.12,0.39])


%-------------left column------------------%
para_xylbl = {'FontSize',11.5,'interpreter','latex'};
para_axis = {'linewidth',1,'FontName','Times New Roman',...
    'FontSize',9.3};
xlb = {'$N_{c}$(cm$^{-3})$';'$\overline{r}(\mu$m)';'$\sigma(\mu$m)';'  '};
nbwd = 'abcdefghijklmn';
% LineStyleList = {'-','--',':'};
ColorList = [0, 0, 128;
    70, 130, 180;
    255, 0, 255]/255;
icol = 1;
for irow = 1:4
    i_p=icol+2*irow-2;
   ax = subplot('position',pstn_list(i_p,:));
   eval(['xdat = xdat',num2str(irow),';'])
   for i_line = 1:length(aer_id)
       dati = squeeze(xdat(:,:,i_line));
       hold on;
%        dati(:,sum(~isnan(dati))<=5)=nan;
       h = plot(nanmean(dati),zz,'Color','k','Color',ColorList(i_line,:));
       set(h,'Linewidth',1.5)
       clr = get(h,'color');
       hold on;
       x_shade=[quantile(dati,0.75)';flipud(quantile(dati,0.25)')];
       y_shade=[zz;flipud(zz)];
       sh = fill(x_shade(~isnan(x_shade)),y_shade(~isnan(x_shade)),...
           'm','FaceColor',clr,'FaceAlpha',0.15,...
           'EdgeColor','none','handlevisibility','off');
   end
   set(ax,para_axis{:},'ytick',0:0.5:5)
   ylabel('Altitude (km)',para_xylbl{:})
   xlabel(cell2mat(xlb(irow)),para_xylbl{:})
   box('on')
%    grid('on')
   ylim([0.6,3.5])
   title(['(',nbwd(i_p),') '],'fontsize',15,...
       'units','normalized','position',[-0.204,1.006,0]);
   if i_p==1
        set(gca, 'XMinorGrid','off', 'xscale','log','xtick',10.^[0:5])
   end
   if i_p == 5
       set(gca, 'xtick', 0:3:15, 'xlim', [0,13])
   end
   if i_p==7
       xlim([0.15,0.52])
   end
end
lgwd = {'N_a=50';'N_a=5,000';'N_a=50,000'};
legend(lgwd,...
    'position',[0.190,0.656,0.178,0.062],...
    'units','normalized','FontSize',6.5)
legend(ax,'boxoff')
% print('-dpng',B,OutFigName,'-r450')
%%
function ax1 = func_fig_shade(x_down,x_up,ydat,i_p,xlbl,ylbl)
global pstn_list
para_xylbl = {'FontSize',11.5,'interpreter','latex'};
para_axis = {'linewidth',1,'FontName','Times New Roman','FontSize',9.3};
nbwd = 'abcdefghijklmn';
clr_list = [223,122,94;227 216 183;130 178 154;60 64 91]/255;
subplot('position',pstn_list(i_p,:))
[~,~,n_line]=size(ydat);
x = x_down;
for i_line = 1:n_line
    dati = squeeze(ydat(:,:,i_line));
    y=nanmean(dati);
    hold on; plot(x,y,'-o','Color',clr_list(i_line,:),'LineWidth',1.7,'markersize',3.4);
    y_shade=[quantile(dati,0.75)';flipud(quantile(dati,0.25)')];
    x_shade = [x';flipud(x')];
    hold on; fill(x_shade,y_shade,...
       'm','FaceColor',clr_list(i_line,:),'FaceAlpha',0.15,...
       'EdgeColor','none','handlevisibility','off');
end
ax1 = gca;
box('on')
ax_range = axis;
[xscl,yscl] = deal('linear');
if ax_range(2)-ax_range(1)>1000
    xscl = 'log';
end
if ax_range(4)-ax_range(3)>1000
    yscl = 'log';
end
set(ax1,para_axis{:},'xscale',xscl,'yscale',yscl,...
   'XTick',10.^[0:5])
ylabel(ax1,ylbl,para_xylbl{:})
xlabel(ax1,xlbl,para_xylbl{:})
ax1.XLim = [40,60000];
title(['(',nbwd(i_p),')'],'fontsize',15,...
    'unit','normalized','position',[-0.204,1.006,0])

if ~isempty(x_up)
    ax2 = axes('Position', ax1.Position);
    ax2.XLim = ax1.XLim;
    xtklb = arrayfun(@(x) num2str(x, '%.0f'), x_up, 'UniformOutput', false);
%     xtklb([2,4,6,8,10]) = {''};
    set(ax2,para_axis{:}, 'XAxisLocation', 'top', ...
               'YAxisLocation', 'left', 'Color', 'none', 'xscale', 'log', ...
               'LineWidth', 1, ...
               'YTick', [], 'XTick', x_down, ...
               'XMinorTick', 'off', 'XTickLabel', xtklb,...
               'XMinorGrid','off','FontSize',6.0,...
               'XTickLabelRotation',20);%8.0
      text(1.029680904522613,0.978219178082192,'$N_{c}$(cm$^{-3}$)','unit','normalized',...
          'interpreter','latex','FontSize',8)%0.7935,1.17
    grid('on')
end

end

function [dat_eps,dat_rbar,dat_nd,dat_sig] = func_af_range_average(path1,mnt_tick,ia,lwcad)
global af_lmt expname
g=9.81;
nregion = 4;
nmnt = length(mnt_tick);
mnt_interval = mnt_tick(2)-mnt_tick(1);
[nx,~,nz] = size(double(ncread([cell2mat(path1),expname,'/wrfbin_d01_0001-01-01_00:00:00'],'ALT')));
[dat_nd,dat_rbar,dat_eps,dat_sig]=deal(nan(nmnt,1,nregion));
im = 0;
for mnt = mnt_tick
    im = im+1;
    ncfile = [cell2mat(path1),expname,'/wrfbin_d01_0001-01-01_0',...
        num2str(floor(mnt/60),'%01d'),':',num2str(mod(mnt,60),'%02d'),':00'];
    [qc,eps0,rbar0,sig0,nd0] = func_get_grid_properties(ncfile);

    ph = double(ncread(ncfile,'PHB'))+double(ncread(ncfile,'PH'));
    z0 = (ph(:,:,2:end)+ph(:,:,1:end-1))/2000/g; % unit: km
    rho=1/double(ncread(ncfile,'ALT'));% kg/m3, or 10^-6 kg/cm3
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

function [qc,eps,rbar,sig,nd] = func_get_grid_properties(ncfile)
r = 2.^(1:1/3:35/3); % bin radius
m = 4/3*pi.*r.*r.*r/10^9; % bin mass
nbin = 33;
qc=double(ncread(ncfile,'QCLOUD'));
[nx,~,nz] = size(qc);
rho=1/double(ncread(ncfile,'ALT'));
icloud = qc>10^-5;
ncbin = zeros(nbin,nx,nx,nz);
for ibin=1:nbin
    ff1ibb = double(ncread(ncfile,['ff1i',num2str(ibin,'%02d')]));
    ncbin(ibin,:,:,:) = ff1ibb/m(ibin).*rho; % cm-3
end
nd = sum(ncbin,1);
rbar = sum(repmat(r',1,nx,nx,nz).*ncbin,1)./nd;
rbar(~icloud) = nan;
sig = sum(repmat(r'.^2,1,nx,nx,nz).*ncbin,1)./nd - rbar.^2;
sig(~icloud) = nan;
sig = sqrt(sig);
eps = sig./rbar;
end

function [dat_eps,dat_rbar,dat_nd,dat_sig] =func_cloudy_region_profile(path1,mnt_tick)
global expname
r = 2.^(1:1/3:35/3); % bin radius
m = 4/3*pi.*r.*r.*r/10^9; % bin mass
nbins = 33;
nmnt = length(mnt_tick);
[nx,~,nz] = size(double(ncread([cell2mat(path1),expname,'/wrfbin_d01_0001-01-01_00:00:00'],'ALT')));
[dat_nd,dat_rbar,dat_eps,dat_sig]=deal(nan(nz,nmnt));
im = 0;
for mnt = mnt_tick
    im = im+1;
    ncfile = [cell2mat(path1),expname,'/wrfbin_d01_0001-01-01_0',...
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
