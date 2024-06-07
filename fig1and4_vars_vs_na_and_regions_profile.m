clear; clc; close all
paths = {'~/WRFV4.5.1/a50';'~/WRFV4.5.1/a100';'~/WRFV4.5.1/a200';'~/WRFV4.5.1/a500';...
'~/WRFV4.5.1/a1000';'~/WRFV4.5.1/a2000';'~/WRFV4.5.1/a5000';...
'~/WRFV4.5.1/a10000';'~/WRFV4.5.1/a20000';'~/WRFV4.5.1/a50000'};
aer_tick= [50,100,200,500,1000,2000,5000,10000,20000,50000];
aer_id = [2,7,10];
expname = 'exp1';

global af_lmt pstn_list
mnt_interval = 3;
mnt_tick = 75:mnt_interval:135;
nmnt = length(mnt_tick);
InputFileName = ['lwcad_',expname,'_75to135.mat'];
InterFileName = 'vars_vs_na_and_regions_profile.mat';
OutFigName = 'nc_rbar_eps_na';
len_aer = length(aer_tick);
af_lmt=[0,0.4,0.7];
nregion = length(af_lmt);
pstn_list = [0.090    0.7250    0.241    0.236;
   0.4068    0.7250    0.241    0.236;
   0.7236    0.7250    0.241    0.236;
   0.090    0.4000    0.241    0.236;
   0.4068    0.4000    0.241    0.236;
   0.7236    0.4000    0.241    0.236;
   0.090    0.0840    0.241    0.236;
   0.4068    0.0840    0.241    0.236;
   0.7236    0.0840    0.241    0.236];
%% prepare data for figure
g=9.81;
fstpath = [cell2mat(paths(1)),'/wrfbin_d01_0001-01-01_01:54:00'];
ph = double(ncread(fstpath,'PHB'))+double(ncread(fstpath,'PH')); 
zz = squeeze(mean(mean(ph(:,:,2:end)+ph(:,:,1:end-1))/2)/g/1000); % unit: km
load(InputFileName)
[dat_rm,dat_eps,dat_nd,dat_epsa,dat_epsb] =deal(nan(nmnt,nregion,len_aer));
for ia = 1:len_aer
    ia
    [dat_eps(:,:,ia),dat_rm(:,:,ia),dat_nd(:,:,ia)] ...
        = func_vars_vs_na(paths(ia),mnt_tick,ia,lwcad);
end

[epsilon1,rbar1,~] = func_regions_profiles(paths(aer_id(1)),mnt_tick,aer_id(1),lwcad);
[epsilon2,rbar2,~] = func_regions_profiles(paths(aer_id(2)),mnt_tick,aer_id(2),lwcad);
[epsilon3,rbar3,~] = func_regions_profiles(paths(aer_id(3)),mnt_tick,aer_id(3),lwcad);
save(InterFileName,'epsilon1','epsilon2','epsilon3','rbar1','rbar2','rbar3','zz','dat_eps','dat_rm','dat_nd')
%% Figure 1
B=figure('position',[488,100.2,765.8,625.8]);
load(InterFileName)
[Y1,Y2,Y3] = deal(zeros(len_aer,3,nregion));
for ir=1:nregion
    for ia = 1:len_aer
        dat = squeeze(dat_nd(:,ir,:));
        Y1(ia,:,ir) = quantile(dat(:,ia),[0.25,0.5,0.75]);
        
        dat = squeeze(dat_rm(:,ir,:));
        Y2(ia,:,ir) = quantile(dat(:,ia),[0.25,0.5,0.75]);
        
        dat = squeeze(dat_eps(:,ir,:));
        Y3(ia,:,ir) = quantile(dat(:,ia),[0.25,0.5,0.75]);
    end
end
func_fig_shade(aer_tick(1:10),Y1(1:end,:,:),1,[50,50000],'$N_{a}$ (cm$^{-3}$)','$N_{c}$ (cm$^{-3}$)')
lgwd = {'AF$>$0';'AF$>$0.4';'AF$>$0.7'}; % legend content
legend(lgwd,'position',[0.1886,0.7223,0.1329,0.0924],'fontsize',11,'interpreter','latex')
legend('boxoff')
func_fig_shade(aer_tick(1:10),Y2(1:end,:,:),2,[50,50000],'$N_{a}$ (cm$^{-3}$)','$\overline{r} (\mu$m)')
func_fig_shade(aer_tick(1:10),Y3(1:end,:,:),3,[50,50000],'$N_{a}$ (cm$^{-3}$)','ε')

%-------------figure settings------------------%
para_xylbl = {'FontSize',14,'interpreter','latex'};
para_axis = {'linewidth',1,'FontName','Times New Roman',...
    'FontSize',12};
xlb = {'$\overline{r} (\mu$m)';'ε'};
var = {'rbar';'epsilon'};
ttwd = {'AF>0.0';'AF>0.4';'AF>0.7'};
nbwd = 'abcdefghijklmn';
for irow = 2:3
for icol = 1:3
    i_p=icol+3*irow-3;
   ax = subplot('position',pstn_list(i_p,:));
   for ia = 1:length(aer_id)
       eval(['dat = ',cell2mat(var(irow-1)),num2str(ia),';'])
       hold on;
       h = plot(dat(:,2,icol),zz);
       set(h,'Linewidth',1.5)
       clr = get(h,'color');
       hold on;
       x=[dat(:,3,icol);flipud(dat(:,1,icol))];
       y=[zz;flipud(zz)];
       sh = fill(x(~isnan(x)),y(~isnan(x)),...
           'm','FaceColor',clr,'FaceAlpha',0.15,...
           'EdgeColor','none','handlevisibility','off');
   end
   if irow ==3
       set(ax,para_axis{:},'xtick',0:0.05:1)
       if icol==1
           xlim([0.18,0.44])
       elseif icol==2
           xlim([0.15,0.41])
       elseif icol==3
           xlim([0.09,0.35])
       end
   else
       set(ax,para_axis{:})
   end
   ylabel('Altitude (km)',para_xylbl{:})
   xlabel(cell2mat(xlb(irow-1)),para_xylbl{:})
   box('on')
   grid('on')
   ylim([0.5,3.1])
   title(['(',nbwd(i_p),') ',cell2mat(ttwd(icol))],'fontsize',15,...
       'units','normalized','position',[0.276,0.833,0],'FontWeight','bold');
   if i_p==4
       lgwd = {'100 cm$^{-3}$';'5000 cm$^{-3}$';'50000 cm$^{-3}$'};
        legend(lgwd,'interpreter','latex','units','normalized','position',[0.174634943849569,0.514924544582934,0.1547,0.0847],'FontSize',10)
        legend(ax,'boxoff')
   end
end
end
% print('-dpng',B,'vars_vs_na_and_epsilon_profile','-r450')
%% Figure 4
 C=figure('position',[488,414.6,380,310]);
xdat = [squeeze(Y2(:,2,1))';
    squeeze(Y2(:,2,2))';
    squeeze(Y2(:,2,3))'];
func_fig_shade(xdat,Y3,1,[3,21],'$\overline{r} (\mu$m)','ε')
set(gca,'position',[0.156,0.164,0.737,0.75],'xtick',0:3:30)
ylim([0.15,0.37])
lgwd = {'AF$>$0';'AF$>$0.4';'AF$>$0.7'}; % legend content
legend(lgwd,'fontsize',12,'interpreter','latex','location','northwest')
legend('boxoff')
% print('-dpng',C,'eps_rbar','-r450')

%%
function func_fig_shade(xdat,ydat,i_p,xlmt,xlbl,ylbl)
global pstn_list
para_xylbl = {'FontSize',14,'interpreter','latex'};
para_axis = {'linewidth',1,'FontName','Times New Roman','FontSize',12,'xtick',[10,100,1000,10^4,10^5,10^6]};
nbwd = 'abcdefghijklmn';
clr_list = [60 64 91;223,122,94;130 178 154]/255;
subplot('position',pstn_list(i_p,:))
[~,~,n_line]=size(ydat);
[n_x,~]=size(xdat);
for i_l = 1:n_line
    dati = squeeze(ydat(:,:,i_l));
    if n_x==1
        x=xdat;
    elseif n_x==n_line
        x=xdat(i_l,:);
    end
hold on; plot(x,dati(:,2),'-o','Color',clr_list(i_l,:),'LineWidth',2,'markersize',4);
y_shade=[dati(:,3);flipud(dati(:,1))];
x_shade = [x';flipud(x')];
hold on; fill(x_shade,y_shade,...
   'm','FaceColor',clr_list(i_l,:),'FaceAlpha',0.15,...
   'EdgeColor','none','handlevisibility','off');
end
box('on')
xlim(xlmt)
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
xlabel(xlbl,para_xylbl{:})
ylabel(ylbl,para_xylbl{:})
grid('on')
title(['(',nbwd(i_p),')'],'fontsize',15,'FontWeight','bold','unit','normalized','position',[0.102,0.833,0])
end

function [epsilon_dat,rbar_dat,sigma_dat] =func_regions_profiles(path1,mnt_tick,ia,lwcad)
global af_lmt
nmnt = length(mnt_tick);
mnt_interval = mnt_tick(2)-mnt_tick(1);
[nx,~,nz] = size(double(ncread([cell2mat(path1),'/wrfbin_d01_0001-01-01_00:00:00'],'ALT')));
im = 0;
nregion = length(af_lmt);
[eps_profile,r_profile,sig_profile] = deal(nan(nmnt,nz,nregion)); 
[epsilon_dat,rbar_dat,sigma_dat] = deal(nan(nz,3,nregion));
for mnt = mnt_tick
    im = im+1;
    pathbin = [cell2mat(path1),'/wrfbin_d01_0001-01-01_0',...
        num2str(floor(mnt/60),'%01d'),':',num2str(mod(mnt,60),'%02d'),':00'];
    [qc0,eps0,r0,sig0] = func_get_pro(pathbin,nx,nz);
    %---------------lwcmax---------------%
%     af = qc./repmat(max(max(qc,[],1),[],2),nx,nx);
    %---------------lwcad----------------%
    rho=1/double(ncread(pathbin,'ALT'));% m3/kg
    lwc = qc0.*rho;
    af = lwc./repmat(permute(lwcad(:,(mnt-mnt_tick(1)+mnt_interval)/mnt_interval,ia),[3,2,1]),nx,nx);
    for ir = 1:nregion
        icloud0 = qc0>10^-5 & af>af_lmt(ir);
        eps_profile(im,:,ir) = func_HorizonAve(eps0,icloud0);
        r_profile(im,:,ir) = func_HorizonAve(r0,icloud0);
        sig_profile(im,:,ir) = func_HorizonAve(sig0,icloud0);
    end
end
for iz=1:nz
    for ir=1:nregion
        eps_mnt=eps_profile(:,iz,ir);
        r_mnt=r_profile(:,iz,ir);
        sig_mnt=sig_profile(:,iz,ir);
        if length(eps_mnt(~isnan(eps_mnt)))<=5
            continue
        end
        epsilon_dat(iz,:,ir) = quantile(eps_mnt(~isnan(eps_mnt)),[0.25,0.5,0.75]); 
        rbar_dat(iz,:,ir) = quantile(r_mnt(~isnan(r_mnt)),[0.25,0.5,0.75]); 
        sigma_dat(iz,:,ir) = quantile(sig_mnt(~isnan(sig_mnt)),[0.25,0.5,0.75]); 
    end
end
end

function [dat_eps,dat_rm,dat_nd] =func_vars_vs_na(path1,mnt_tick,ia,lwcad)
global af_lmt
g=9.81;
nregion = length(af_lmt);
nmnt = length(mnt_tick);
mnt_interval = mnt_tick(2)-mnt_tick(1);
[nx,~,nz] = size(double(ncread([cell2mat(path1),'/wrfbin_d01_0001-01-01_00:00:00'],'ALT')));
[dat_nd,dat_rm,dat_eps]=deal(nan(nmnt,nregion));
im = 0;
for mnt = mnt_tick
    im = im+1;
    ncpath = [cell2mat(path1),'/wrfbin_d01_0001-01-01_0',...
        num2str(floor(mnt/60),'%01d'),':',num2str(mod(mnt,60),'%02d'),':00'];
    ph = double(ncread(ncpath,'PHB'))+double(ncread(ncpath,'PH'));
    h0 = squeeze(mean(mean(ph(:,:,2:end)-ph(:,:,1:end-1)))/g); % unit: m
    [qc,eps0,rm0,~,nd0] = func_get_pro(ncpath,nx,nz);
    
    %---------------lwcmax---------------%
%     af = qc./repmat(max(max(qc,[],1),[],2),nx,nx);
    %---------------lwcad----------------%
    rho=1/double(ncread(ncpath,'ALT'));% m3/kg
    nd0 = nd0.*rho;
    lwc = qc.*rho;
    af = lwc./repmat(permute(lwcad(:,(mnt-mnt_tick(1)+mnt_interval)/mnt_interval,ia),[3,2,1]),nx,nx);
    for ir = 1:3
        icloud0 = qc>10^-5 & af>af_lmt(ir); 
        dat_eps(im,ir) = func_ProfileAve(eps0,icloud0,h0);
        dat_rm(im,ir) = func_ProfileAve(rm0,icloud0,h0);
        dat_nd(im,ir) = func_ProfileAve(nd0,icloud0,h0);
    end
end
end

function rsl = func_ProfileAve(var3dim,icloud0,h0)
    profile = func_HorizonAve(var3dim,icloud0);
    h = h0;
    h(isnan(profile)) = [];
    profile(isnan(profile)) = [];
    rsl = sum(h.*profile)./sum(h);
end

function rsl = func_HorizonAve(var3dim,icloud0)
    var3dim(~icloud0) = 0;
    numicloud = squeeze(sum(sum(icloud0)));
    rsl = squeeze(sum(sum(var3dim)))./numicloud;
end

function [qc,eps,rm,sig,nd] = func_get_pro(ncpath,nx,nz)
r = 2.^(1:1/3:35/3); % bin radius
m = 4/3*pi.*r.*r.*r/10^9; % bin mass
nbins = 33;
[qc,nd,rm,sig] = deal(zeros(nx,nx,nz));
for ibin=1:nbins
    q=double(ncread(ncpath,['ff1i',num2str(ibin,'%02d')]));
    qc = qc+q;
    npdr=q./m(ibin); % 10^6 kg-1
    rm = rm+npdr.*r(ibin);
    nd = nd+npdr;
end
rm = rm./nd;
icloud00 = qc>10^-5;
for ibin=1:nbins
    q=double(ncread(ncpath,['ff1i',num2str(ibin,'%02d')]));
    q(~icloud00) = nan;
    npdr=q./m(ibin);
    sig = sig+npdr.*(rm-r(ibin)).^2;
end
sig = sqrt(sig./nd);
eps = sig./rm;
end
