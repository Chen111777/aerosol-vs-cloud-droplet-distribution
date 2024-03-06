clear; clc; close all
% --------- change these ------------%
% paths = {'~/WRFV4.5.1/a50_exp2';...
% '~/WRFV4.5.1/a2000_exp2';...
% '~/WRFV4.5.1/a10000_exp2';...
% '~/WRFV4.5.1/a50000_exp2'};
% aer_tick=[50,2000,10000,50000];
paths = {'~/WRFV4.5.1/a50';'~/WRFV4.5.1/a100';'~/WRFV4.5.1/a200';'~/WRFV4.5.1/a500';...
'~/WRFV4.5.1/a1000';'~/WRFV4.5.1/a2000';'~/WRFV4.5.1/a5000';...
'~/WRFV4.5.1/a10000';'~/WRFV4.5.1/a20000';'~/WRFV4.5.1/a50000'};
aer_tick= [50,100,200,500,1000,2000,5000,10000,20000,50000];
InputFileName = 'lwcad_10cases_exp1_75to135.mat';
InterFileName = 'dat_fig2_10cases_exp1.mat';
OutFigName = 'fig2_exp1';
mnt_interval = 6;
mnt_tick = 78:mnt_interval:132;

len_aer = length(aer_tick);
global m r varls 
r = 2.^(1:1/3:35/3); % bin radius
m = 4/3*pi.*r.*r.*r/10^9; % bin mass
varls = []; % variable name of bin mixing ratio
for ibin = 1:33
    varls = [varls;['ff1i',num2str(ibin,'%02d')]];
end
%%
load(InputFileName)
g=9.81;
fstpath = [cell2mat(paths(1)),'/wrfbin_d01_0001-01-01_01:54:00'];
phb = double(ncread(fstpath,'PHB'));
php = double(ncread(fstpath,'PH'));
ph = phb+php; % because ph merely varies with time and N_a
[nx,~,nz]=size(ph);
nz=nz-1;
hh = squeeze(mean(mean(ph(:,:,2:end)-ph(:,:,1:end-1)))/g);
hh=repmat(hh,1,len_aer);
clear phb php ph

% get four-region and domain-mean epsilon and r_mean profiles
[dat_eps,dat_rm] =deal(nan(nz,5,len_aer));
for ia = 1:len_aer
    ia
    [dat_eps(:,:,ia),dat_rm(:,:,ia)] ...
        = func_regions_2pro(paths(ia),mnt_tick,ia,nx,nz,lwcad);
end
save(InterFileName,'hh','dat_eps','dat_rm')
%%
load(InterFileName)
C2=[50,50,50;78.4,124.8,204;140.8,186.4,183.2;224.4,168.96,178.64;204.5,82.6,77.0]/255;
yh = squeeze(dat_eps(:,1,:));
[~,len_aer]=size(yh);
for i=1:len_aer
     [~,epsmin(i)]=min(yh(1:end*3/5,i));
end
for ii=1:5
[X(ii,:),Y(ii,:)] = func_weighted_above(squeeze(dat_rm(:,ii,:)),squeeze(dat_eps(:,ii,:)),hh,epsmin);
[X(ii+5,:),Y(ii+5,:)] = func_weighted_below(squeeze(dat_rm(:,ii,:)),squeeze(dat_eps(:,ii,:)),hh,epsmin);
end
X=X';
Y=Y';
B = figure('position',[488,358.6,810,303.4]);
func_dash_and_solid(repmat(aer_tick',1,10),Y, 1 ,[10,100000],'$N_{a} (cm^{-3})$','log',[0.1,0.46],C2)
set(gca,'xtick',[10,100,1000,10^4,10^5])
func_dash_and_solid(X,Y, 2 ,[0,30],'$\mathbf{\overline{r} (\mu m)}$','linear',[0.1,0.46],C2)

lgwd = {'0.0~1.0  , Above';'0.0~0.1  , Above';'0.1~0.5  , Above';'0.5~0.85, Above';'0.85~1.0, Above';...
    '0.0~1.0  , Below';'0.0~0.1  , Below';'0.1~0.5  , Below';'0.5~0.85, Below';'0.85~1.0, Below'}; % legend content
legend(lgwd,'Location','eastoutside','fontsize',11,'position',[0.775,0.216,0.199,0.672])
legend('boxoff')
text(1.23,0.947636983613718,'AF Ranges','unit','normal',...
    'fontsize',11,'FontName','Times New Roman')
print('-dpng',B,OutFigName,'-r450')
%%
function [dat_eps,dat_rm] =func_regions_2pro(path1,mnt_tick,ia,nx,nz,lwcad)
nmnt = length(mnt_tick);
mnt_interval = mnt_tick(2)-mnt_tick(1);
[epsmean,rmmean] = deal(nan(nz,nmnt,5));
ii = 0;
for im = mnt_tick
    ii = ii+1;
    ncpath = [cell2mat(path1),'/wrfbin_d01_0001-01-01_0',...
        num2str(floor(im/60),'%01d'),':',num2str(mod(im,60),'%02d'),':00'];
    [qc,eps,rm,~,~] = func_get_pro(ncpath,nx,nz);
    
    %---------------lwcmax---------------%
%     af = qc./repmat(max(max(qc,[],1),[],2),nx,nx);
    %---------------lwcad----------------%
    rho=1/double(ncread(ncpath,'ALT'));% m3/kg
    lwc = qc.*rho;
    af = lwc./repmat(permute(lwcad(:,(im-mnt_tick(1)+mnt_interval)/mnt_interval,ia),[3,2,1]),nx,nx);

    icloud0 = qc>10^-5; 
    epsmean(:,ii,1) = func_HorizonAve(eps,icloud0);
    rmmean(:,ii,1) = func_HorizonAve(rm,icloud0);
    icloud0 = qc>10^-5 & af<=0.1;
    epsmean(:,ii,2) = func_HorizonAve(eps,icloud0);
    rmmean(:,ii,2) = func_HorizonAve(rm,icloud0);
    icloud0 = qc>10^-5 & af>0.1 & af<=0.5;
    epsmean(:,ii,3) = func_HorizonAve(eps,icloud0);
    rmmean(:,ii,3) = func_HorizonAve(rm,icloud0);
    icloud0 = qc>10^-5 & af>0.5 & af<=0.85;
    epsmean(:,ii,4) = func_HorizonAve(eps,icloud0);
    rmmean(:,ii,4) = func_HorizonAve(rm,icloud0);
    icloud0 = qc>10^-5 & af>0.85;
    epsmean(:,ii,5) = func_HorizonAve(eps,icloud0);
    rmmean(:,ii,5) = func_HorizonAve(rm,icloud0);
end
dat_rm = squeeze(nanmean(rmmean,2));
dat_eps = squeeze(nanmean(epsmean,2));
end

function [qc,eps,rm,sig,N] = func_get_pro(ncpath,nx,nz)
global m r varls
nbins = length(varls);
qc = zeros(nx,nx,nz);
for ibin=1:nbins
    q=double(ncread(ncpath,varls(ibin,:)));
    qc = qc+q;
end
icloud00 = qc>10^-5;
[N,rm,sig] = deal(zeros(size(icloud00)));
for ibin=1:nbins
    q=double(ncread(ncpath,varls(ibin,:)));
    q(~icloud00) = nan;
    npdr=q./m(ibin); % 10^6 kg-1
    rm = rm+npdr.*r(ibin);
    N = N+npdr;
end
rm = rm./N;
for ibin=1:nbins
    q=double(ncread(ncpath,varls(ibin,:)));
    q(~icloud00) = nan;
    npdr=q./m(ibin);
    sig = sig+npdr.*(rm-r(ibin)).^2;
end
sig = sqrt(sig./N);
eps = sig./rm;
end

function rsl = func_HorizonAve(var3dim,icloud0)
shadeflag=0; % get no edges of shade
if shadeflag
    [~,~,nz]=size(icloud0);
    rsl = nan(nz,3);
    for iz = 1:nz
        icloud1 = icloud0(:,:,iz);
        if all(~icloud1)
            continue
        end
        data = var3dim(:,:,iz);
        data = data(icloud1);
        rsl(iz,1) = mean(data);
% type 1: quartiles
        rsl(iz,2:3) = quantile(data,[0.25,0.75]);
% type 2: standard deviation
     %   sig = std(data);
     %   rsl(iz,2) = rsl(iz,1)-sig;
     %   rsl(iz,3) = rsl(iz,1)+sig;
    end
else
    var3dim(~icloud0) = 0;
    numicloud = squeeze(sum(sum(icloud0)));
    rsl = squeeze(sum(sum(var3dim)))./numicloud;
end
end

function func_dash_and_solid(X,Y,i_p,xlmt,xlbl,xscl,ylmt,clr_list)
%-------------figure settings------------------%
para_xylbl = {'FontSize',16,'FontWeight','bold'};
para_axis = {'linewidth',1,'FontName','Times New Roman'};
titlewd = 'abcdefghijklmn';
pstn_list=[0.07,0.19,0.32,0.77;
0.45,0.19,0.32,0.77];

subplot('position',pstn_list(i_p,:))
for i=1:5
plot(X(:,i),Y(:,i),'-o','Markersize',5,...
    'LineWidth',1.5,'color',clr_list(i,:));
hold on
end
for i=1:5
plot(X(:,i+5),Y(:,i+5),'--*','Markersize',4,....
    'LineWidth',1.5,'color',clr_list(i,:));
hold on
end
box('on')
xlim(xlmt)
ylim(ylmt)
set(gca,para_axis{:},'xscale',xscl,...
   'XMinorGrid','off','ytick',0:0.1:1,'FontSize',12)
xlabel(xlbl,'interpreter','latex',para_xylbl{:})
ylabel('\epsilon',para_xylbl{:})
grid('on')
title(['(',titlewd(i_p),')'],para_xylbl{:},'unit','normalized','position',[0.075075085009182,0.865695792880259,0])
end

function [xm,ym] =func_weighted_above(xh,yh,H,epsmin)
[~,na]=size(xh);
for i=1:na
    b=epsmin(i);
    xh(1:b,i)=0;
    yh(1:b,i)=0;
    H(1:b,i)=0;
end
H(isnan(xh))=0;
yh(isnan(xh))=0;
xh(isnan(xh))=0;
xm=sum(H.*xh)./sum(H);
ym=sum(H.*yh)./sum(H);
end

function [xm,ym] =func_weighted_below(xh,yh,H,epsmin)
[~,na]=size(xh);
for i=1:na
    b=epsmin(i);
    xh(b+1:end,i)=0;
    yh(b+1:end,i)=0;
    H(b+1:end,i)=0;
end
H(isnan(xh))=0;
yh(isnan(xh))=0;
xh(isnan(xh))=0;
xm=sum(H.*xh)./sum(H);
ym=sum(H.*yh)./sum(H);
end
