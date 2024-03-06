clear; clc; close all
% --------- change these ------------%
paths = {'~/WRFV4.5.1/a50';'~/WRFV4.5.1/a2000';'~/WRFV4.5.1/a50000'};
lgwd = {'N_a=50';'N_a=2,000';'N_a=50,000'}; % legend content
InputFileName = 'lwcad_3cases_exp1_75to135.mat';
InterFileName = 'dat_fig1_3cases_exp1.mat';
OutFigName = 'fig1_exp1';
mnt_interval = 3;
mnt_tick = 75:mnt_interval:135;

len_aer = length(paths);
global m r varls zz
r = 2.^(1:1/3:35/3); % bin radius
m = 4/3*pi.*r.*r.*r/10^9; % bin mass
varls = []; % variable name of bin mixing ratio
for ibin = 1:33
    varls = [varls;['ff1i',num2str(ibin,'%02d')]];
end
%%
g=9.81;
fstpath = [cell2mat(paths(1)),'/wrfbin_d01_0001-01-01_01:54:00'];
phb = double(ncread(fstpath,'PHB'));
php = double(ncread(fstpath,'PH'));
ph = phb+php; % because ph merely varies with time and aerosol
[nx,~,nz]=size(ph);
nz=nz-1;
ph1d=squeeze(mean(mean(ph,1),2))/g/1000;
zz = (ph1d(1:end-1)+ph1d(2:end))/2;
clear phb php ph

% get three-region r_mean, epsilon and sigma profiles
load(InputFileName)
[dat1,dat2,dat3,dat4,dat5,dat6,...
    dat7,dat8,dat9] =deal(nan(nz,3,len_aer));
for ia = 1:len_aer
    ia
    [dat1(:,:,ia),dat2(:,:,ia),dat3(:,:,ia),...
        dat4(:,:,ia),dat5(:,:,ia),dat6(:,:,ia),...
        dat7(:,:,ia),dat8(:,:,ia),dat9(:,:,ia)] = func_regions_3pro(paths(ia),mnt_tick,ia,nx,nz,lwcad);
end
save(InterFileName,'dat1','dat2','dat3','dat4','dat5','dat6','dat7','dat8','dat9','zz')
%%
load(InterFileName)
%-------------figure settings------------------%
para_xylbl = {'FontSize',16,'FontWeight','bold','FontName','Times New Roman'};
para_axis = {'linewidth',1,'FontName','Times New Roman',...
    'FontSize',10};
pst = [0.1100    0.7060    0.252    0.255;
   0.4068    0.7060    0.252    0.255;
   0.7036    0.7060    0.252    0.255;
   0.1100    0.3900    0.252    0.255;
   0.4068    0.3900    0.252    0.255;
   0.7036    0.3900    0.252    0.255;
   0.1100    0.0740    0.252    0.255;
   0.4068    0.0740    0.252    0.255;
   0.7036    0.0740    0.252    0.255];
xlb = {'$\overline{r} (\mu m)$';'$\epsilon$';'$\sigma (\mu m)$'};
ttlwd = {'$\mathbf{AF>0.85}$';'$\mathbf{0.5<AF\leq 0.85}$';'$\mathbf{AF\leq 0.5}$'};
nbwd = {'(a1)';'(a2)';'(a3)';'(b1)';'(b2)';'(b3)';'(c1)';'(c2)';'(c3)'};

B = figure('position',[288,70,850,680]);
shadeflag = 1; % draw shade
for ip = 1:9
   ax = subplot('position',pst(ip,:));
   for ia = 1:len_aer
       eval(['dat = dat',num2str(ip),';'])
       hold on;
       h = plot(dat(:,1,ia),zz);
       set(h,'Linewidth',1)
       if shadeflag
           clr = get(h,'color');
           hold on;
           x=[dat(:,3,ia);flipud(dat(:,2,ia))];
           y=[zz;flipud(zz)];
           sh = fill(x(~isnan(x)),y(~isnan(x)),...
               'm','FaceColor',clr,'FaceAlpha',0.15,...
               'EdgeColor','none','handlevisibility','off');
       end
   end
   set(ax,para_axis{:})
   if mod(ip,3)==1
       ylabel('Altitude (km)','FontSize',16)
       text(-0.36,1.07058824931874,ttlwd(ceil(ip/3)),'interpreter','latex','unit','normal',...
           'FontSize',15)
   end
   if ip>=7
       xlabel(xlb(ip-6),'interpreter','latex',para_xylbl{:})
   end
   if ip==1
       legend(lgwd,'Location','northeast','FontSize',11)
       legend(ax,'boxoff')
   end
   box('on')
   grid('on')
   ylim([0.5,3.8])
   title(nbwd(ip),'fontsize',14,...
       'units','normalized','position',[0.1,0.85,0],'FontWeight','bold');
end
print('-dpng',B,OutFigName,'-r450')

function [dat1,dat2,dat3,dat4,dat5,dat6,dat7,dat8,dat9] =func_regions_3pro(path1,mnt_tick,ia,nx,nz,lwcad)
% datn is data for subplot n. 
nmnt = length(mnt_tick);
mnt_interval = mnt_tick(2)-mnt_tick(1);
[epsmean0,smean0,sigmean0,...
    epsmean4,smean4,sigmean4,...
    epsmean8,smean8,sigmean8,] = deal(nan(nz,3,nmnt)); 
% column 1: mean; column 2: lower quartile; column 3: upper quatile
ii = 0;
for im = mnt_tick
    ii = ii+1;
    pathbin = [cell2mat(path1),'/wrfbin_d01_0001-01-01_0',...
        num2str(floor(im/60),'%01d'),':',num2str(mod(im,60),'%02d'),':00'];
    [qc,eps,rm,sig] = func_get_pro(pathbin,nx,nz);
    %---------------lwcmax---------------%
%     af = qc./repmat(max(max(qc,[],1),[],2),nx,nx);
    %---------------lwcad----------------%
    rho=1/double(ncread(pathbin,'ALT'));% m3/kg
    lwc = qc.*rho;
    af = lwc./repmat(permute(lwcad(:,(im-mnt_tick(1)+mnt_interval)/mnt_interval,ia),[3,2,1]),nx,nx);
    icloud0 = qc>10^-5 & af<=0.5;
    epsmean0(:,:,ii) = func_HorizonAve(eps,icloud0);
    sigmean0(:,:,ii) = func_HorizonAve(sig,icloud0);
    smean0(:,:,ii) = func_HorizonAve(rm,icloud0);
    icloud0 = qc>10^-5 & af>0.5 & af<=0.85;
    epsmean4(:,:,ii) = func_HorizonAve(eps,icloud0);
    sigmean4(:,:,ii) = func_HorizonAve(sig,icloud0);
    smean4(:,:,ii) = func_HorizonAve(rm,icloud0);
    icloud0 = qc>10^-5 & af>0.85;
    epsmean8(:,:,ii) = func_HorizonAve(eps,icloud0);
    sigmean8(:,:,ii) = func_HorizonAve(sig,icloud0);
    smean8(:,:,ii) = func_HorizonAve(rm,icloud0);
end

dat7 = nanmean(smean0,3);
dat8 = nanmean(epsmean0,3);
dat9 = nanmean(sigmean0,3);
dat4 = nanmean(smean4,3);
dat5 = nanmean(epsmean4,3);
dat6 = nanmean(sigmean4,3);
dat1 = nanmean(smean8,3);
dat2 = nanmean(epsmean8,3);
dat3 = nanmean(sigmean8,3);
end

function [qc,eps,rm,sig] = func_get_pro(pathbin,nx,nz)
global m r varls
nbins = length(varls);
qc = zeros(nx,nx,nz);
for ibin=1:nbins
    q=double(ncread(pathbin,varls(ibin,:)));
    qc = qc+q;
end
icloud00 = qc>10^-5;
[N,rm,sig] = deal(zeros(size(icloud00)));
for ibin=1:nbins
    q=double(ncread(pathbin,varls(ibin,:)));
    q(~icloud00) = nan;
    npdr=q./m(ibin); % 10^6 kg-1
    rm = rm+npdr.*r(ibin);
    N = N+npdr;
end
rm = rm./N;
for ibin=1:nbins
    q=double(ncread(pathbin,varls(ibin,:)));
    q(~icloud00) = nan;
    npdr=q./m(ibin);
    sig = sig+npdr.*(rm-r(ibin)).^2;
end
sig = sqrt(sig./N);
eps = sig./rm;
end

function rsl = func_HorizonAve(var3dim,icloud0)
shadeflag = 1; % get edges of shade
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
    %    sig = std(data);
    %    rsl(iz,2) = rsl(iz,1)-sig;
    %    rsl(iz,3) = rsl(iz,1)+sig;
    end
else
    var3dim(~icloud0) = 0;
    numicloud = squeeze(sum(sum(icloud0)));
    rsl = squeeze(sum(sum(var3dim)))./numicloud;
end
end
