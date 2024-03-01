clear; clc; close all
% --------- change these ------------%
paths = {'~/WRFV4.5.1/a50_exp3',...
'~/WRFV4.5.1/a2000_exp3',...
'~/WRFV4.5.1/a10000_exp3','~/WRFV4.5.1/a50000_exp3'};
aer_tick=[50,2000,10000,50000];%[50,100,200,500,1000,2000,5000,10000,20000,50000];
InputFileName = 'lwcad_4cases_exp3.mat';
InterFileName = 'dat_fig3_4cases_exp2.mat';
OutFigName = 'fig3_exp2';
nz = 180; nx = 100;

len_aer = length(aer_tick);
nbins = 33;
varls = []; % variable name of bin mixing ratio
for ibin = 1:nbins
    varls = [varls;['ff1i',num2str(ibin,'%02d')]];
end
a=2.53E12; % coefficients for S
b=5.42E3;
mnt_range = 72:6:132;
nmnt = length(mnt_range);
s_bnd=[-21:0.01:7]; % boundaries of S bins, unit: %
n_sbin = length(s_bnd)-1;
%% get S pdf
ii = 0;
for accn=aer_tick
    ii = ii+1
    nmnt=zeros(1,4);
    for im = mnt_range
    ncpath = [cell2mat(paths(ii)),'/wrfbin_d01_0001-01-01_0',...
          num2str(floor(im/60),'%01d'),':',num2str(mod(im,60),'%02d'),':00'];
    qc = zeros(nx,nx,nz);
    for ibin=1:nbins
        q=double(ncread(ncpath,varls(ibin,:)));
        qc = qc+q;
    end
    qv = double(ncread(ncpath,'QVAPOR'));
    pp = double(ncread(ncpath,'P'));
    pb = double(ncread(ncpath,'PB'));
    P = pp+pb;
    TH=double(ncread(ncpath,'T'))+300;
    T=TH.*(P./10^5).^0.286;
    es=a/10*exp(-b/T);
    ew=qv.*P./(0.622+0.378*qv);
    s0 = ew./es*100-100;
    %---------------lwcmax---------------%
%     af = qc./repmat(max(max(qc,[],1),[],2),nx,nx);
    %---------------lwcad----------------%
    rho=1/double(ncread(ncpath,'ALT'));% m3/kg
    lwc = qc.*rho;
    af = lwc./repmat(permute(lwcad(:,(im-66)/6,ii),[3,2,1]),nx,nx);
for iaf=1:4
    if iaf==1 
    ss = s0(qc>10^-5 & af<=0.1);
    elseif iaf==2
    ss = s0(qc>10^-5 & af>0.1 & af<=0.5);
    elseif iaf==3
    ss = s0(qc>10^-5 & af>0.5 & af<=0.85);
    elseif iaf==4
    ss = s0(qc>10^-5 & af>0.85);
    end
    if isempty(ss)
        continue % no this region, no S, so skip the next step.
    else
    nmnt(iaf)=nmnt(iaf)+1;
    for ix = 1:n_sbin
        n=length(ss(ss>s_bnd(ix) & ss<s_bnd(ix+1)));
        dat(ii,ix,iaf) = dat(ii,ix,iaf)+n/length(ss); % distribute S into bins
    end
    end
end
    
    end
for iaf=1:4
    dat(ii,:,iaf) = dat(ii,:,iaf)/nmnt(iaf); % ratio of S at a specific bin to all S, averaged over time
end 
end
save(InterFileName,'dat')
%% statistics for s_bar and epss_s
load(InterFileName)
s_tick=(s_bnd(2:end)+s_bnd(1:end-1))/2;
[sbar,sbar_trun,epss,epss_trun,skew_trun] = deal(ones(len_aer,4));
for ip = 1:len_aer
    ip
    for iaf=1:4
    pdf=dat(ip,:,iaf);
    rt(ip,iaf) = sum(pdf(s_tick<0))/sum(pdf)*100;
    id = 0;
    pdf(s_tick<0) = []; % exclude negative value
    s_tick(s_tick<0) = [];
    s_tick(pdf==0) = [];
    pdf(pdf==0) = [];
    sbar(ip,iaf)=sum(s_tick.*pdf)/sum(pdf);
    std=sqrt(sum((s_tick-sbar(ip,iaf)).^2.*pdf)/sum(pdf));
    epss(ip,iaf)=std/sbar(ip,iaf);
    while skew_trun(ip,iaf)>0.05
        pdf(end-id) = 0; % truncate from right to left
        sbar1=sum(s_tick.*pdf)/sum(pdf);
        std=sqrt(sum((s_tick-sbar1).^2.*pdf)/sum(pdf));
        skew_trun(ip,iaf)=sum(((s_tick-sbar1)/std).^3.*pdf)/sum(pdf);
        id = id+1;
    end
    epss_trun(ip,iaf)=std/sbar1;
    sbar_trun(ip,iaf)=sbar1;
end
end
num2str(rt','%7.3f') % display negative S ratio - Table 1
% num2str(epss','%7.3f')
% num2str(epss_trun','%7.3f')
% num2str(sbar_trun','%7.3f')
B = figure('position',[488,358.6,810,303.4]);
C2=[50,50,50;78.4,124.8,204;140.8,186.4,183.2;224.4,168.96,178.64;204.5,82.6,77.0]/255;
X=repmat(aer_tick',1,4);
func_dash_and_solid(X,[epss(:,3:4),epss_trun(:,3:4)], 1 ,'$N_{a} (cm^{-3})$','$\epssilon_{S} (S>0)$',C2(end-1:end,:))
func_dash_and_solid(X,[sbar(:,3:4),sbar_trun(:,3:4)], 2 ,'$N_{a} (cm^{-3})$','$\mathbf{\overline{S}} (\%, S>0)$',C2(end-1:end,:))
lgwd = {'0.5~0.85, Original','0.5~0.85, Zero-skewed','0.85~1.0, Original','0.85~1.0, Zero-skewed'};
legend(lgwd,'Location','eastoutside','fontsize',11,'position',[0.783,0.396,0.200,0.289])
set(gca,'position',[0.455,0.19,0.31,0.77])
legend('boxoff')
text(1.203885350318471,0.686438356164383,'AF Ranges','unit','normal',...
    'fontsize',11,'FontName','Times New Roman')
print('-dpng',gcf,OutFigName,'-r450')

function func_dash_and_solid(X,Y,i_p,xlbl,ylbl,clr_list)
%-------------figure settings------------------%
para_xylbl = {'FontSize',16,'FontWeight','bold'};
para_axis = {'linewidth',1,'FontName','Times New Roman'};
titlewd = 'abcdefghijklmn';
pstn_list=[0.07,0.19,0.31,0.77;
0.455,0.19,0.31,0.77];
subplot('position',pstn_list(i_p,:))
for i=1:2
plot(X(:,i),Y(:,i),'-o','Markersize',5,...
    'LineWidth',2,'color',clr_list(i,:));
hold on
plot(X(:,i+2),Y(:,i+2),'--*','Markersize',5,....
    'LineWidth',2,'color',clr_list(i,:));
hold on
end
box('on')
xlim([10,100000])
set(gca,para_axis{:},'xscale','log',...
   'XMinorGrid','off','FontSize',12,'xtick',[10,100,1000,10000,100000])
xlabel(xlbl,'interpreter','latex',para_xylbl{:})
ylabel(ylbl,'interpreter','latex',para_xylbl{:})
grid('on')
title(['(',titlewd(i_p),')'],para_xylbl{:},'unit','normalized','position',[0.075075085009182,0.865695792880259,0])
end