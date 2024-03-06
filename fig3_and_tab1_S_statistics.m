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
InterFileName = 'dat_fig3_10cases_exp1.mat';
OutFigName = 'fig3_exp1';
mnt_interval = 6;
mnt_tick = 78:mnt_interval:132;
nz = 180;%80 
nx = 100;%50;

len_aer = length(aer_tick);
nbins = 33;
varls = []; % variable name of bin mixing ratio
for ibin = 1:nbins
    varls = [varls;['ff1i',num2str(ibin,'%02d')]];
end
a=2.53E12; % coefficients for S
b=5.42E3;
s_bnd=[0:0.01:7]; % boundaries of S bins (%), exclude negative values
n_sbin = length(s_bnd)-1;
%% get S pdf
load(InputFileName)
dat = zeros(len_aer,n_sbin,4);
for ia = 1:len_aer
    ia
    n_exist=zeros(1,4);
    for im = mnt_tick
    ncpath = [cell2mat(paths(ia)),'/wrfbin_d01_0001-01-01_0',...
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
    af = lwc./repmat(permute(lwcad(:,(im-mnt_tick(1)+mnt_interval)/mnt_interval,ia),[3,2,1]),nx,nx);
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
    n_exist(iaf)=n_exist(iaf)+1;
    for ix = 1:n_sbin
        n=length(ss(ss>s_bnd(ix) & ss<s_bnd(ix+1)));
        dat(ia,ix,iaf) = dat(ia,ix,iaf)+n/length(ss); 
        % distribute S into bins; SS in denominator includes negative values
    end
    end
end
    
    end
for iaf=1:4
    dat(ia,:,iaf) = dat(ia,:,iaf)/n_exist(iaf); 
    % ratio of S at a specific bin to all S, averaged over time
end 
end
save(InterFileName,'dat')
%% statistics for s_bar and epss_s
load(InterFileName)
s_tick=(s_bnd(2:end)+s_bnd(1:end-1))/2;
[sbar,sbar_trun,epss,epss_trun,skew_trun] = deal(ones(len_aer,4));
for ia = 1:len_aer
    for iaf=1:4
    pdf=dat(ia,:,iaf);
    rt(ia,iaf) = 100-sum(pdf)*100;
    id = 0;
    sbar(ia,iaf)=sum(s_tick.*pdf)/sum(pdf);
    sigs=sqrt(sum((s_tick-sbar(ia,iaf)).^2.*pdf)/sum(pdf));
    epss(ia,iaf)=sigs/sbar(ia,iaf);
    while skew_trun(ia,iaf)>0.05
        pdf(end-id) = 0; % truncate from right to left
        sbar1=sum(s_tick.*pdf)/sum(pdf);
        sigs=sqrt(sum((s_tick-sbar1).^2.*pdf)/sum(pdf));
        skew_trun(ia,iaf)=sum(((s_tick-sbar1)/sigs).^3.*pdf)/sum(pdf);
        id = id+1;
    end
    epss_trun(ia,iaf)=sigs/sbar1;
    sbar_trun(ia,iaf)=sbar1;
end
end
num2str(rt','%7.3f') % display negative S ratio - Table 1
% num2str(epss','%7.3f')
% num2str(epss_trun','%7.3f')
% num2str(sbar_trun','%7.3f')
B = figure('position',[488,358.6,810,303.4]);
C2=[50,50,50;78.4,124.8,204;140.8,186.4,183.2;224.4,168.96,178.64;204.5,82.6,77.0]/255;
X=repmat(aer_tick',1,4);
func_dash_and_solid(X,[epss(:,3:4),epss_trun(:,3:4)], 1 ,'$N_{a} (cm^{-3})$','$\epsilon_{S} (S>0)$',C2(end-1:end,:))
func_dash_and_solid(X,[sbar(:,3:4),sbar_trun(:,3:4)], 2 ,'$N_{a} (cm^{-3})$','$\mathbf{\overline{S}} (\%, S>0)$',C2(end-1:end,:))
lgwd = {'0.5~0.85, Original';'0.5~0.85, Zero-skewed';'0.85~1.0, Original';'0.85~1.0, Zero-skewed'};
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
