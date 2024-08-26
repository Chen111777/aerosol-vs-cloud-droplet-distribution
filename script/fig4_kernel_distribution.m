clear; clc; close all
paths = {'~/WRFV4.5.1/a50';'~/WRFV4.5.1/a100';'~/WRFV4.5.1/a200';'~/WRFV4.5.1/a500';...
'~/WRFV4.5.1/a1000';'~/WRFV4.5.1/a2000';'~/WRFV4.5.1/a5000';...
'~/WRFV4.5.1/a10000';'~/WRFV4.5.1/a20000';'~/WRFV4.5.1/a50000'};
aer_tick= [50,100,200,500,1000,2000,5000,10000,20000,50000];

mnt_interval = 3;
mnt_tick = 75:mnt_interval:135;
nmnt = length(mnt_tick);
InputFileName = ['lwcad_exp1_75to135.mat'];
InterFileName = 'kernel_distribution.mat';
OutFigName = 'kernel_distribution';
len_aer = length(aer_tick);
flp_hot_clmp=flipud(colormap('hot'));
%%
load(InputFileName)
[dat_eps,dat_rbar,dat_z,dat_af] =func_data_point(paths,mnt_tick,lwcad);
length(dat_rbar)
save(InterFileName,'dat_eps','dat_rbar','dat_z','dat_af')
%%
% close all
clc
load(InterFileName)
pstn_list=[0.085,0.55,0.34,0.370;...
    0.500,0.55,0.34,0.370;...
    0.085,0.100,0.34,0.370;...
    0.500,0.100,0.34,0.370];
para_xylbl = {'FontSize',14,'interpreter','latex'};
para_axis = {'LineWidth',1,'FontName','Times New Roman'};
nbwd = 'abcdefghijklmn';
figure('position',[157.8,64.2,710,597.8])
i_p=0;
for ia=[2,6] % Na=100 and 2,000 cm-3
x=dat_rbar(:,ia);
y=dat_eps(:,ia);
z=dat_z(:,ia);
af=dat_af(:,ia);
x(x==0) = [];
y(y==0) = [];
z(z==0) = [];
af(af==0) = [];
deltaf=0.2;
af_edge=0:deltaf:0.8;
deltz=0.2;
hacb_edge=0:deltz:3;
deltx=(ceil(max(x))-floor(min(x)))/25;
x_edge=floor(min(x)):deltx:ceil(max(x));
delty=0.025;
y_edge=floor(min(y)):delty:ceil(max(y));
len_x=length(x_edge)-1;
len_y=length(y_edge)-1;
len_hacb=length(hacb_edge)-1;
len_af=length(af_edge)-1;
[hacbsum_bin,bin_count,afsum_bin]=deal(zeros(len_y,len_x));
for i=1:len_x
    for j=1:len_y
        bin_count(j,i) = sum(x<x_edge(i+1) & x>=x_edge(i) & y<y_edge(j+1) & y>=y_edge(j));
        hacbsum_bin(j,i) = sum(z(x<x_edge(i+1) & x>=x_edge(i) & y<y_edge(j+1) & y>=y_edge(j)));
        afsum_bin(j,i) = sum(af(x<x_edge(i+1) & x>=x_edge(i) & y<y_edge(j+1) & y>=y_edge(j)));
    end
end
[hacb_count,x_hacb,y_hacb]=deal(zeros(1,len_hacb));
for i=1:len_hacb
        hacb_count(i) = sum(z<hacb_edge(i+1) & z>=hacb_edge(i));
        x_hacb(i) = sum(x(z<hacb_edge(i+1) & z>=hacb_edge(i)));
        y_hacb(i) = sum(y(z<hacb_edge(i+1) & z>=hacb_edge(i)));
end
[af_count,x_af,y_af]=deal(zeros(1,len_af));
for i=1:len_af
        af_count(i) = sum(af<af_edge(i+1) & af>=af_edge(i));
        x_af(i) = sum(x(af<af_edge(i+1) & af>=af_edge(i)));
        y_af(i) = sum(y(af<af_edge(i+1) & af>=af_edge(i)));
end
pdf=bin_count/delty/deltx/sum(sum(bin_count));
[X,Y]=meshgrid((x_edge(1:end-1)+x_edge(2:end))/2,(y_edge(1:end-1)+y_edge(2:end))/2);
i_p=i_p+1;
subplot('position',pstn_list(i_p,:))
h=pcolor(X,Y,pdf); 
shading flat
caxis([min(min(pdf)),max(max(pdf))])
set(h,'handlevisibility','off');
hold on
contour(X,Y,hacbsum_bin./bin_count,'k','showtext','on','LevelStep',0.2,'LineWidth',1)
b=corrcoef(x,y);
text(0.644,0.777,{['CC=',num2str(b(1,2),'%.2f')];['     ¦Å=',num2str(mean(y),'%.2f')],},'FontName','Times New Roman',...
    'FontSize',13,'units','normalized')
hold on
scatter(x_hacb./hacb_count,y_hacb./hacb_count,36,'b','filled','MarkerEdgeColor','k');
legend('HFCB (km)','location','NorthEast')
ylim([0.02,0.7])
v=axis;
xlim([3,v(2)])
if i_p<=2
    xlim([7,v(2)])
end
set(gca,para_axis{:})
xlabel('$\overline{r} (\mu$m)',para_xylbl{:});
ylabel('¦Å',para_xylbl{:});
title(['(',nbwd(i_p),')'],'fontsize',16,'FontWeight','bold',...
    'unit','normalized','position',[0.07,0.89,0])

i_p=i_p+1;
subplot('position',pstn_list(i_p,:))
h=pcolor(X,Y,pdf); 
shading flat
caxis([min(min(pdf)),max(max(pdf))])
set(h,'handlevisibility','off');
hold on
contour(X,Y,afsum_bin./bin_count,'k','showtext','on','LevelStep',0.2,'LineWidth',1);
text(0.644,0.777,{['CC=',num2str(b(1,2),'%.2f')];['     ¦Å=',num2str(mean(y),'%.2f')],},'FontName','Times New Roman',...
    'FontSize',13,'units','normalized')
hold on
scatter(x_af./af_count,y_af./af_count,36,'g','filled','MarkerEdgeColor','k');
legend('AF','location','NorthEast')
ylim([0.02,0.7])
v=axis;
xlim([3,v(2)])
if i_p<=2
    xlim([7,v(2)])
end
set(gca,para_axis{:})
xlabel('$\overline{r} (\mu$m)',para_xylbl{:});
ylabel('¦Å',para_xylbl{:});
title(['(',nbwd(i_p),')'],'fontsize',16,'FontWeight','bold',...
    'unit','normalized','position',[0.07,0.89,0])

colormap(flp_hot_clmp)
cb=colorbar('position',[0.872,pstn_list(i_p,2),0.0252,0.369]);
set(get(cb,'ylabel'),'string','PDF (\mum^-^1)')
end
% print('-dpng',gcf,OutFigName,'-r450')

function [dat_eps,dat_rbar,dat_hacb,dat_af] =func_data_point(path,mnt_tick,lwcad)
g=9.81;
nmnt = length(mnt_tick);
mnt_interval = mnt_tick(2)-mnt_tick(1);
[nx,~,nz] = size(double(ncread([cell2mat(path(1)),'/wrfbin_d01_0001-01-01_00:00:00'],'ALT')));
[dat_rbar,dat_eps,dat_hacb,dat_af]=deal(zeros(2500*nmnt,10));
im = 0;
for ia = 1:10
    ia
    i_fill = 1;
    for mnt = mnt_tick
        im = im+1;
        ncfile = [cell2mat(path(ia)),'/wrfbin_d01_0001-01-01_0',...
            num2str(floor(mnt/60),'%01d'),':',num2str(mod(mnt,60),'%02d'),':00'];
        [qc,eps0,rbar0,~,~] = func_get_grid_properties(ncfile,nx,nz);

        ph = double(ncread(ncfile,'PHB'))+double(ncread(ncfile,'PH'));
        z0 = (ph(:,:,2:end)+ph(:,:,1:end-1))/2000/g; % unit: km
        izb = find(any(any(qc>10^-5,1),2),1)-1; % cloud base index
        hacb=z0-repmat(z0(:,:,izb),1,1,nz); % height above cloud base
        rho=1/double(ncread(ncfile,'ALT'));% kg/m3, or 10^-6 kg/cm3
        lwc = qc.*rho;
        af = lwc./repmat(permute(lwcad(:,(mnt-mnt_tick(1)+mnt_interval)/mnt_interval,ia),[3,2,1]),nx,nx);
        
        icloud = qc>10^-5 & hacb>0.2;
        dat_afi=af(icloud);
        dat_epsi = eps0(icloud);
        dat_rbari = rbar0(icloud);
        dat_hacbi = hacb(icloud);
        i_fill_s = i_fill;
        i_fill = i_fill+length(dat_epsi);
        dat_eps(i_fill_s:i_fill-1,ia) = dat_epsi;
        dat_af(i_fill_s:i_fill-1,ia) = dat_afi;
        dat_rbar(i_fill_s:i_fill-1,ia) = dat_rbari;
        dat_hacb(i_fill_s:i_fill-1,ia) = dat_hacbi;
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
