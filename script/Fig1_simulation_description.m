clear; clc; 
paths = {'~/submit_v1_wrfrsl/a50_';...
'~/submit_v1_wrfrsl/a100_';...
'~/submit_v1_wrfrsl/a200_';...
'~/submit_v1_wrfrsl/a500_';...
'~/submit_v1_wrfrsl/a1000_';...
'~/submit_v1_wrfrsl/a2000_';...
'~/submit_v1_wrfrsl/a5000_';...
'~/submit_v1_wrfrsl/a10000_';...
'~/submit_v1_wrfrsl/a20000_';...
'~/submit_v1_wrfrsl/a50000_'};
expname = 'tke100';
lwcad_file = ['lwcad_',expname,'.mat'];
InterFileName = ['SimulationDescription_',expname,'.mat'];
mnt_interval = 1;
mnt_tick = 0:mnt_interval:135;
nmnt = length(mnt_tick);
len_aer = length(paths);
nx=50; nz=80;
% clmp = flipud(colormap('hot'));
%%
[largest_volume_ratio,cloud_coverage] = deal(nan(len_aer,nmnt));
nbins = 33;
for ia = 1:len_aer
    ia
im = 0;
for mnt = mnt_tick
    im = im+1;
    ncfile = [cell2mat(paths(ia)),expname,'/wrfbin_d01_0001-01-01_0',...
        num2str(floor(mnt/60),'%01d'),':',num2str(mod(mnt,60),'%02d'),':00'];
    qc=double(ncread(ncfile,'QCLOUD'));
    icloud = qc>10^-5;
    cloud_coverage(ia,im)=sum(sum(sum(icloud,3)>0))*0.1*0.1;
    try
            [largest_volume_ratio(ia,im), ~] = func_cloud_counter(ncfile);
    catch
        largest_volume_ratio(ia,im) = nan;
    end
end
end

ia=5;
ncfile = [cell2mat(paths(ia)),expname,'/wrfbin_d01_0001-01-01_01:57:00'];
U=double(ncread(ncfile,'U'));
V=double(ncread(ncfile,'V'));
W=double(ncread(ncfile,'W'));
U=(U(2:end,:,:)+U(1:end-1,:,:))/2;
V=(V(:,2:end,:)+V(:,1:end-1,:))/2;
W=(W(:,:,2:end)+W(:,:,1:end-1))/2;
u_2d=squeeze(U(:,nx/2,:))';
w_2d=squeeze(W(:,nx/2,:))';
X=repmat(0.05:0.1:4.95,nz,1);
ph = double(ncread(ncfile,'PHB'))+double(ncread(ncfile,'PH'));
ph=(ph(:,:,2:end)+ph(:,:,1:end-1))/2;
zz = squeeze(mean(mean(ph/9.81/1000))); % unit: km
Z=repmat(zz,1,nx);
load(lwcad_file)
rho=1/double(ncread(ncfile,'ALT'));% kg/m3, or 10^-6 kg/cm3
qc=double(ncread(ncfile,'QCLOUD'));
qc(qc<10^-5)=nan;
lwc = qc.*rho;
im = 120;
idx_in_lwcad = (im-mnt_tick(1)+mnt_interval)/mnt_interval;
af = lwc./repmat(permute(lwcad(:,idx_in_lwcad),[3,2,1]),nx,nx);
af_2d = squeeze(af(:,nx/2,:))';
save(InterFileName,...
    'mnt_tick','cloud_coverage','largest_volume_ratio',...
   'X','Z','af_2d','u_2d','w_2d')
%%
load(InterFileName)
figure('position',[488,342,696.2,440])
pst_list = [0.08,0.56,0.215,0.32;...
    0.40,0.56,0.20,0.32;...
    0.70,0.56,0.215,0.32;...
    0.08,0.10,0.215,0.32;...
    0.40,0.10,0.215,0.32;...
    0.70,0.10,0.20,0.32];
dm_area=nx*nx*0.1*0.1;
cloud_fraction_time=mean(cloud_coverage)/dm_area;
largest_volume_ratio_time=mean(largest_volume_ratio);

i_p = 4;
subplot(2,3,i_p)
plot(mnt_tick,cloud_fraction_time,'b','LineWidth',1.5)
v=axis;
hold on; plot([75,75],v(3:4),'--k','LineWidth',1)
axis(v)
xlabel('Time (min)')
ylabel('Cloud Fraction')
xlim([0,135])
ylim([0,0.15])

yyaxis right
ylim([0,0.15] .* dm_area)
ylabel('Cloud Coverage (km^2)')

set(gca,'LineWidth',1,'FontName','Times New Roman','XTick',0:20:140,...
    'position',pst_list(i_p,:),'YColor','k')
title('(c)','fontsize',15,'unit','normalized','position',[-0.112,1.042,0])

i_p = 5;
subplot(2,3,i_p)
plot(mnt_tick,largest_volume_ratio_time,'b','LineWidth',1.5)
v=axis;
hold on; plot([75,75],v(3:4),'--k','LineWidth',1)
% hold on; plot(minmax(mnt_tick),[0.94,0.94],'LineWidth',1)
xlabel('Time (min)')
ylabel('Largest Volume Ratio')
xlim([0,135])
% ylim([0.9,1])
set(gca,'LineWidth',1,'FontName','Times New Roman','XTick',0:20:140,...
    'position',pst_list(i_p,:))
title('(d)','fontsize',15,'unit','normalized','position',[-0.112,1.042,0])

i_p=6;
subplot(2,3,i_p)
pcolor(X,Z,af_2d); shading flat
on=zeros(size(X)); on(1:3:end,1:2:end)=1; on=on==1;
hold on; quiver(X(on),Z(on),u_2d(on),w_2d(on),'k','LineWidth',0.3)
%axis([0,5,0,3])
xlabel('x (km)')
ylabel('Altitude (km)')
cb = colorbar('units','normalized');
cb.Label.String = 'AF';
set(gca,'LineWidth',1,'FontName','Times New Roman',...
    'position',pst_list(i_p,:))
% colormap(clmp); 
title('(e)','fontsize',15,'unit','normalized','position',[-0.112,1.042,0])


%%%
i_p = 1;
subplot(2,3,i_p)
wrfin = 'wrfinput_d01';
ph = double(ncread(wrfin,'PHB'))+double(ncread(wrfin,'PH'));
ph=(ph(:,:,2:end)+ph(:,:,1:end-1))/2;
height = squeeze(mean(mean(ph(:,:,:)))/9.81/1000); % unit: km
P = double(ncread(wrfin,'P'))+double(ncread(wrfin,'PB'));
TH = double(ncread(wrfin,'T'))+300;
T=TH.*(P./10^5).^0.286-273.15;
T = squeeze(T(1,1,:));
qv = double(ncread(wrfin,'QVAPOR')); % (g/kg)
qv = squeeze(qv(1,1,:));
ax1 = gca;
plot(T, height, 'r', 'LineWidth', 1.5);
xlabel(ax1, 'Temperature (¡ãC)');
ylabel(ax1, 'Altitude (km)');   
ax1.XColor = 'r';  
ax2 = axes('Position', ax1.Position, 'XAxisLocation', 'top', ...
           'YAxisLocation', 'left', 'Color', 'none');
hold on; plot(ax2, qv, height, 'b--', 'LineWidth', 2);
xlabel(ax2, 'Water Mixing Ratio (g/kg)');
ax2.XColor = 'b';
legend(ax1, 'Temperature','box','off');
legend(ax2, 'Mixing Ratio','box','off');
set([ax1,ax2],'LineWidth',1,'FontName','Times New Roman',...
    'position',pst_list(i_p,:),'YLim',minmax(height'))
title('(a)','fontsize',15,'unit','normalized','position',[-0.112,1.042,0])
ax2.YTick = ax1.YTick; 
ax2.XTick = 0.003:0.003:0.03;

i_p = 2;
subplot(2,3,i_p)
wrfout = 'wrfbin_d01_0001-01-01_00_01_00';
info = ncinfo(wrfout);
r_l = 2.^(1:1/3:35/3); % bin radius
m_l = 4/3*pi.*r_l.*r_l.*r_l/10^9; % bin mass
m_ccn = m_l(2).*2.0.^[-42:0];
r_ccn = (3.*m_ccn./4./pi.*10^9).^(1/3);
Z = zeros(length(height),43);
rho=1/double(ncread(wrfout,'ALT'));
for ibin=1:43
    qccn_bin=double(ncread(wrfout,['ff8i',num2str(ibin,'%02d')]));
    Z(:,ibin) = squeeze(mean(mean(qccn_bin.*rho)))/10^6*log(2)/3; % cm-3
end
delt_lnr = log(2)/3;
plot(r_ccn,Z(1,:)./delt_lnr,'b','LineWidth',1.5)
xlabel('r (\mum)')
ylabel('dN_a/dlnr (cm^-^3)')
grid('on')
set(gca,'LineWidth',1,'FontName','Times New Roman',...
    'position',pst_list(i_p,:),'xscale','log',...
    'xtick',[0.0001,0.001,0.01,0.1,1],'xminorgrid','off')
title('(b)','fontsize',15,'unit','normalized','position',[-0.112,1.042,0])
% sum(Z(1,:))
% qnccn = double(ncread(wrfout,'QNCCN'));
% qnccn(1,1,1)*rho(1,1,1)/10^6

% i_p = 3;
% subplot(2,3,i_p)
% plot(sum(Z,2),height,'b','LineWidth',1.5)
% xlabel('N_a (cm^-^3)')
% ylabel('Altitude (km)')
% grid('on')
% set(gca,'LineWidth',1,'FontName','Times New Roman',...
%     'position',pst_list(i_p,:),'xscale','log','xminorgrid','off')
% title('(c)','fontsize',15,'unit','normalized','position',[-0.112,1.042,0])


% print('-dpng',gcf,'SimulationDescription','-r450')

function [largest_cloud_ratio,num_clouds] = func_cloud_counter(ncfile)
qc = double(ncread(ncfile, 'QCLOUD'));
threshold = 1e-5;  % Threshold for cloud water content
cloud_mask = qc > threshold;  % Create a mask for cloud presence

% Use bwconncomp function to identify 3D connected components
CC = bwconncomp(cloud_mask,6);  % 6 indicates the connectivity in 3D space

% Count the total number of clouds
num_clouds = CC.NumObjects;

% Calculate the volume of each cloud (number of voxels in each connected component)
cloud_volumes = cellfun(@numel, CC.PixelIdxList);
% Calculate the total cloud volume (sum of all cloud volumes)
total_cloud_volume = sum(cloud_volumes);
% Find the largest cloud volumes
sorted_volumes = sort(cloud_volumes, 'descend');  % Sort volumes in descending order
largest_cloud_volume = sorted_volumes(1);
% Calculate the proportions of total cloud volume for the largest and second largest clouds
largest_cloud_ratio = largest_cloud_volume / total_cloud_volume;
end