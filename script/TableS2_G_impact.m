clear; clc; close all
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
len_aer = length(paths);
mnt_interval = 3;
mnt_tick = 75:mnt_interval:135;
idx = 0;
[std1,std2,std3,num_clouds,largest_cloud_ratio] = deal(zeros(len_aer,length(mnt_tick)));
for ia = 1:len_aer
    ia
    im = 0;
    for mnt = mnt_tick
        im = im+1;
        ncfile = [cell2mat(paths(ia)),expname,'/wrfbin_d01_0001-01-01_0',...
            num2str(floor(mnt/60),'%01d'),':',num2str(mod(mnt,60),'%02d'),':00'];
        [std1(ia,im),std2(ia,im),std3(ia,im)] = func_std_3types(ncfile);
        [largest_cloud_ratio(ia,im), ~] = func_cloud_counter(ncfile);
    end
end
nanmean(std1(:))
nanmean(std2(:))
nanmean(std3(:))
nanmean(largest_cloud_ratio,2)

function [std1,std2,std3] = func_std_3types(ncfile)
P = double(ncread(ncfile,'P'))+double(ncread(ncfile,'PB'));
TH = double(ncread(ncfile,'T'))+300;
T = TH.*(P./10^5).^0.286;
qc = double(ncread(ncfile,'QCLOUD'));
TZERO = 273.150;
PZERO = 1.013E6;
D0 = 0.221;
Dv = D0.*(PZERO./P).*(T./TZERO).^1.94;
CF_MY = 2.4E3;
Lv = 2.5E10;
Rv = 461.5E4;
a=2.53E12;
b=5.42E3;
es=a/10*exp(-b./T);
qv = double(ncread(ncfile,'QVAPOR'));
ew=qv.*P./(0.622+0.378*qv);
s = ew./es*100-100;

ind = s>0;
% ind = qc>10^-5;

s = s(ind);
FD = Rv*T./Dv./es;
FK = (Lv/(Rv*T)-1.)*Lv/CF_MY./T;
G = FD+FK;
G = G(ind);
std1 = std(s./G);
std2 = std(s./mean(G));
std3 = std(mean(s)./G);
end


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