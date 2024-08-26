clear; clc
% --------- change these ------------%
paths = {'~/WRFV4.5.1/a50_exp3',...
'~/WRFV4.5.1/a2000_exp3',...
'~/WRFV4.5.1/a10000_exp3',...
'~/WRFV4.5.1/a50000_exp3'};
OutFileName = 'lwcad_exp3_a50_2000_10000_50000';

n_na = length(paths);
g=9.81;
global m r a b varls zz
a=2.53E12;
b=5.42E3;
r = 2.^(1:1/3:35/3); % bin radius
m = 4/3*pi.*r.*r.*r/10^9; % bin mass
varls = []; % variable name of bin mixing ratio
for ibin = 1:33
    varls = [varls;['ff1i',num2str(ibin,'%02d')]];
end
fstpath = [cell2mat(paths(1)),'/wrfbin_d01_0001-01-01_01:54:00'];
phb = double(ncread(fstpath,'PHB'));
php = double(ncread(fstpath,'PH'));
ph = phb+php; % because ph merely varies with time
[nx,~,nz]=size(ph);
nz=nz-1;
zz = squeeze(mean(mean((ph(:,:,1:end-1)+ph(:,:,2:end))/2/g/1000,1),2));
clear phb php ph
mnt_range = 72:6:132; % time range
nbins = 33;
%% save lwcad 
 nmnt = length(mnt_range);
 lwcad = nan(nz,nmnt,n_na);
 for ia=1:n_na
     ia
 ii = 0;
 for im = mnt_range
     ii = ii+1;
     pathbin = [cell2mat(paths(ia)),'/wrfbin_d01_0001-01-01_0',...
       num2str(floor(im/60),'%01d'),':',num2str(mod(im,60),'%02d'),':00'];
     pp = double(ncread(pathbin,'P'));
     pb = double(ncread(pathbin,'PB'));
     P = pp+pb;
     TH=double(ncread(pathbin,'T'))+300;
     T=TH.*(P./10^5).^0.286;
     qc = zeros(nx,nx,nz);
     for ibin=1:nbins
         q=double(ncread(pathbin,varls(ibin,:)));
         qc = qc+q;
     end
     rho = 1./double(ncread(pathbin,'ALT'));
     lwc = rho.*qc;
     lwc(qc<10^-5) = 0;
     lwcad(:,ii,ia) = func_inputs_for_LWCad(lwc,P,T);
 end
 end
 save(OutFileName,'lwcad')
%%
function lwcad = func_inputs_for_LWCad(lwc,P,T)
global zz
lwc(isnan(lwc)) = 0;
lwcmax = max(max(lwc));
izb = find(any(any(lwc,1),2),1)-1; % cloud base index
izt = find(any(any(lwc,1),2),1,'last'); % cloud top index
Zb = zz(izb)*1000; % unit: m
lwcb = lwc(:,:,izb);
lwcmax = lwcmax(izb);
T = T(:,:,izb);
P = P(:,:,izb);
Tb = T(lwcb == lwcmax); Tb = Tb(1);
Pb = P(lwcb == lwcmax)/100; Pb = Pb(1); % unit: hPa
lwcad = nan(length(zz),1);
[lwc_result,~,~,~,~] = func_adiabatic_lwc(Zb,Tb,Pb,zz(izb:izt)*1000);
lwcad(izb:izt) = lwc_result/1000; % unit: kg m-3
end

function [ lwc_result,qv_result,ql_result,T_result ,delta_qv] = func_adiabatic_lwc(Zb,Tb,Pb,Zz )
%   input:   Zb: cloud base altitude (m);  Tb: cloud base temperature (K);  
%               Pb: cloud base pressure (hPa);  Zz: aimed altitudes (m);
%   output: qv,ql: g/kg;  LWC: g/m^-3; the same matrix size as Zz. 

nz = length(Zz);
g=9.8; %m s-2
cp=1005; %J K-1 kg-1
Rd=287;%J K-1 kg-1
Rv=461 ;%J K-1 kg-1
T1=Tb;
Z1=Zb;
P1=Pb;
% get cloud base qv
LL=[2.603E+03,2.575E+03,2.549E+03,2.525E+03,2.501E+03,2.489E+03,2.477E+03,2.466E+03,2.453E+03,2.442E+03,2.430E+03,2.418E+03,2.406E+03];
TT=[233,243,253,263,273,278,283,288,293,298,303,308,313];
for i=1:12
    if (Tb>TT(i)&&Tb<=TT(i+1))
        Latent=LL(i)+(LL(i+1)-LL(i))/(TT(i+1)-TT(i))*(Tb-TT(i));
        Latent=Latent*1000;
    end
end

esb=6.11*10^(3.45e-6*Latent*((Tb-273.15)/Tb));
qvb=0.622*esb/(Pb-esb) ;
qv1=qvb;
lwc_result(1)=0;
qv_result(1)=qvb*1000;
ql_result(1)=0;
delta_qv(1)=0;

delta_z=0.001;
iz = 2;
while 1
    for i=1:12
        if (T1>TT(i)&&T1<=TT(i+1))
            Latent=LL(i)+(LL(i+1)-LL(i))/(TT(i+1)-TT(i))*(T1-TT(i));
            Latent=Latent*1000;
        end
    end
    Ra=(1+0.61*qv1)*Rd;
    c1=Latent*qv1/(Ra*T1);
    c2=Latent^2*qv1/(cp*Rv*T1^2);
    moisture_lap_rate=g/cp*(1+c1)/(1+c2);
    T2=T1-moisture_lap_rate*delta_z;

    c=g*P1/(Ra*T1);
    P2=P1-c*delta_z;
    Z2=Z1+delta_z;

    if abs(Z2-Zz(iz))>0.002
        T1=T2;
        Z1=Z2;
        P1=P2;
        for i=1:12
            if (T1>TT(i)&&T1<=TT(i+1))
                Latent=LL(i)+(LL(i+1)-LL(i))/(TT(i+1)-TT(i))*(T1-TT(i));
                Latent=Latent*1000;
            end

        es1=6.11*10^(3.45e-6*Latent*((T1-273.15)/T1));
        qv1=0.622*es1/(P1-es1) ;
        end
    else
        T_result=T2;
        Z_result=Z2;
        P_result=P2;

        for i=1:12
            if (T_result>TT(i)&&T_result<=TT(i+1))
                Latent=LL(i)+(LL(i+1)-LL(i))/(TT(i+1)-TT(i))*(T_result-TT(i));
                Latent=Latent*1000;
            end
        end

        es=6.11*10^(3.45e-6*Latent*((T_result-273.15)/T_result));
        qv=0.622*es/(P_result-es) ;

        ql=qvb-qv;
        e=qv/0.622*P_result/(qv/0.622+1);
        row=P_result*100/(Rd*T_result)*(1-e/P_result*(1-0.622));
        lwc=ql*row;
        delta_qv=qv-qvb;

        lwc_result(iz)=lwc*1000;
        qv_result(iz)=qv*1000;
        ql_result(iz)=ql*1000;
        delta_qv(iz)=delta_qv*1000;
        if iz == nz
            break
        end
        iz = iz+1;
    end
end
end
