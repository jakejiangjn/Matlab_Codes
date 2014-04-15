close all;clear all;clc;
% parameter setup
% water layer
param.H1 = 100;         % depth of water（m）
param.c1 = 1600;        % speed in the water（m/s）
param.rho1 = 1;         % rho in the water 密度(g/cm^3)
param.alpha1 = 0;       % attention in the water

% sand layer
param.H2 = 110;          % depth of sand
param.c2_1 = 1600;      % speed in the sand
param.c2_2 = 1700;
param.rho2 = 1.9;       % rho in the sand
param.alpha2 = 0.1;     % attention in the sand

param.range = [950 1050];    % experiment

freq = 1000;             % 信号频率（Hz）

sensor.M = 20;
sensor.zr = linspace(0.1 ,param.H1-0.1, sensor.M);    %sensor location 阵元位置
sensor.xr = zeros(size(sensor.zr ));

num_source = 1;
%source.zs = 0.3+(param.H1-0.6)*rand(num_source,1);  % [0.3~param.H1-0.3]区间均匀分布
%source.xs = param.range(1)+(param.range(end)-param.range(1))*rand(num_source,1);      % [param.range(1)+1~param.range(end)-1]区间均匀分布
source.zs = 40;
source.xs = 1000;
% kraken
NMedia = 2;
%---------------------------------------layer--------------------------------------------------------%
Z_axs = 80;
Z1 = linspace(0,param.H1,1000);
SSP1 = 1500*(1+0.00737.*(2*(Z1-Z_axs)./Z_axs-1+exp(-2*(Z1-Z_axs)./Z_axs)));
% Z1 = [0 param.H1];
% SSP1 =  [param.c1 param.c1];
Z2 =  [param.H1 param.H2];
SSP2 =  [param.c2_1 param.c2_2];
Layer(1) = struct('Depth',param.H1,'Z',Z1,'SSP',SSP1,'rho',param.rho1,'alpha',param.alpha1);
Layer(2) = struct('Depth',param.H2,'Z',Z2,'SSP',SSP2,'rho',param.rho2,'alpha',param.alpha2);
%---------------------------------------bottom-----------------------------------------------------------%
%Bottom = struct('Depth',param.H2,'c',param.c3,'rho',param.rho3,'alpha',param.alpha3);
Bottom = struct('Depth',param.H2,'cp',1900,'rho',2,'alpha',0.3);
%--------------------------------------------------------------------------------------------------%
titleenv = 'pekeris';
titleflp = 'field';
p_true = 0;
for ii=1:num_source
    NPROF = 1;
    RD = sensor.zr;           % the receiver depth
    NRD = length(RD);         % the number of the receiver depth
    SD = source.zs(ii);       % the source depth
    NSD = length(SD);         % the number of source depth
    R = [0:source.xs(ii)]/1000;   % the receiver ranges (km)
    RPROF = 0;
    NR = length(R);           % number of receiver ranges

    WriteFieldflp( titleflp,NPROF,RPROF,NR, R, NSD, SD, NRD, RD );             % field file
    Writeenvfile( titleenv, freq, NMedia, Layer,Bottom, NSD,SD,RD); % env file
    kraken(titleenv);
    [ PlotTitle, PlotType, freq, atten, Pos, p_true1 ] = read_shd( [titleenv,'.shd']) ;
    p_true  = p_true+squeeze( p_true1 );
end
R_x = p_true(:,end)*p_true(:,end)';  % 真实的数据协方差矩阵
plotshd('pekeris.shd')
%% 匹配场搜索
dz = 1; 
dx = 1;
SD = 1: dz: param.H1-1;         % the source depth
NSD = length(SD);               % the number of source depth
R = [param.range(1) :dx: param.range(end)]/1000;                % the receiver ranges (km)
NR = length(R);                 % number of receiver ranges

WriteFieldflp( titleflp,NPROF,RPROF, NR, R, NSD, SD, NRD, RD );             % field file
Writeenvfile( titleenv, freq, NMedia, Layer,Bottom, NSD,SD,RD); % env file
kraken(titleenv);
[ PlotTitle, PlotType, freq, atten, Pos, p_simu ] = read_shd( [titleenv,'.shd']) ;
p_simu = reshape( p_simu,[ NSD, sensor.M , NR ]); 

BartMFP = zeros(NR,NSD);
tic
for ii = 1:NSD
    G = squeeze(p_simu(ii,:,:));
    n_G  = sqrt(sum(abs(G).*abs(G)));
    G = G*diag(1./n_G);                     % normlize Phi
    BartMFP(:,ii) = diag(abs( G'*R_x*G ));
end
toc
BartMFP = 10*log10(BartMFP./max(max(BartMFP)) );
find(BartMFP ==max(max(BartMFP)) )
figure(2); 
pcolor( R*1000 ,SD  , (BartMFP') );
shading interp; colorbar;set(colorbar,'box','off');
xlabel('Range (m)','FontSize',10); ylabel('Depth (m)','FontSize',10);
axis ij; box off;set(gca,'FontSize',10);

for ii = 1:num_source
    hold on;
    p_marker = plot( source.xs(ii), source.zs(ii), 'o');
    set( p_marker,'Color','k','MarkerSize',10 );
end
hold off;title('Bartlett MFP');
%% 时反处理
param.range = [990 1010];    % experiment

sensor.zr = linspace(10 ,param.H1, 10);    %sensor location 阵元位置
sensor.xr = 0;
sensor.M = length(sensor.zr );

num_source = 20;
source.zs = linspace(0.1 ,param.H1-0.1, num_source);
source.xs = [param.range(1):1:param.range(2)];

p_re = 1;

for ii = 1:num_source
    RD = sensor.zr;           % the receiver depth
    NRD = length(RD);         % the number of the receiver depth
    SD = source.zs(ii);       % the source depth
    NSD = length(SD);         % the number of source depth
    R = source.xs/1000;       % the receiver ranges (km)
    NR = length(R);           % number of receiver ranges

    WriteFieldflp( titleflp,NPROF,RPROF, NR, R, NSD, SD, NRD, RD );             % field file
    Writeenvfile( titleenv, freq, NMedia, Layer,Bottom, NSD,SD,RD); % env file
    kraken(titleenv);
    [PlotTitle, PlotType, freq, atten, Pos, p_re1] = read_shd( [titleenv,'.shd']) ;
    p_re  = p_re+squeeze( p_re1 ).*conj(p_true(ii,end));
end
%%
% aa = abs(p_re)./max(max(abs(p_re))); 
p_re1 = 20*log10(abs(p_re)./max(max(abs(p_re))));
figure(1); 
pcolor( source.xs,sensor.zr,p_re1);
shading interp; colorbar;set(colorbar,'box','off');
xlabel('Range (m)','FontSize',10); ylabel('Depth (m)','FontSize',10);
axis ij; box off;set(gca,'FontSize',10);
title('TRM');colorbar
