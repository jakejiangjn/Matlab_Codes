clear all;close all;clc;
%%
z = linspace(0,2000,201);
ssp = 1500.0*(1.0+ 0.00737*(2/1300*(z-1300)-1+exp(2/1300*(z-1300))) ); % Munk Profile
% ssp = linspace( 1500,1600,length(z) ); % negative gradient
%%
titleflp = 'field'; % to generate field.flp
SD = 10; % Source Elements' Depth in m
RD = linspace(0.0,200,2001); % Each Source Receive Aarriy Elements' Depth in m
R = linspace(0.0,20,2001); % Source Receive Aarries' Location in km
Pos = struct('s',struct('depth',SD),'r',struct('depth',RD,'range',R));
write_fieldflp( 'RA', Pos );
%%
envfil = 'MunkKraken'; % to generate MunkKraken.env
freq = 6e3; % frequency
Bdry = struct('Top',struct('Opt',{'NVW'},'HS',struct),'Bot',struct('Opt',{'A'},'HS',struct));
Bdry.Bot.HS = struct('alphaR',2000,'betaR',0,'rho',2.0,'alphaI',0.3,'betaI',0);
%%
NMedia = 1;
Layer(1) = struct( 'z',z,'alphaR',ssp,'betaR',zeros(1,length(z)),'rho',ones(1,length(z)),...
    'alphaI',zeros(1,length(z)),'betaI',zeros(1,length(z)) );
SSP = struct( 'NMedia',NMedia,'N',[2500],'sigma',zeros(1,1),'depth',[0 2000],'raw',Layer );
Beam = struct();
cInt = struct('Low',1400,'High',1600);  RMax = 0;  varargin = [];
%%
if exist([envfil '.env'], 'file')
    delete([envfil '.env']);
end;
write_env( envfil, 'KRAKEN', 'Munk Profile', freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
%write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
%%
kraken( envfil ); % krakenc( envfil );
[PlotTitle,PlotType,~,alten,~,p] = read_shd( [envfil '.shd'] );
% Struct 'Pos' includes 3 structs 'Pos.theta', 'Pos.s'('Pos.s.depth' array of NSD*1) and
% 'Pos.r'('Pos.r.depth' array of NRD*1;'Pos.r.range' array of NR*1)
% size(p) = length(SD)*length(RD)*length(R)
Green = squeeze(p);  Green = abs(Green);  Green = 10*log10(Green/max(max(Green)));
%%
figure(1); pcolor(R,RD,Green); % pcolor(Pos.r.range,Pos.r.depth,Green);
xlabel('range/km');  ylabel('depth/m');
shading interp; % 在网格片内采用颜色插值处理，得出表面图显得最光滑
axis ij; % 设置坐标原点在左上角
caxis([-30 -0]); % 设置颜色的范围为指定的最小值和最大值 对应的最大最小值映射到该区间上 
colorbar;
%%
figure(2);
plotshd( [envfil '.shd'] ); colorbar;