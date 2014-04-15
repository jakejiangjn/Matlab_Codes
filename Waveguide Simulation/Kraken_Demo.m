clear all;close all;clc;
%% Generate field.flp
titleflp = 'field'; % to generate field.flp
SD = 10; % Source Elements' Depth in m
RD = linspace(0.0,200,201); % Each Source Receive Aarriy Elements' Depth in m
R = linspace(0.0,0.2,201); % Source Receive Aarries' Location in km
Pos = struct('s',struct('depth',SD),'r',struct('depth',RD,'range',R));
write_fieldflp( 'RA', Pos );
%% Generate env file
envfil = 'KrakenDemo'; % to generate MunkKraken.env
freq = 1e3; % frequency
% Boundary Codition Setting
Bdry.Top.Opt = 'NVW'; % Options of Top Boundary
Bdry.Top.Hs.alphaR = [];
Bdry.Top.Hs.betaR = [];
Bdry.Top.Hs.rho = [];
Bdry.Top.Hs.alphaI = [];
Bdry.Top.Hs.betaI = [];
Bdry.Bot.Opt = 'A'; % Options of Bottom Boundary
Bdry.Bot.HS.alphaR = 2000; % CPT -> P-wave Sound Speed of Bottom Boundary
Bdry.Bot.HS.betaR = 0; % CST -> S-wave Sound Speed of Bottom Boundary
Bdry.Bot.HS.rho = 2.0; % RHOT -> density of Bottom Boundary
Bdry.Bot.HS.alphaI = 0.3; % APT -> Attenuation of P-wave Bottom Boundary
Bdry.Bot.HS.betaI = 0; % AST -> Attenuation of S-wave Bottom Boundary
% Layer Media Setting
SSP.NMedia = 1; % Media Numbers between Top and Rigid Bottom 
Layer(1).z = linspace(0,200,201);
L = length( Layer(1).z );
Layer(1).alphaR = linspace( 1460,1450, L );
Layer(1).betaR = zeros( 1, L );
Layer(1).rho = ones( 1, L );
Layer(1).alphaI = zeros( 1, L );
Layer(1).betaI = zeros( 1, L );
SSP.N = 25000; % NMesh of each Layer -> integer array
SSP.sigma = zeros(1,1); % RMS roughness at the interface
SSP.depth = [0 200]; % Starting Depth of each Media Layer and the last element is the depth of Bottom
SSP.raw = Layer;
Beam = struct();
cInt = struct('Low',1400,'High',1600);
RMax = 0;
clear Layer L;
if exist([envfil '.env'], 'file')
    delete([envfil '.env']);
end;
write_env( envfil, 'KRAKEN', envfil, freq, SSP, Bdry, Pos, Beam, cInt, RMax, [] );
%write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
%% Run Kraken
kraken( envfil ); % krakenc( envfil );
%% Post processing
[PlotTitle,PlotType,~,alten,~,p] = read_shd( [envfil '.shd'] );
% Struct 'Pos' includes 3 structs 'Pos.theta', 'Pos.s'('Pos.s.depth' array of NSD*1) and
% 'Pos.r'('Pos.r.depth' array of NRD*1;'Pos.r.range' array of NR*1)
% size(p) = length(SD)*length(RD)*length(R)
Green = squeeze(p);  Green = abs(Green);  Green = 10*log10(Green/max(max(Green)));
%%
figure(1); pcolor(R,RD,Green); % pcolor(Pos.r.range,Pos.r.depth,Green);
xlabel('range/km');  ylabel('depth/m');
shading interp; % interpolation for smooth display
axis ij; % Set the Origin Point to the up-left corner
caxis([-30 -0]); % Set the display range
colorbar;
%%
figure(2);
plotshd( [envfil '.shd'] ); colorbar;
