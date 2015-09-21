fclose('all');   clear all;	close all;	clc;
%% Generate field.flp
titleflp = 'field'; % to generate field.flp
SD = 10; % Source Elements' Depth in m
RD = linspace(0.0,100,101); % Each Source Receive Aarry Elements' Depth in m
R = linspace(0.0,2,201); % Source Receive Aarries' Location in km
Pos.theta = []; Pos.s.depth = SD;
Pos.r.depth = RD; Pos.r.range = R;
rProf = [0.0 1.0 2.0];
NProf = length(rProf);
write_flp( 'RC', NProf, rProf, Pos );
%% Generate env file
envfil = 'KrakenDemo_rd'; % to generate KrakenDemo_rd.env
freq = 800; % frequency
% Boundary Codition Setting
Bdry.Top.Opt = 'NVW'; % Options of Top Boundary
Bdry.Top.Hs.alphaR = []; Bdry.Top.Hs.betaR = [];
Bdry.Top.Hs.rho = [];
Bdry.Top.Hs.alphaI = []; Bdry.Top.Hs.betaI = [];
Bdry.Bot.Opt = 'A'; % Options of Bottom Boundary
Bdry.Bot.HS.alphaR = 2000; % CPT -> P-wave Sound Speed of Bottom Boundary(m/s)
Bdry.Bot.HS.betaR = 0; % CST -> S-wave Sound Speed of Bottom Boundary(m/s)
Bdry.Bot.HS.rho = 2.0; % RHOT -> density of Bottom Boundary(g/cm^3)
Bdry.Bot.HS.alphaI = 0.3; % APT -> Attenuation of P-wave Bottom Boundary
Bdry.Bot.HS.betaI = 0; % AST -> Attenuation of S-wave Bottom Boundary
% Layer Media Setting
SSP.NMedia = 2; % Media Numbers between Top and Rigid Bottom
Layer(1).z = linspace(0,100,101); L = length( Layer(1).z );
Layer(1).alphaR = linspace( 1430,1435, L );
Layer(1).betaR = zeros( 1, L );
Layer(1).rho = ones( 1, L );
Layer(1).alphaI = zeros( 1, L );
Layer(1).betaI = zeros( 1, L );
Layer(2).z = linspace(100,120,21); L = length( Layer(2).z );
Layer(2).alphaR = linspace( 1600,1600, L );
Layer(2).betaR = zeros( 1, L );
Layer(2).rho = ones( 1, L )*1.9;
Layer(2).alphaI = ones( 1, L )*0.5;
Layer(2).betaI = zeros( 1, L );
SSP.N = [8000 2000]; % NMesh of each Layer -> integer array
% NMesh(n) = ( >10 ) * Layer(n).z(end) * max(freq) / min(Layer(n).alphaR)
SSP.sigma = zeros(1,2); % RMS roughness at the interface
SSP.depth = [0 100 120]; % Starting Depth of each Media Layer and the last element is the depth of Bottom
SSP.raw = Layer;
Beam = struct();
RMax = 0.0;
cInt.Low = 1400; cInt.High = 2000;
clear Layer L;
if exist([envfil '.env'], 'file')
    delete([envfil '.env']);
end;
for n = 1:NProf
	SSP.raw(1).alphaR = SSP.raw(1).alphaR + 50;
	write_env( envfil, 'KRAKENC', ['Range Dependent #' int2str(n)], freq, SSP, Bdry, Pos, Beam, cInt, RMax, 'a' );
end;
% syntax: write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
%% Run Kraken
% generate KrakenDemo.shd and KrakenDemo.mod
krakenc( envfil );
%% Post processing
[~,~,~,atten,~,p] = read_shd( [envfil '.shd'] );
% Struct 'Pos' includes 3 structs 'Pos.theta', 'Pos.s'('Pos.s.depth' array of NSD*1) and
% 'Pos.r'('Pos.r.depth' array of NRD*1;'Pos.r.range' array of NR*1)
% size(p) = length(SD)*length(RD)*length(R)
Green = squeeze(p); Green = abs(Green); Green = 10*log10(Green/max(max(Green)));
%%
figure(1); pcolor(R,RD,Green); % pcolor(Pos.r.range,Pos.r.depth,Green);
xlabel('range/km'); ylabel('depth/m');
shading interp; % interpolation for smooth display
axis ij; % Set the Origin Point to the up-left corner
caxis([-30 0]); % Set the display range
colorbar;
%%
figure(2);
plotshd( [envfil '.shd'] ); colorbar;
%%
for n = 1:NProf
    Modes(n) = read_modes( [envfil '.mod'] );
end;
