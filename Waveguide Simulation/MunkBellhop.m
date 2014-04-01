clear all;close all;clc;
%%
z = linspace(0,2000,201);
ssp = 1500.0*( 1.0+ 0.00737*(2/1300*(z-1300)-1+exp(2/1300*(z-1300))) ); % Munk Profile
% ssp = linspace( 1500,1600,length(z) ); % negative gradient
%% Generate field.flp
SD = 10; % Source Elements' Depth in m
RD = 50; % Each Source Receive Aarriy Elements' Depth in m
R = 2; % Source Receive Aarries' Location in km
Pos = struct('s',struct('depth',SD),'r',struct('depth',RD,'range',R));
write_fieldflp( 'RA', Pos );
%%
envfil = 'MunkBellhop';
freq = 2e3; % frequency
Bdry = struct('Top',struct('Opt',{'SVW'},'HS',struct),'Bot',struct('Opt',{'A'},'HS',struct));
Bdry.Bot.HS = struct('alphaR',2000,'betaR',0,'rho',2.0,'alphaI',0.3,'betaI',0);
%%
NMedia = 1;
Layer(1) = struct( 'z',z,'alphaR',ssp,'betaR',zeros(1,length(z)),'rho',ones(1,length(z)),...
    'alphaI',zeros(1,length(z)),'betaI',zeros(1,length(z)) );
SSP = struct( 'NMedia',NMedia,'N',[2500],'sigma',zeros(1,1),'depth',[0 2000],'raw',Layer );
Beam = struct( 'RunType',{'R'}, 'Nrays',20, 'Nbeams',20, 'Ibeams',[],'alpha',[-40 40],...
    'deltas',0, 'Box',struct('z',2500,'r',R(end)), 'epmult',[], 'rLoop',[], 'Nimage',[], 'Ibwin',[]  );
cInt = struct('Low',1400,'High',1600);  RMax = 0;  varargin = [];
%% Generate MunkBellhop.env
if exist([envfil '.env'], 'file')
    delete([envfil '.env']);
end;
write_env( envfil, 'BELLHOP', 'Munk Profile', freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
%write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
%%
bellhop( envfil );  plotray( envfil );
%%
Beam.RunType = 'A';
write_env( envfil, 'BELLHOP', 'Munk Profile', freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
bellhop( envfil );
[ Arr, ~ ] = read_arrivals_asc( [envfil '.arr'] );