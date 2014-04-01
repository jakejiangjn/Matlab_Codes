% A Bellhop simulation for the Mogan Lake Experiment on Dec 2013.
clear all;close all;clc;
%%
z = linspace(0,20,201);
ssp = [linspace(1460,1459.48,20) linspace(1459.48,1459.73,181)]; ssp = round( ssp*1000 ) / 1000;
%% Generate field.flp
SD = 10; % Source Elements' Depth in m
RD = [10 20]; % Each Source Receive Aarriy Elements' Depth in m
R = linspace(0,0.3,301); % Source Receive Aarries' Location in km
Pos = struct('s',struct('depth',SD),'r',struct('depth',RD,'range',R));
write_fieldflp( 'RA', Pos );
%% Environment Preparation
envfil = 'LakeB'; % env file name
freq = 6e3; % frequency
% Boundary Structure
Bdry = struct('Top',struct('Opt',{'NVW'},'HS',struct),'Bot',struct('Opt',{'A'},'HS',struct));
Bdry.Bot.HS = struct('alphaR',2000,'betaR',0,'rho',2.0,'alphaI',0.5,'betaI',0);
% Media Structure
Layer = struct( 'z',z,'alphaR',ssp,'betaR',zeros(1,length(z)),'rho',ones(1,length(z)),...
    'alphaI',zeros(1,length(z)),'betaI',zeros(1,length(z)) );
SSP = struct( 'NMedia',1,'N',2500,'sigma',0,'depth',[0 20],'raw',Layer );
Beam = struct( 'RunType',{'R'}, 'Nrays',50, 'Nbeams',50, 'Ibeams',[],'alpha',[-50 50],...
    'deltas',0, 'Box',struct('z',25,'r',R(end)), 'epmult',[], 'rLoop',[], 'Nimage',[], 'Ibwin',[]  );
cInt = struct('Low',1400,'High',2000);  RMax = 0;
clear Layer;
%% Generate LakeB.env
if exist([envfil '.env'], 'file')
    delete([envfil '.env']);
end;
write_env( envfil, 'BELLHOP', 'MoganLakeX', freq, SSP, Bdry, Pos, Beam, cInt, RMax, [] );
% varargin =[] indicates to delete/overwrite the existing file,
% otherwise it will concatenate the existing one.
%%
bellhop( envfil );  plotray( envfil );
%%
Beam.RunType = 'A';
if exist([envfil '.env'], 'file')
    delete([envfil '.env']);
end;
write_env( envfil, 'BELLHOP', 'MoganLakeX', freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
bellhop( envfil );
[ Arr, ~ ] = read_arrivals_asc( [envfil '.arr'] );
% [ Arr, Pos ] = read_arrivals_asc( [envfil '.arr'],50 )