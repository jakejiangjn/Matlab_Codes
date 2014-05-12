clear all;close all;clc;
%% Generate field.flp
titleflp = 'field'; % to generate field.flp
SD = 10; % Source Elements' Depth in m
RD = 50; % Each Source Receive Aarriy Elements' Depth in m
R = 2; % Source Receive Aarries' Location in km
Pos = struct('s',struct('depth',SD),'r',struct('depth',RD,'range',R));
write_fieldflp( 'RA', Pos );
%% Generate env file
envfil = 'BellhopDemo';
freq = 1e3; % frequency
% Boundary Codition Setting
Bdry.Top.Opt = 'SVW'; % Options of Top Boundary
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
Layer(1).z = linspace(0,200,201);	L = length( Layer(1).z );
Layer(1).alphaR = linspace( 1460,1450, L );
Layer(1).betaR = zeros( 1, L );
Layer(1).rho = ones( 1, L );
Layer(1).alphaI = zeros( 1, L );
Layer(1).betaI = zeros( 1, L );
SSP.N = 2500; % NMesh of each Layer -> integer array( size(SSP.N)=1,SSP.NMedia )
SSP.sigma = zeros(1,1); % RMS roughness at the interface
SSP.depth = [0 200]; % Starting Depth of each Media Layer and the last element is the depth of Bottom( size(SSP.N)=1,SSP.NMedia+1 )
SSP.raw = Layer;
clear Layer L;
Beam.RunType = 'R'; % RunType
Beam.Nrays = 20; % just same as Beam.Nbeams for safety
Beam.Nbeams = 20; % the numbers of Rays/Beams to trace (use 0 to have the program calculate a value automatically)
Beam.Ibeams = [];
Beam.alpha = [-40 40]; % Beam angles (negative angles toward surface)
Beam.deltas = 0; % the step size used for tracing the rays(m)
Beam.Box.z = 250; % the maxium depth to trace a ray(m)
Beam.Box.r = R(end); % the maxium range to trace a ray(km)
% the following 4 parameters are for Cerveny-style Gaussian beams when RunType contains 2 letters->Beam.RunType(2:2) = 'G'/'B'/'S'
Beam.epmult = []; Beam.rLoop = []; Beam.Nimage = []; Beam.Ibwin = [];
cInt = struct('Low',1400,'High',1600);	RMax = 0;
clear Layer L;
if exist([envfil '.env'], 'file')
    delete([envfil '.env']);
end;
%% Generate MunkBellhop.env
if exist([envfil '.env'], 'file')
    delete([envfil '.env']);
end;
%%
if Beam.RunType == 'R' % generates a ray file to plot
	write_env( envfil, 'BELLHOP', 'Munk Profile', freq, SSP, Bdry, Pos, Beam, cInt, RMax, [] );
	% write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
	bellhop( envfil );  plotray( envfil );
else
	Beam.RunType = 'A'; % generates an complex amplitude-delay(sec) file for usage
	write_env( envfil, 'BELLHOP', 'Munk Profile', freq, SSP, Bdry, Pos, Beam, cInt, RMax, [] );
	bellhop( envfil );
	[ Arr, ~ ] = read_arrivals_asc( [envfil '.arr'] );
end;
