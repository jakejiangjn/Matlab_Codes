fclose('all');   clear all;	close all;	clc;
%% Generate field.flp
titleflp = 'field'; % to generate field.flp
SD = 10; % Source Elements' Depth in m
RD = linspace(0.0,100,101); % Each Source Receive Aarry Elements' Depth in m
R = 4; % Source Receive Aarries' Location in km
Pos.theta = []; Pos.s.depth = SD;
Pos.r.depth = RD; Pos.r.range = R;
rProf = [0.0 1.0 2.0];
NProf = length(rProf);
write_flp( 'RC', NProf, rProf, Pos );
%% Generate env file
envfil = 'BellhopDemo'; % to generate KrakenDemo.env
freq = 800; % frequency
% Boundary Codition Setting
Bdry.Top.Opt = 'QVW'; % Options of Top Boundary
Bdry.Top.Hs.alphaR = []; Bdry.Top.Hs.betaR = [];
Bdry.Top.Hs.rho = [];
Bdry.Top.Hs.alphaI = []; Bdry.Top.Hs.betaI = [];
Bdry.Bot.Opt = 'A*'; % Options of Bottom Boundary
Bdry.Bot.HS.alphaR = 2000; % CPT -> P-wave Sound Speed of Bottom Boundary(m/s)
Bdry.Bot.HS.betaR = 0; % CST -> S-wave Sound Speed of Bottom Boundary(m/s)
Bdry.Bot.HS.rho = 2.0; % RHOT -> density of Bottom Boundary(g/cm^3)
Bdry.Bot.HS.alphaI = 0.3; % APT -> Attenuation of P-wave Bottom Boundary
Bdry.Bot.HS.betaI = 0; % AST -> Attenuation of S-wave Bottom Boundary
% Layer Media Setting
SSP.NMedia = 1; % Media Numbers between Top and Rigid Bottom
Layer(1).z = linspace(0,100,101); L = length( Layer(1).z );
Layer(1).alphaR = linspace( 1480,1485, L );
Layer(1).betaR = zeros( 1, L );
Layer(1).rho = ones( 1, L );
Layer(1).alphaI = zeros( 1, L );
Layer(1).betaI = zeros( 1, L );
SSP.N = [8000 2000]; % NMesh of each Layer -> integer array
% NMesh(n) = ( >10 ) * Layer(n).z(end) * max(freq) / min(Layer(n).alphaR)
SSP.sigma = zeros(1,1); % RMS roughness at the interface
SSP.depth = [0 100]; % Starting Depth of each Media Layer and the last element is the depth of Bottom
SSP.raw = Layer;
Beam.RunType = 'R'; % RunType R / A
Beam.Nrays = 9; % just same as Beam.Nbeams for safety
Beam.Nbeams = 9; % the numbers of Rays/Beams to trace (use 0 to have the program calculate a value automatically)
Beam.Ibeams = [];
Beam.alpha = [-45 45]; % Beam angles (negative angles toward surface)
Beam.deltas = 0; % the step size used for tracing the rays(m)
Beam.Box.z = 1.1*RD(end); % the maxium depth to trace a ray(m)
Beam.Box.r = 1.01*R; % the maxium range to trace a ray(km)
% the following 4 parameters are for Cerveny-style Gaussian beams when RunType contains 2 letters->Beam.RunType(2:2) = 'G'/'B'/'S'
Beam.epmult = []; Beam.rLoop = []; Beam.Nimage = []; Beam.Ibwin = [];
RMax = 0.0;
cInt.Low = 1400; cInt.High = 15000;
c = zeros( L, NProf );
clear Layer L;
if exist([envfil '.env'], 'file')
    delete([envfil '.env']);
end;
write_env( envfil, 'BELLHOP', 'Range Dependent', freq, SSP, Bdry, Pos, Beam, cInt, RMax );
c(:,1) = SSP.raw(1).alphaR;
for n = 2:NProf
    c(:,n) = c(:,n-1) + 50;
end;
writessp( [envfil, '.ssp'], rProf, c );
writebty( [envfil, '.bty'], 'L', [rProf;50,75,100]' );
% syntax: write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
%%
bellhop( envfil );
%%
plotray( envfil );
%%
% [ Arr, Pos ] = read_arrivals_asc( [envfil '.arr'] );
