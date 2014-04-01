% A Karken based Reverberation Simulation 
% (Environment Paramenters are set based on Mogan Lake Experiment on Dec 2013.)
% This program is based on
% Ellis, Dale D, A shallow-water normal-mode reverberation model, JASA, 97(5):2804-2814, 1995.
clear all;close all;clc;
%%
z = linspace(0,20,201);
ssp = [linspace(1460,1459.48,20) linspace(1459.48,1459.73,181)]; ssp = round( ssp*1000 ) / 1000;
%% Generate field.flp
SD = 10; % Source Elements' Depth in m
RD = 20; % Each Source Receive Aarriy Elements' Depth in m
R = linspace(0,0.3,301); % Source Receive Aarries' Location in km
Pos = struct('s',struct('depth',SD),'r',struct('depth',RD,'range',R));
% write_fieldflp( 'RA', Pos );
%% Environment Preparation
envfil = 'LakeK'; % to generate MunkKraken.env
freq = 6e3; % frequency
% Boundary Structure
Bdry = struct('Top',struct('Opt',{'NVW'},'HS',struct),'Bot',struct('Opt',{'A'},'HS',struct));
Bdry.Bot.HS = struct('alphaR',2000,'betaR',0,'rho',2.0,'alphaI',0.3,'betaI',0);
% Media Structure
NMedia = 2;
Layer(1) = struct( 'z',z,'alphaR',ssp,'betaR',zeros(1,length(z)),'rho',ones(1,length(z)),...
    'alphaI',zeros(1,length(z)),'betaI',zeros(1,length(z)) );
Layer(2) = struct( 'z',[20 22],'alphaR',[1460 1460],'betaR',[0 0],'rho',[1.5 1.5],...
    'alphaI',[0.7 0.7],'betaI',[0 0] );
SSP = struct( 'NMedia',NMedia,'N',[2500 2500],'sigma',zeros(1,2),'depth',[0 20 22],'raw',Layer );
Beam = struct();
cInt = struct('Low',1400,'High',2000);  RMax = 0;
clear Layer;
%% Generate LakeK.env
% if exist([envfil '.env'], 'file')
%     delete([envfil '.env']);
% end;
% write_env( envfil, 'KRAKEN', 'MoganLakeX', freq, SSP, Bdry, Pos, Beam, cInt, RMax, [] );
% syntax: write_env( envfil, model, TitleEnv, freq, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
%% Run Kraken
% kraken( envfil ); % krakenc( envfil );
% [status,result] = system( ['kraken <' envfil '.env '] );
% if ( status ~= 0 )
%     disp( 'kraken error!' );
% else
%     disp( 'kraken done!' );
% end;
%% Read shade file for pressure field
% [PlotTitle,PlotType,~,atten,~,p] = read_shd( [envfil '.shd'] );
% [ PlotTitle, PlotType, freq, atten, Pos, p ] = read_shd( [envfil '.shd'] );
% plotshd( [envfil '.shd'] );
% Notice:
% Struct 'Pos' includes 3 structs 'Pos.theta', 'Pos.s'('Pos.s.depth' array of NSD*1) and
% 'Pos.r'('Pos.r.depth' array of NRD*1;'Pos.r.range' array of NR*1)
% size(p) = length(SD)*length(RD)*length(R)
% Notice END.
plotssp( envfil );
%% Read Modes
Modes = read_modes( [envfil '.mod'] );
% plotmode( [envfil '.mod'] );
%% Reverberation Generation
zs = find( Modes.z == SD );
zr = find( Modes.z == z(end) );
miu = 10^(-18/10); % The Scatter Strength of mud.
cp = 2*pi * freq ./ real(Modes.k); % Modes.k = column vector
cp(cp<ssp(end)) = ssp(end);
cg = ssp(end)^2 ./ cp; % must be in lossless medium
theta_g = acos( ssp(end)./cp ); % grazing angle
S = miu * sin( theta_g ) *  sin( theta_g.' );
S = S.^(0.5);
I0 = 1; % Stimulation Intensity
T = 10e-3; % Singal Pulse Duration
fs = 50e3;
R = R * 1000;
%% Varialbes Memory Pre Allocation
rever_coh = zeros( size(R) );   rever_inc = zeros( size(R) ); % Pressure value at time/range
% Rever_coh = zeros( size(R) );   Rever_inc = zeros( size(R) ); % Intensity at time/range
for n = 9:length(R)
%     p = 1j*pi*Modes.phi(zs,:).*Modes.phi(zr,:).*besselh(0,1,abs(imag(Modes.k')*R(n)));
    p = 1j*pi*Modes.phi(zs,:).*Modes.phi(zr,:).*besselh(0,1,(Modes.k')*R(n));
    rever_coh(n) = (p*S*p.') * sqrt( pi*R(n) * T*ssp(end) * I0 );
    rever_inc(n) = sqrt( (abs(p).^2)*(S.^2)*(abs(p').^2)  * pi*R(n) * T*ssp(end) * I0 );
%     rever_inc(n) = (p*S*p') *  sqrt( pi*R(n) * T*ssp(end) * I0 );
%     Rever_coh(n) = abs( (p*S*p.')^2 ) * pi*R(n) * T*ssp(end) * I0;
%     Rever_inc(n) = (abs(p).^2)*(S.^2)*(abs(p').^2)  * pi*R(n) * T*ssp(end) * I0;
end;
%% Calculate the delay points
delay_coh = -atan2( imag(rever_coh), real(rever_coh) );
delay_coh(delay_coh<0) = delay_coh(delay_coh<0) + 2*pi;
delay_coh = round( (delay_coh/2/pi/freq+(R.^2+100).^0.5*2/ssp(end))*fs );
delay_inc = -atan2( imag(rever_inc), real(rever_inc) );
delay_inc(delay_inc<0) = delay_inc(delay_inc<0) + 2*pi;
delay_inc = round( (delay_inc/2/pi/freq+(R.^2+100).^0.5*2/ssp(end))*fs );
% Amplitude
rever_coh = abs( rever_coh );   rever_inc = abs( rever_inc );
%%
s = 0:(1/fs):(T-1/fs);  s = cos( 2*pi*freq*s) .* hanning( length(s) ).';
r = [s zeros( 1, max( delay_coh(end), delay_inc(end) ) + length(s) )];
rec_coh = zeros( size(r) ); rec_inc = zeros( size(r) );
for n = 1:length(R)
    rec_coh = rec_coh + circshift( r, [0 delay_coh(n)]) * rever_coh(n);
    rec_inc = rec_inc + circshift( r, [0 delay_inc(n)]) * rever_inc(n);
end;
flt_co = fir1( 128, [5e3 11e3]/fs*2 );
rec_coh = ( filtfilt(flt_co,1,rec_coh.') ).';
rec_inc = ( filtfilt(flt_co,1,rec_inc.') ).';
save( 'Rever_coh.mat',  'rec_coh', 'delay_coh' );
save( 'Rever_inc.mat',  'rec_inc', 'delay_inc' );
%% Plot Work
figure(1);
subplot(2,1,1);plot( (0:length(r)-1)/fs, rec_coh );
title('Reverberation simulated by Coheret Calculation'); xlabel('Time(s)'); ylabel('Amplitude');
subplot(2,1,2);plot( (0:length(r)-1)/fs, rec_inc );
title('Reverberation simulated by Incoheret Calculation'); xlabel('Time(s)'); ylabel('Amplitude');
figure(2);
subplot(2,1,1);spectrogram(rec_coh,512,500,1024,fs);colorbar;
title('Reverberation simulated by Coheret Calculation'); 
subplot(2,1,2);spectrogram(rec_inc,512,500,1024,fs);colorbar;
title('Reverberation simulated by Incoheret Calculation'); 
figure(3);
plot( R,10*log10( abs(rever_coh).^2 ),'b',  R,10*log10( abs(rever_inc).^2 ),'r--');
legend('Coherent','Incoherent');
title('Bottom Reverberation for Shallow Water');btitle();
xlabel('Range(m)');ylabel('Reverberation Level(dB)');