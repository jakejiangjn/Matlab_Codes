clear all;close all;clc;
%% Basic Parameters
fs = 20e3;  fl = 2e3;	fh = 6e3;   B = fh-fl; % Sample rate; Low/High cutoff frequency; Bandwith
T = 10e-3;	t = ( 0:1/fs:T-1/fs ).';   N = length(t); % Time axis
beg = fix( N/8 );
k = B / t(N-2*beg+1);
r_tag = 2.5;  d_tag = 75;
%% Signal Generation
m = 0.5*t(2:N-2*beg+1).^2;
m = [ zeros(beg,1); m; m(end)+(1:beg).'/fs*t(N-2*beg+1) ];
N_pre = fix( N/4 );
% Window
w = hann(N_pre*2);   w = [w(1:N_pre);ones(N-N_pre*2,1);w(N_pre+1:end)];
% LFM
s = cos( 2*pi* (fl*t+k*m) ) .* w; % real( exp( 2j**pi* (fl*t+k*s) ) );
clear m w N_pre beg N;
%% env File Configuration
z = linspace(0,200,201);
ssp = 1480.0*(1.0+ 0.00737*((z-50)/50-1+exp((50-z)/50)) ); % Munk Profile
SD = 59; % Source Elements' Depth in m
RD = linspace(1,199,199); % Each Source Receive Aarriy Elements' Depth in m
R = linspace(0.5,5.5,201); % Source Receive Aarries' Location in km
Pos = struct('s',struct('depth',SD),'r',struct('depth',RD,'range',R));
r_index = find( abs(R-r_tag) <= min(abs(R-r_tag)), 1 );
d_index = find( abs(RD-d_tag) <= min(abs(RD-d_tag)), 1 );
envfil = 'Local'; % to generate env file
Bdry = struct('Top',struct('Opt',{'SVW'},'HS',struct),'Bot',struct('Opt',{'A'},'HS',struct));
Bdry.Bot.HS = struct('alphaR',2000,'betaR',0,'rho',2.0,'alphaI',0.3,'betaI',0);
NMedia = 1;
Layer(1) = struct( 'z',z,'alphaR',ssp,'betaR',zeros(1,length(z)),'rho',ones(1,length(z)),...
    'alphaI',zeros(1,length(z)),'betaI',zeros(1,length(z)) );
SSP = struct( 'NMedia',NMedia,'N',[2500],'sigma',zeros(1,1),'depth',[0 200],'raw',Layer );
Beam = struct( 'RunType',{'A'}, 'Nrays',20, 'Nbeams',20, 'Ibeams',[],'alpha',[-45 45],...
    'deltas',0, 'Box',struct('z',210,'r',R(end)), 'epmult',[], 'rLoop',[], 'Nimage',[], 'Ibwin',[]  );
cInt = struct('Low',1400,'High',1600);  RMax = 0;  varargin = 0;
%% Time-Delay & Attenation Generation
% % write_fieldflp( 'RA', Pos );
% % if exist([envfil '.env'], 'file')
% %     delete([envfil '.env']);
% % end;
% % write_env( envfil, 'BELLHOP', 'Localization', fl, SSP, Bdry, Pos, Beam, cInt, RMax, varargin );
% % bellhop( envfil );
% ( Nrr, Narrmx, Nrd, Nsd )
clear z ssp RMax varargin;
%% Read Generated Para
[ Arr, ~ ] = read_arrivals_asc( [envfil '.arr'], 20 );
A = Arr.A(r_index,:,d_index)*100;   delay = Arr.delay(r_index,:,d_index);
thres = 0.01*max(abs(A));
delay = delay( abs(A) > thres );    A = A( abs(A) > thres );
delay = round( ( delay - min(delay) ) * fs );
%% Received Signal Generation by Circle Shift
% s_temp = [ s ; zeros( round(max(delay)+length(s)/2), 1 ) ];
% r = zeros( size(s_temp) );
% for n = 1:length(A)
%     r = r + circshift( real( hilbert(s_temp)*A(n) ), delay(n) );
% end;
% clear s_temp;
% r = [ zeros( round( length(s)/2), 1 ) ; r ];
% figure(1); plot( (0:(length(r)-1))/fs, r );
%% Received Signal Generation by Convolution - Recommended
h_temp = zeros( round(max(delay)+length(s)/2), 1 );
h_temp(delay+1) = A;
h_temp = [zeros( round( length(s)/2), 1 ) ; h_temp];
r = real( conv( hilbert(s), h_temp ) );
clear h_temp;
figure(1); plot( (0:(length(r)-1))/fs, r );
%% AWGN
E = norm(s) * max( abs(A) ) / length(s);    SNR = 20;
w = sqrt( E/10^(SNR/10) ) * randn( size(r) );
rec = r + w;
figure(2); plot( (0:(length(r)-1))/fs , rec );
%% Matched Filter
result = conv( rec, flipud(s) );
result = abs( hilbert(result) );
result_temp = result/100/norm(s)^2;
result_temp(result_temp < 0.05*max(result_temp)) = 0;
tao = find( diff(sign(diff(result_temp))) == -2 ) + 1;
atten = result_temp(tao);
tao = ( tao ) - min(tao);
figure(3); plot( (1:(length(result)))/fs , result );
%% 
match = zeros( length(R), length(RD) );
matched = zeros( size(match) );
hw = waitbar( 0, 'Please wait ...');
for m = 1:length(R)
    for n = 1:length(RD)
        A_temp = Arr.A(m,:,n)*100;   delay_temp = Arr.delay(m,:,n);
        delay_temp = delay_temp( abs(A_temp) > thres );    A_temp = A_temp( abs(A_temp) > thres );
        delay_temp = round( ( delay_temp - min(delay_temp) ) * fs );
        h_temp = zeros( round(max(delay_temp)+length(s)/2), 1 );
        h_temp(delay_temp+1) = conj( A_temp );
        temp = real( conv( hilbert(rec), flipud(h_temp) ) );
        match(m,n) = norm( temp, Inf ) / norm( temp, 1 );
        index = find( abs(h_temp) > 0 );    h_temp = abs( h_temp(index) ); index = index - 1;
%         f_index = find( ismember(index,tao) == 1 );
%         b_index = find( ismember(tao,index) == 1 );
        [f_index b_index] = matchsqn( index, tao, round(length(s)/10) );
        if (~isempty(f_index)) && (~isempty(b_index))
            matched(m,n) = sum( abs(h_temp(f_index) ./...
                (abs(index(f_index)-tao(b_index))+1).^2 ./ (abs(h_temp(f_index)-atten(b_index))+1)) );
        end;
        waitbar( (m+n) / (length(R)+length(RD)) );
    end;
end;
close(hw);
clear hw h_temp temp A_temp delay_temp;
clear index f_index b_index;
[r_hat d_hat] = find( match == max(max(match)) );
[red_hat ded_hat] = find( matched == max(max(matched)) );
match = 10 * log10( match/max(max(match)) );
matched = 10 * log10( matched/max(max(matched)) );
%%
figure(4); pcolor( R, RD, match.' );	colormap(jet);    colorbar;
xlabel('Range (km)');  ylabel('Depth (m)');
shading interp; axis ij;
caxis([-10 0]);
figure(5); pcolor( R, RD, matched.' );	colormap(hot);    colorbar;
xlabel('Range (km)');  ylabel('Depth (m)');
shading interp; axis ij;
caxis([-20 0]);