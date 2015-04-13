function [Ambg,delay,doppler] = AmbgFunc_WideBand( x, s, fs, c, v )
%Wide-band Signal Ambiguity Function Calculation
% Two row signal vectors are expected as well as their sample frequency fs.
%Inputs:
%		x : received signal ( A long vector is expected )
%		s : orginal signal	( a short vector is expected )
%		fs: sample frequency ( Hertz )
%		c : propagation speed ( m/s )
%		v : velocity search interval ( a vector is expected. Unit: m/s )
% This functions returns the  Ambiguity Function's absolute value - Ambg	( size(Ambg) = [length(doppler), length(delay)])
%                             Delay axis vector (timestamp) (s)   - delay 	( X-axis, row of Ambg at one doppler)
%                             Doppler axis value (Hz)             - doppler ( Y-axis, column of Ambg at one delay)
%		Ambg	( size(Ambg) = [length(doppler), length(delay)])
%       delay 	( X-axis, row of Ambg at one doppler)
%       doppler ( Y-axis, column of Ambg at one delay)
% Delay interval is 1/fs, and Doppler is determined by v.
% input signals are Hilbert Transfomed for distortless.
%   Usage Demo:
% % [Ambg,delay,doppler] = AmbgFunc_WideBand( x, s, fs, c, v );
% % Ambg = 20*log10( Ambg/max(max(Ambg)) ); 
% % surf( delay, doppler, Ambg); shading interp; % contour
% % colorbar; caxis([-30 0]); % Setting diplay color range
% % xlabel('Delay \tau (sec)');  ylabel('f_{Doppler} (Hz)');
if nargin < 5
	error('Wrong Input Arguments! Usage: [Ambg,delay,doppler] = AmbgFunc_WideBand( x, s, fs, c, v ); For more, see HELP.\n\r');
end;
x = check(x);	s = check(s);	v = check(v);
%% Wideband Ambiguity Generation
L_x = length(x);	L_s = length(s);	L_v = length(v);
% N_c = length( resample(s, c, min(v)*2+c ) );
N_c = length( resample(s, c-min(v), min(v)+c ) );
N = N_c + L_x - 1;
Ambg = zeros( L_v, N );
for n = 1:L_v
%     temp = resample( s,c,c+2*v(n) ) * sqrt( 1 + 2*v(n)/c );
    temp = resample( s,c-v(n),c+v(n) ); temp = temp / norm(temp);
    temp = [ temp zeros( 1, N_c-length(temp) ) ];
    temp = abs( hilbert(conv( x, fliplr(temp) )) );
    Ambg(n,:) = temp; % revised
    % l = length( temp ); % original
    % A(n,N-l+1:end) = temp; % original
end;
delay = ( (1:N)-N_c ) / fs;
doppler = v/c;

function result = check( input )
if size( input, 1 )~=1 && size( input, 2 )~=1
    error('Vector is expected!');
end;
result = input(:).';
