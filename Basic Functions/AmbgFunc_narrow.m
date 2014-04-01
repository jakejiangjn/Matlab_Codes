function [Ambg,delay,doppler] = AmbgFunc_narrow( x1, x2, fs )
%Narrow-band Signal Ambiguity Function Calculation
% Two row signal vectors are expected as well as their sample frequency fs.
% This functions returns the  Ambiguity Function's absolute value - Ambg
%                             Delay axis vector (timestamp) (s)   - delay
%                             Doppler axis value (Hz)             - doppler
% Delay interval is 1/fs, and Doppler interval is fs/length(x).
% input signals are Hilbert Transfomed for distortless.
%   Usage Demo:
% % [Ambg,delay,doppler] = AmbgFunc_narrow( s, fs );
% % Ambg = 20*log10( Ambg/max(max(Ambg)) ); 
% % surf( delay, doppler, Ambg);%contour
% % colorbar; caxis([-30 0]); % Setting diplay color range
% % xlabel('\tau');  ylabel('f_{D}');

% Checking
[m1,n1] = size( x1 );
if nargin == 2
    if m1 > n1
        x1 = x1.';
    end;
    [m1,n1] = size( x1 );
    if m1 ~= 1
        error( 'A ROW VECTOR is expected!' );
    end;
    fs = x2;    x2 = x1;
elseif nargin == 3
    [m2,n2] = size( x2 );
    if (m1 ~= 1 || m2 ~= 1) || (n1 ~= n2)
        error( 'Two equal lengthed ROW VECTORs are expected!' );
    end;
else
    error('Wrong Input Arguments');
end;
% Axis Setting
L = n1;
t = (0:L-1)/fs;
delay = (1-L):(L-1);
doppler = -fs/L*20:0.2*fs/L:fs/L*20;
% Ambiguity Function Calculation
m = length( doppler );  n = length( delay );
Ambg = zeros( m, n );
x1 = hilbert( x1 );
x2 = hilbert( x2 ); x2 = conj( x2 ); x2 = fliplr(x2);
for k = 1:m
    Ambg(k,:) = conv( x1.*exp(2j*pi*doppler(k)*t), x2 );
    Ambg(k,:) = abs( Ambg(k,:) );
end;