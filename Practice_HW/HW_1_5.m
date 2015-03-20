clear all; close all;clc;
% Author: jakejiangjn
% Date 2015-3-18
% All Rights Reserved
%%
fc = 10e3;
fs = 50e3;
int = 1e-2; % bit rate = 1/int
T = 11*int;
SNR = 0;
seq = [+1 +1 +1 -1 -1 -1 +1 -1 -1 +1 -1];
temp = (0:1/fs:int-1/fs).';L_s = numel(temp);
carrier_I = cos( 2*pi*fc*temp );
carrier_Q = sin( 2*pi*fc*temp );
x = carrier_I*seq;  L_x = numel(x);
x = reshape( x, 1, L_x );
temp = exp( 1j*unwrap(rand(1,L_x)-0.5)*pi ); % Phase Disturbance; a -2*pi~2*pi disturbance will highly cause severe Phase Ambiguity
x = real( 10^(SNR/20)*x.*temp + [1,1j]*randn(2,L_x)/sqrt(2) );
temp = (0:L_x-1)/L_x;
figure(1);plot( temp*T, x );title('Received Signals');xlabel('Time axis(s)');ylabel('Amplitude');
figure(2);plot( temp*fs, abs(fft(x)) );xlim([0 fs/2]);
title('Spectral Density of Received Signals');xlabel('Frequency axis(Hz)');ylabel('Magnitude');
%% IQ DeModulation
temp = 0:1/fs:T-1/fs;
carrier_I = 2*cos( 2*pi*fc*temp );
carrier_Q = 2*sin( 2*pi*fc*temp );
x_iq = x.*carrier_I - x.*1j.*carrier_Q; % Notice!
b = fir1(48,1e3*2/fs);
x_iq = filtfilt( b, 1, x_iq );
result_iq = reshape( real(x_iq), L_s, 11 );
result_iq = sign(sum(sign( result_iq )));
%% Hilbert DeModulation
x_h = hilbert(x);
x_h = x_h .*exp(-2j*pi*fc*temp);
x_h = filtfilt( b, 1, x_h );
result_h = reshape( real(x_h), L_s,11 );
result_h = sign(sum(sign( result_h )));
%% Figure Plot;
figure(3)
plot( temp,real(x_iq),'r', temp,real(x_h),'b' );legend( 'IQ DeModulation', 'Hilbert DeModulation' );
title('DeModulated Signals');xlabel('Time axis(s)');ylabel('Amplitude');
ylim([-1.5, 1.5]);
figure(4)
plot( temp,angle(x_iq)/pi,'r', temp,angle(x_h)/pi,'b' );legend( 'IQ DeModulation', 'Hilbert DeModulation' );
title('DeModulated Signals');xlabel('Time axis(s)');ylabel('Phase(\pi)');
figure(5);
subplot(3,1,1);stem(seq);title('Original Bits');
subplot(3,1,2);stem(result_iq);title('IQ DeMod Bits');
subplot(3,1,3);stem(result_h);title('Hilbert DeMod Bits');
