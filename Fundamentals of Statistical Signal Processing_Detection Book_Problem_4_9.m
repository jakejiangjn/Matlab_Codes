clear all;  close all;  clc;
%%
f_0 = 1/13;
f = (0:0.01:2);
cmd = ['r', 'b', 'k', 'g'];
figure(1);  hold on;
for n = 2 : 5
    N = 10^n;
    t = 0:(N-1);
    s = cos( 2*pi * f_0 * t );
	s = fliplr( s ); % impulse response of matched filter h[n] = s[N-1-n]
	
    %% using FFT
%     S = fft( s )/sqrt( N );
%     f = t / N * 2;
%     f(f>1) = f(f>1) - 2;
	
    %% DIY 
    [t_0,t_1] = meshgrid( t, f );
    S = s * exp( -1i*pi*(t_0').*(t_1') );
	
    %% Plot work
    plot( f, abs(S), cmd(n-1) );
end;
hold off;
%%
xlabel('Frequency(\omega/\pi)');  ylabel('Amplitude');
legend('N=1e2','N=1e3','N=1e4','N=1e5');
title('Spectrum Magnitude|H(f)| v.s. Sample Length(N)');
