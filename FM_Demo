% Frequency Modulation Demo
clear all;  close all;  clc;
%%
fs = 50e3; % Sample Frequency 采样频率
T = 10e-3; % Signal Duration 信号时长
t = (0:1/fs:T).';
L = length( t );
B = 4e3; % Bandwidth 带宽
k = B; % Modulation Index 调制指数
fc = 6e3; % starting frequency 起始频率
wav = [ 0;triang( L-2 );0 ]; % modulation waveform 待调制波形
wav = k * cumsum( wav / fs ); % Integral 积分
s = cos( 2*pi*(fc*t + wav) ); % transmitted signal(Modulated waveform) 最终发射信号(已调波)
%%
figure(1);  plot( t, s );
figure(2);  spectrogram( s, 64, 60, 64, fs );   colorbar;
figure(3);  plot( abs(fft(s)) );
