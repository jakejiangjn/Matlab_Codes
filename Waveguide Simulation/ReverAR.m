clear all;close all;clc;
%%
load( 'Rever_coh.mat' ); % For 'rec_coh' and 'delay_coh'
load( 'Rever_inc.mat' ); % For 'rec_inc' and 'delay_inc'
L = length( rec_coh );
beg = find( rec_coh ~= 0, 1 );
data = rec_coh(beg:beg+999).';
%%
x = zscore( data ); % returns a centered, scaled version of the input data
% i.e. z = (x–mean(x))./std(x)
z0 = iddata( data, [], 1/50e3 ); % creates an iddata object
%确定模型阶数
range = 80; % the search range of the order of AR model
test_factor =zeros( range, 2 );
opt = arOptions('Approach', 'yw'); % AR options
tic
for p=1:80
    m = ar( z0, p, opt ); % Estimate parameters of AR model for scalar time series
    AIC = aic( m ); % Akaike Information Criterion for estimated model
    test_factor(p,:) = [p  AIC];
end
[AICvalue,AICloc] = min( test_factor(:,2) )
p_best = AICloc;
%拟合------
m = ar( z0, p_best, opt );
compare( z0, m, 1 );
toc
% ar_cof = m.a;
% save( 'AR_matlab.mat', 'ar_cof' );
%%
tic
% [ar_c En AIC] = myRealAR( data );
[ar_c En AIC] = myComplexAR( data );
toc
%%
load( 'AR_matlab.mat' );
%%
flt_result = ( filtfilt(ar_c,1,rec_coh.') ).'; % Flitering out
% flt_result = ( filtfilt(ar_cof,1,rec_coh.') ).'; % Flitering out
% flt_result = ( filter(ar_c(2:end),1,data.') ).'; % Estimated Waveform
% flt_result = ( filter(ar_cof(2:end),1,data.') ).';
plot(flt_result,'b');
% mm = m; mm.a = ar_c.'; mm.NoiseVariance = En;
% compare( z0, mm, 1 );
% hold on;plot(data,'r'); hold off;
%%
T = 10e-3; % Singal Pulse Duration
fs = 50e3;  freq = 6e3;
s = 0:(1/fs):(T-1/fs);  s = cos( 2*pi*freq*s) .* hanning( length(s) ).';
flt_result = ( filtfilt(ar_c,1,s.') ).';
figure(2);plot(flt_result,'r');
% hold on;plot(s,'b');hold off;