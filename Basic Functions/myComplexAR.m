function [AR_vec En AIC] = myComplexAR( data )
%An AR model estimator - can handle complex input data
% When given a specific input time-domain data vector,
% this function returns an appropriate estimated AR coefficient vector
% (based on AIC criterion), its estimation error(Loss Function)
% and its AIC information value.
% Notes:
%   Akaike's Information Criterion(AIC) = 
%       log( En ) + 2*( m-1 ) / L; 
%   or MDL Criterion can be applied, MDL =
%       L*log( En ) + ( m-1 )*log( L )
%   L = length( data );m = length(ar_vec)
%
%    based on The Levinson-Durbin Algorithm
%   Reference: Stoica P, Moses R L. Spectral analysis of signals[M]. Pearson/Prentice Hall, 2005.
%               Section 3.5.1
[m,n] = size( data );
if m ~= 1 && n ~= 1
    error( 'Input data must be a vector' );
end;
if m > n
    data = data.';
end;
L = length( data );
R = xcorr( data ) / L;
R = R(L:end);
% Initialization
ar_1 = -R(2)/R(1);  En_1 = R(1) - R(2)*conj( R(2) )/R(1);   b_1 = ar_1;
n = 1;AIC_1 = log( En_1 ) + 2*n/L;
range = min( max( L-1, L/10 ), 80 ); % Searching range of the order 
flag = (1);
% Iteration
while( flag )
    if n > range
        flag = (0); break;
    end;
    n = n+1;
    % Preserve the lastest result
    En_0 = En_1;    ar_0 = ar_1;    AIC_0 = AIC_1;	b_0 = b_1;
    b_1 = -( R(n+1) + conj( fliplr(R(2:n)) )*ar_0 ) / En_0;
    En_1 = En_0 * ( 1 - b_1*b_1' );
    ar_1 = [ar_0;0] + b_1*[conj( flipud(ar_0) );1];
    AIC_1 = log( En_1 ) + 2*n/L;
    if AIC_0 < AIC_1
        AR_vec = [1;ar_0];  En = En_0;  AIC = AIC_0;
        flag = (0); break;
    end;
end;
% AIC  = L*log(En) + 2*m;
% AICc = L*log(En) + (2*m * L)./(L-m-1);
% GIC  = L*log(En) +  m * v;
% BIC  = L*log(En) +  m * log(L);