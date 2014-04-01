function [AR_vec En AIC] = myRealAR( data )
%An AR model estimator - can only handle real input data
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
%    based on The Delsarte-Genin Algorithm
%   Reference: Stoica P, Moses R L. Spectral analysis of signals[M]. Pearson/Prentice Hall, 2005.
%               Section 3.5.2
if size( find( imag(data) ~= 0 ) ~= 0 )
    error( 'Input data must be a real vector' );
end;
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
D_0 = 1;	b_0 = R(1);	y_0 = R(2);
D_1 = [1;1];	b_1 = R(1) + R(2);	y_1 = R(2) + R(3);
ar_1 = [1;-R(2)/R(1)];  En_1 = R(1) - R(2)^2/R(1);
n = 1;AIC_1 = log( En_1 ) + 2*n/L;
range = min( max( L-1, L/10 ), 80 ); % Searching range of the order 
flag = (1);
% Iteration
while( flag )
    if n > range
        flag = (0); break;
    end;
    n = n+1;
    En_0 = En_1;    ar_0 = ar_1;    AIC_0 = AIC_1; % preserve the lastest result
    a = ( b_1 - y_1 ) / ( b_0 - y_0 );
    b_2 = 2*b_1 - a*b_0;
    D_2 = [D_1;0] + [0;D_1] - a*[0;D_0;0];
    y_2 = R(2:(n+2)) * D_2;
    En_1 = b_2 - b_2 /b_1 * y_1;
    ar_1 = D_2 - b_2/b_1*[0;D_1];
    AIC_1 = log( En_1 ) + 2*n/L;
    if AIC_0 < AIC_1
        AR_vec = ar_0;  En = En_0;  AIC = AIC_0;
        flag = (0); break;
    end;
    % Update
    y_0 = y_1;  y_1 = y_2;
    b_0 = b_1;  b_1 = b_2;
    D_0 = D_1;  D_1 = D_2;
end;
% AIC  = L*log(En) + 2*m;
% AICc = L*log(En) + (2*m * L)./(L-m-1);
% GIC  = L*log(En) +  m * v;
% BIC  = L*log(En) +  m * log(L);