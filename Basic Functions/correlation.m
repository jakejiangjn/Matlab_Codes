function [coefficient, index] = correlation( x, varargin )
%correlation Cross-correlation cCalculation.
%   [coefficient, index] = correlation( x, y )
%   , where x and y are length M vector and length N vector(M,N>1)
%   , returns the length M+N-1 cross-correlation sequence C.
%   [coefficient, index] = correlation( x, y, flag )  
%   flag only indicates whether an envelop is needed.
%
%
%   See also XCOV, CORRCOEF, CONV, CCONV, COV and XCORR2.

%   Author(s): Jake Jingning Jiang jakejiangjn
flag = 0; % coefficient becomes an envelop form
switch nargin,
    case 1,
        y = x;
    case 2,
        if isnumeric(varargin{1})
            if length(varargin{1}) == 1
                y = x; flag = 1;
            else
                y = varargin{1};
            end;
        else
            error(message('Unkown Input!'));
        end;
    case 3,
        if isnumeric(varargin{2})
            if length(varargin{1}) > 1
                y = varargin{1};flag = 1;
            else
                error(message('Unkown Input!'));
            end;
        else
            error(message('Unkown Input!'));
        end;
end;
x = check( x ); y = check( y );
y = fliplr( conj( y ) );
matched = conv( x, y );
coefficient = conv( abs(x).^2, ones( size(y) ) );
coefficient = matched ./ sqrt( y*y' * coefficient );
coefficient( isnan(coefficient) ) = 0;
if flag == 1
    coefficient = abs( hilbert(coefficient) );
end;
index = ( 1:length(x)+length(y)-1 ) - length(y);

function result = check( input )
result = input;
if size( input, 1 ) > size( input, 2 )
    result = result.';
end;
if min( size( input ) ) > 1
    result = result(1,:);
end;
