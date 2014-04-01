%Read data from .fig file or .png file
% clear all;  close all;  clc;
%% Function Mode
function [x,y] = fig2data( filename )
%%
open( filename );
data = findall( gca, 'type', 'line' );
data_x = get( data, 'xdata' );  data_y = get( data, 'ydata' );
close all;
if iscell( data_x ) == 1
    x = data_x{1,1};    y = data_y{1,1};
else
    x = data_x;         y = data_y;
end;
clear data data_x data_y;