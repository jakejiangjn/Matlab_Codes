function [ Modes ] = read_modes_rd( filename, NProf )
% Read the modes produced by KRAKEN
% usage:
%    [ Modes ] = read_modes( filename, NProf )
% filename should include the extension
% iProf is the profile number
% modes is an optional vector of mode indices

% mbp, May 2001
% revised by jakejiangjn@zju.edu.cn

% identify the file type

switch filename(1:6)
    case 'MODFIL'
        FileType = 'mod';
    case 'MOAFIL'
        FileType = 'moa';
    otherwise
        endchar = length( filename );
        if ( endchar >= 4 )
            FileType = lower( filename( endchar-2 : endchar ) );
        end
end

% read the modal data
if strcmp(FileType, 'mod') % binary format
    
    fid   = fopen( filename, 'r' );
    if ( fid == -1 )
        errordlg( 'Mode file does not exist', 'read_modes_bin' )
        error(    'read_modes_bin: Mode file does not exist' )
    end
    iRecProfile = 1;   % (first time only)
    lrecl = 4 * fread( fid, 1, 'long' );
    for n = 1:NProf
        [mode,iRecProfile] = read_modes_bin_rd( fid, iRecProfile, lrecl );
        % calculate wavenumbers in halfspaces
        if ( mode.Top.bc == 'A' )   % top
            mode.Top.k2     = ( 2 * pi * mode.freq / mode.Top.cp )^2;
            gamma2           = mode.k .^ 2 - mode.Top.k2;
            mode.Top.gamma  = PekerisRoot( gamma2 );
            mode.Top.phi    = mode.phi( 1, : );   % mode value at halfspace
        end
        
        if ( mode.Bot.bc == 'A' )   % bottom
            mode.Bot.k2    = ( 2 * pi * mode.freq / mode.Bot.cp )^2;
            gamma2          = mode.k .^ 2 - mode.Bot.k2;
            mode.Bot.gamma = PekerisRoot( gamma2 );
            mode.Bot.phi   = mode.phi( end, : );   % mode value at halfspace
        end
        Modes(n) = mode;
    end
else
    warndlg( 'Unrecognized file extension', 'Warning' )
end
fclose(fid);
%%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [ Modes,iRecProfile ] = read_modes_bin_rd( fid, iRecProfile, lrecl )
% Read the modes produced by KRAKEN
% usage:
%    [ Modes ] = read_modes_bin( filename, modes )
% filename is without the extension, which is assumed to be '.moA'
% modes is an optional vector of mode indices
%
% This version can keep reading modesets from the same modefile.
% On the first call, it opens the mode file and sets the record pointer to
% the beginning of the file.
% From then on it reads sequentially.

% derived from readKRAKEN.m    Feb 12, 1996 Aaron Thode
%
% Modes.M          number of modes
% Modes.k          wavenumbers
% Modes.z          sample depths for modes
% Modes.phi        modes
%
% Modes.Top.bc
% Modes.Top.cp
% Modes.Top.cs
% Modes.Top.rho
% Modes.Top.depth
%
% Modes.Bot.bc
% Modes.Bot.cp
% Modes.Bot.cs
% Modes.Bot.rho
% Modes.Bot.depth
%
% Modes.N          Number of depth points in each medium
% Modes.Mater      Material type of each medium (acoustic or elastic)
% Modes.Nmedia     Number of media
% Modes.depth      depths of interfaces
% Modes.rho        densities in each medium

rec = iRecProfile - 1;
fseek( fid, rec * lrecl + 4, -1 );

Modes.title  = fread( fid, 80, '*char' )';
disp( Modes.title )
Modes.freq   = fread( fid,  1, 'float' );
Modes.Nmedia = fread( fid,  1, 'long'  );
Ntot         = fread( fid,  1, 'long'  );
NMat         = fread( fid,  1, 'long'  );

if Ntot < 0, return; end

% N and Mater
rec   = iRecProfile;
fseek( fid, rec * lrecl, -1 );
for Medium = 1 : Modes.Nmedia
    Modes.N(     Medium ) = fread( fid, 1, 'long' );
    Modes.Mater( Medium, : ) = fread( fid, 8, '*char' )';
end

% Top
rec = iRecProfile + 1;
fseek( fid, rec * lrecl, -1 );
Modes.Top.bc    = fread( fid, 1, '*char' );
cp              = fread( fid, [ 2, 1 ], 'float' );
Modes.Top.cp    = complex( cp( 1 ), cp( 2 ) );
cs              = fread( fid, [ 2, 1 ], 'float' );
Modes.Top.cs    = complex( cs( 1 ), cs( 2 ) );
Modes.Top.rho   = fread( fid, 1, 'float' );
Modes.Top.depth = fread( fid, 1, 'float' );

% Bottom
Modes.Bot.bc    = char( fread( fid, 1, 'char' )' );
cp              = fread( fid, [ 2, 1], 'float' );
Modes.Bot.cp    = complex( cp( 1 ), cp( 2 ) );
cs              = fread( fid, [ 2, 1], 'float' );
Modes.Bot.cs    = complex( cs( 1 ), cs( 2 ) );
Modes.Bot.rho   = fread( fid, 1, 'float' );
Modes.Bot.depth = fread( fid, 1, 'float' );

% depth and rho
rec = iRecProfile + 2;
fseek( fid, rec * lrecl, -1 );
bulk        = fread( fid, [ 2, Modes.Nmedia ], 'float' );
Modes.depth = bulk( 1, : );
Modes.rho   = bulk( 2, : );

% m
rec = iRecProfile + 3;
fseek( fid, rec * lrecl, -1 );
Modes.M = fread( fid, 1, 'long' );
% Lrecl   = fread( fid, 1, 'long' );

% z
rec = iRecProfile + 4;
fseek( fid, rec * lrecl, -1 );
Modes.z = fread( fid, Ntot, 'float' );

% read in the modes
modes = 1 : Modes.M;    % read all modes

rec = iRecProfile + 5;
fseek( fid, rec * lrecl, -1 );

% if there are modes, read them
if ( Modes.M == 0 )
    Modes.phi = [];
    Modes.k   = [];
else
    Modes.phi = zeros( NMat, length( modes ), 'single' );   %number of modes
    
    for ii = 1: length( modes )
        rec = iRecProfile + 4 + modes( ii );
        fseek( fid, rec * lrecl, -1 );
        phi = single( fread( fid, [ 2, NMat ], 'single' )' ); %Data is read columwise
        Modes.phi( :, ii ) = phi( :, 1 ) + 1i * phi( :, 2 );
    end
    
    % following is faster, but only valid if all the record sizes are the same
    % phitmp    = fread( fid, [ 2, Ntot * Modes.M ], 'float' ); %Data is read columwise
    % phitmp    = phitmp( 1, : ) + 1i * phitmp( 2, : );
    % Modes.phi = reshape( phitmp, Ntot, Modes.M );
    %
    % read in the wavenumbers
    
    % Ifirst = 1;
    % cktot = [];
    
    %for I = 1 : ( 1 + ( 2 * m - 1 ) / Lrecl ),
    %   rec = 5 + m + I;
    %   fseek( fid, rec * lrecl, -1 );
    %   Ilast = min( [ m Ifirst + Lrecl / 2 - 1 ] );
    %   ck = fread( fid, [ 2, Ilast - Ifirst + 1 ], 'float' )';
    %   cktot = [ cktot; ck ];
    %   Ifirst = Ilast + 1;
    %end
    %ck = cktot( modes, 1 ) + i * cktot( modes, 2 );
    
    rec = iRecProfile + 5 + Modes.M;
    fseek( fid, rec * lrecl, -1 );
    
    k       = zeros( 2, Modes.M, 'single' );   % pre-allocate single
    Modes.k = zeros( 1, Modes.M, 'single' );
    
    k    = single( fread( fid, [ 2, Modes.M ], 'float' ) );
    Modes.k = ( k( 1, : ) + 1i * k( 2, : ) ).';   % column vector
end

iRecProfile = iRecProfile + 6 + Modes.M + 1 + floor( ( 2 * Modes.M - 1 ) / lrecl );   % advance to next profile