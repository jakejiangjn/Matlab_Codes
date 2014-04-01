function data = sra_read(filename, arr_num )
% Read the .dat file recorded by NI instruments.
% e.g.: filename = 'Binarydata.dat';
%       arr_num  = array number; (with default 16channels).
if nargin == 1
    arr_num = 16;
end;
fid = fopen( filename, 'r' );
r = fread( fid, 'double' );
sam = 1000; % sample to read;
r_temp = reshape( r, sam, length(r)/sam );
clear r;
r_temp = r_temp.';
for time1 = 1:arr_num
    temp = ( downsample(r_temp,arr_num,time1-1) ).';
    temp1(time1,:) = reshape(temp,1,[]);                                                                                                                                                  
end
r=temp1;
% figure;
% for time1=1:arr_num
%     plot((r(time1,:)-(time1-1)));    
%     hold on;
% end
data = r;
fclose( fid );
