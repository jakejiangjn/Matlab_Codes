function mat2labn(data,filename)
%data:data number*element number;
%for example:filename='s*';

% data_for_write=data.';
% data=data_for_write(1,:);
% fid=fopen(filename,'wt');
%
% for time1=1:size(data_for_write,2)
%         fprintf(fid,'%f\t',data_for_write(:,time1));
%         fprintf(fid,'\n');
% end

fid=fopen(filename,'wt');

for time1=1:size(data,2)
    fprintf(fid,'%f\t',data(:,time1));
    fprintf(fid,'\n');
end

% data_for_write=data.';
% data=data_for_write(1,:);
% fid=fopen(filename,'wt');
%
% for time1=1:size(data_for_write,1)
%     for time2=1:(length(data)-1)
%         fprintf(fid,'%f\t',data_for_write(time1,time2));
%     end
%         fprintf(fid,'%f\n',data_for_write(time1,length(data)));
% end
fclose all;