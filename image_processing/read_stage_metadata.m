function [date,metadata]=read_stage_metadata(file_name)

fid=fopen(file_name);
i=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if i==0
        date=datevec(sscanf(tline,'%s'),'dddd,mmmdd,HH:MM:SS.FFF');
    else
        c=sscanf(tline,'%f');
        metadata(i,:)=c';
    end
    i=i+1;
end

metadata(:,1)=metadata(:,1)-metadata(1,1);
fclose(fid);



end

