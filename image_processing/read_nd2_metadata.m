function timestamp = read_nd2_metadata(file_name)
    
    fid = fopen(file_name);
    
    while 1
        
        tline = fgetl(fid);
        
        if ~ischar(tline) 
            break, 
        end

        p = strfind(tline,'SizeT');
        p2 = strfind(tline,'SizeZ');

        if ~isempty(p)
            SizeT = sscanf(tline(p+5:end),'%f');
            timestamp = zeros(SizeT,1);
            j = 1;
        end

        if ~isempty(p2)
            SizeZ = sscanf(tline(p2+5:end),'%f');
        end

        k = strfind(tline,'timestamp');
        
        if ~isempty(k)
            c = sscanf(tline(k+10:end),'%f');
            timestamp(j) = c(2);
            j = j+1;
        end
 
    end
    
    timestamp = sort(timestamp);

fclose(fid);