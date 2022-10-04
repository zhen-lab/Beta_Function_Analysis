function [xyz_col] = read_motorstage_position_behavior(data)

%     data = bfopen(fname);
    metadata = data{1, 2}; % Metadata is stored in the 2nd variable
    dim = 3; % Dimension for the coordinates
    
%     metadataKeys = metadata.keySet().iterator();
%     for i=1:metadata.size()
%       key = metadataKeys.nextElement();
%       if i==10 % Coordinates information is all stored in the 10th element
% %       if i==3 % Sometimes stored in the 3rd element???????????? WTF
%         value = metadata.get(key);
%       end
%     end
    maxsize = 0; 
    metadataKeys = metadata.keySet().iterator();
    for i=1:metadata.size()
        key = metadataKeys.nextElement();
        value = metadata.get(key);
        valuesize = size(value,2);
        if valuesize>maxsize
            maxsize = valuesize;
            value_target = value;
            number_target = i;
        end
    end
    
    xyz = regexp(value_target, '[-]?\d*[.]\d*', 'match'); % Match the regular expression of numeric values
    fprintf(['Stage position info is stored in #: ' num2str(number_target) '\n']);
    xyz_col = reshape(xyz, dim, [])';
    
end