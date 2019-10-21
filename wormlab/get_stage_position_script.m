[filename, pathname]  = uigetfile({'*.txt'});

fname = [pathname filename];

filetext = fileread(fname);
expr = '[^\n]*Image:[^\n]*';
matches = regexp(filetext,expr,'match');

sz = size(matches);

coordinates = zeros(sz(1,2),3); % variable for holding coordinates

for i = 1:sz(1,2)    
    Str = char(matches(1,i));    
    Str(strfind(Str, '=')) = [];   
    checkerx = contains(Str, 'xpos'); checkery = contains(Str, 'ypos'); checkerz = contains(Str, 'zpos');
    if checkerx||checkery||checkerz
        Keyx = 'xpos';
        Keyy = 'ypos';
        Keyz = 'zpos';
        Indexx = strfind(Str, Keyx);
        Indexy = strfind(Str, Keyy);
        Indexz = strfind(Str, Keyz);
        coordinates(i,1) = sscanf(Str(Indexx(1) + length(Keyx):end), '%g', 1);
        coordinates(i,2) = sscanf(Str(Indexy(1) + length(Keyy):end), '%g', 1);
        coordinates(i,3) = sscanf(Str(Indexz(1) + length(Keyz):end), '%g', 1);
    else
        coordinates(i,:) = NaN;       
    end    
end

% stage coordinates stored in coordinates
