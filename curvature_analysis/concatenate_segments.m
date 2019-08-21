prompt = {'How many segments to concatenate?', 'Total # of frames'};
dlgtitle = '# of segments and frames';
dims = [1 35];
definput = {'2', '1800'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
numseg = str2double(answer{1}); 
numframes = str2double(answer{2});

if numseg < numframes
    
    chunk_curv = cell(numseg,1);
    chunk_dorsal = cell(numseg,1);
    chunk_dorsal_r = cell(numseg,1);
    chunk_ventral = cell(numseg,1);
    chunk_ventral_r = cell(numseg,1);
    chunk_seg = zeros(numseg,2);

    for i = 1:numseg

        filename = uigetfile('.mat', ['Select segment' num2str(i)]);
        load(filename);
        deflength = size(curvdatafiltered,2);
        
        prompt = {['Start of segment ' num2str(i)], ['Length of segment ' num2str(i)],...
            ['Segment ' num2str(i) ' begins from:']};
        dlgtitle = 'Range of segment';
        dims = [1 35];
        definput = {'1', num2str(deflength), ''};
        answer = inputdlg(prompt,dlgtitle,dims,definput);    
        istart = str2double(answer{1});
        ilength = str2double(answer{2});
        framestart = str2double(answer{3});
        actualstart = framestart+istart-1;
        actualend = framestart+ilength-1;
        
        if istart>numframes || ilength>numframes || istart>ilength || actualstart>numframes
            fprintf('Inputs exceed boundary. \n');
        else
            chunk_curv{i,1} = curvdatafiltered(:, istart:ilength);
            chunk_dorsal{i,1} = dorsal_smd(:, istart:ilength);
            chunk_dorsal_r{i,1} = dorsal_smd_r(:, istart:ilength);
            chunk_ventral{i,1} = ventral_smd(:, istart:ilength);
            chunk_ventral_r{i,1} = ventral_smd_r(:, istart:ilength); 
            chunk_seg(i,1) = actualstart;
            chunk_seg(i,2) = actualend;
        end

    end
    
    total_curv = NaN(size(curvdatafiltered,1), numframes);
    total_dorsal = NaN(size(dorsal_smd,1), numframes);
    total_dorsal_r = NaN(size(dorsal_smd_r,1), numframes);
    total_ventral = NaN(size(ventral_smd,1), numframes);
    total_ventral_r = NaN(size(ventral_smd_r,1), numframes);
    
    for j = 1:numseg
        total_curv(:, chunk_seg(j,1):chunk_seg(j,2)) = chunk_curv{j,1};
        total_dorsal(:, chunk_seg(j,1):chunk_seg(j,2)) = chunk_dorsal{j,1};
        total_dorsal_r(:, chunk_seg(j,1):chunk_seg(j,2)) = chunk_dorsal_r{j,1};
        total_ventral(:, chunk_seg(j,1):chunk_seg(j,2)) = chunk_ventral{j,1};
        total_ventral_r(:, chunk_seg(j,1):chunk_seg(j,2)) = chunk_ventral_r{j,1};        
    end
    fprintf('concatenation is completed. \n');
    
    figure; imagesc(total_curv); title('After concatenation');
    fname_reg = regexp(filename, '\w*(?=[_])', 'match'); fname = fname_reg{1,1};
    save([fname '_concat.mat'], 'total_curv', 'total_dorsal', 'total_dorsal_r', ...
        'total_ventral', 'total_ventral_r');
    fprintf('concatenated data are saved. \n');
    
else
    
    fprintf('# of segments cannot exceed # of frames. \n');
    
end