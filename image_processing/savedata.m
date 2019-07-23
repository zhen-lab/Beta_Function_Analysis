if length(questdlg('Save this data? '))==3
    
    clear data;
    clear imagelist;
    clear img_stack;
    [fn, savepathname]= uiputfile('*.mat', 'choose file to save', strcat(fname, '_',num2str(istart),'-',num2str(iend),'.mat'));
    if length(fn) > 1
        fnamemat = strcat(savepathname,fn);
        save(fnamemat);
    end
    
end