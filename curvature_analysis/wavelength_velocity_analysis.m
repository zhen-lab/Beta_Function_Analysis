title('click start and end frames to analyze phase velocity and angle of attacks');

figure; imagesc(curvdatafiltered(:,:)); colormap(cmap); colorbar; caxis([-10 10]); hold on;

t1 = ginput(1);

plot( [1 numcurvpts],[t1(2) t1(2)], '-w');
t2 = ginput(1);

plot( [1 numcurvpts],[t2(2) t2(2)], '-w');

t_start=ceil(t1(2));
t_end=ceil(t2(2));

colormap(jet);
c2 = curvdatafiltered(t_start:t_end,:) > 0;
figure; imagesc(c2(:,:));  hold on;
title('filtered, binary');

maskhead = 0.1;
masktail = 0.1;
minimum_fraction_for_fit = 0.8;

c3 = edge(single(c2),'sobel', 0);

c3(:,round((1-masktail)*numcurvpts):end) = 0;
c3(:,1:round(maskhead*numcurvpts)) = 0;
        
[c4 numlab] = bwlabel(c3);

numcycles2 = 0;

%clear slopedata timedata slopedatatmp timedatatmp okdatatmp curvsigndatatmp curvsigndata;

okdatatmp = zeros(numlab, 1);

normrthresh = 220;

%   draw fit limits

%     
hold on;

for n=1:numlab
    c5 = (c4 == n);
    [y x] = find(c5);

    yshift = 3;
    yshifted = ceil(1+0.5*(1+sign(y-yshift)) .* (y-yshift-1));
    curvshift = zeros(size(x));
            
    for jj=1:length(x)
        curvshift(jj) = curvdata(yshifted(jj), x(jj));
    end
    %         
    %         disp(n);
    %         disp('mean(c2shift)');
    %         disp(mean(curvshift));
    %         
            % exclude points near head and tail for fitting
    tmp = x;
    x=x(logical((tmp>=maskhead * numcurvpts) .* (tmp<=(1-masktail) * numcurvpts)));
    y=y(logical((tmp>=maskhead * numcurvpts) .* (tmp<=(1-masktail) * numcurvpts)));

    if max(x) - min(x) >=  (1-maskhead-masktail)*minimum_fraction_for_fit*numcurvpts
        
        [p,S] = polyfit(x,y,1);
    %             plot(polyval(p,[1:numcurvpts]), 'r'); hold on;    
        %         disp([n , S.normr , S.normr < normrthresh]);
    %             S.normr

        if S.normr < normrthresh
            if mean(curvshift) > 0
                plotcol = '-g';
            else
                plotcol = '--g';
            end
            plot(polyval(p,[1:numcurvpts]), plotcol); hold on;        
            numcycles2 = numcycles2 + 1;
            slopedatatmp(n) = p(1);
            timedatatmp(n) = p(2);
            okdatatmp(n) = 1;
            negshift = (mean(curvshift) > 0);
            curvsigndatatmp(n) = negshift;
            xpos = 5;
            ypos = p(2)-1;
            if p(2)<1
                xpos = numcurvpts/4;
                ypos = 5;
            end
                text(xpos,ypos,num2str([numcycles2 p(1)]), 'Color', 'white'); hold on;
        end
    end
   
end

plot( [0.5+maskhead*numcurvpts 0.5+maskhead*numcurvpts],[1 numframes], ':k');
plot( [0.5+ (1-masktail)*numcurvpts 0.5+ (1-masktail)*numcurvpts],[1 numframes], ':k');
        
title(strcat(num2str(sum(curvsigndatatmp)),' positive, ', num2str(numcycles2-sum(curvsigndatatmp)),...
                                ' negative.  Click on bad fits, press return'));
badfits = ginput;

epsilon = 4;  % how close to fit you need to click
for j=1:size(badfits,1)
    for n=1:numlab
        if okdatatmp(n) 
            if abs(timedatatmp(n) + slopedatatmp(n)*badfits(j,1) - badfits(j,2))<epsilon
                %                 disp(strcat('matches #', num2str(n)));
                okdatatmp(n) = 0;
            end
        end
    end
end

numcycles2 = 0;
c4b = c4;
for n=1:numlab
    if okdatatmp(n)
        numcycles2 = numcycles2+1;
        slopedata(numcycles2) = slopedatatmp(n);
        timedata(numcycles2) = timedatatmp(n);
        curvsigndata(numcycles2)=curvsigndatatmp(n);
        c4b(c4b==n) = numcycles2;
    else
        c4b(c4b==n) = 0;
    end
end

imagesc(c2); hold on;

plot( [0.5+maskhead*numcurvpts 0.5+maskhead*numcurvpts],[1 numframes], ':k');
plot( [0.5+ (1-masktail)*numcurvpts 0.5+ (1-masktail)*numcurvpts],[1 numframes], ':k');

for n=1:numcycles2
    if curvsigndata(n)
        plotcol = '-g';
    else
        plotcol = '--g';
    end

        plot(polyval([slopedata(n) timedata(n)],[1:numcurvpts]), plotcol); hold on;    
        %             text(5,p(2)-1,num2str([n numcycles2 S.normr]), 'Color', 'white'); hold on;
        %             if p(2)<1 
        %                 text(5,2,num2str([n numcycles2 S.normr]), 'Color', 'white'); hold on;
        %             end
end

if mean(curvsigndata) ~= 0.5
    msgbox('Warning: unequal number of positive and negative fits','','error') 
end

title('Press return to continue');
pause;
imagesc(c2(:,:));
    
title('click two points separated in time by N cycles');
t1 = ginput(1);

plot( [1 numcurvpts],[t1(2) t1(2)], '-w');
t2 = ginput(1);

plot( [1 numcurvpts],[t2(2) t2(2)], '-w');
    %   draw fit limits

plot( [0.5+maskhead*numcurvpts 0.5+maskhead*numcurvpts],[1 numframes], ':k');
plot( [0.5+ (1-masktail)*numcurvpts 0.5+ (1-masktail)*numcurvpts],[1 numframes], ':k');

answer = inputdlg('Enter number of cycles');
period = abs(t1(2) - t2(2)) / str2num(answer{1});

title(strcat(answer{1},' cycles, ', num2str(numcycles2), ' fits'), 'Interpreter', 'None');


for n=1:numcycles2
    if curvsigndata(n)
        plotcol = '-g';
    else
        plotcol = '--g';
    end
    
    plot(polyval([slopedata(n) timedata(n)],[1:numcurvpts]), plotcol); hold on;    

        %             text(5,p(2)-1,num2str([n numcycles2 S.normr]), 'Color', 'white'); hold on;
        %             if p(2)<1 
        %                 text(5,2,num2str([n numcycles2 S.normr]), 'Color', 'white'); hold on;
        %             end
end


period_s = period/fps;
disp('mean period (sec)')
disp(period_s);


wavevelocity = mean(1./slopedata)/numcurvpts*fps;
disp('mean wave velocity (L/s)');
disp(wavevelocity);

wavelength = wavevelocity * period_s;
disp('wavelength = velocity*period (L)');
disp(wavelength);



segment1_curv=std(mean(curvdatafiltered(t_start:t_end,1:20),2));
segment2_curv=std(mean(curvdatafiltered(t_start:t_end,21:40),2));
segment3_curv=std(mean(curvdatafiltered(t_start:t_end,41:60),2));
segment4_curv=std(mean(curvdatafiltered(t_start:t_end,61:80),2));
segment5_curv=std(mean(curvdatafiltered(t_start:t_end,81:100),2));

disp('segment amplitude (1-5), from anterior to posterior');
disp([segment1_curv segment2_curv segment3_curv segment4_curv segment5_curv]);





segment1_dorsal_calcium=std(mean(dorsal_brightness_data_filtered(t_start:t_end,16:35),2));
segment2_dorsal_calcium=std(mean(dorsal_brightness_data_filtered(t_start:t_end,36:55),2));
segment3_dorsal_calcium=std(mean(dorsal_brightness_data_filtered(t_start:t_end,56:75),2));
segment4_dorsal_calcium=std(mean(dorsal_brightness_data_filtered(t_start:t_end,76:90),2));




segment1_ventral_calcium=std(mean(ventral_brightness_data_filtered(t_start:t_end,16:35),2));
segment2_ventral_calcium=std(mean(ventral_brightness_data_filtered(t_start:t_end,36:55),2));
segment3_ventral_calcium=std(mean(ventral_brightness_data_filtered(t_start:t_end,56:75),2));
segment4_ventral_calcium=std(mean(ventral_brightness_data_filtered(t_start:t_end,76:90),2));


disp('ventral calcium level (1-4), from anterior to posterior');
disp([segment1_ventral_calcium,segment2_ventral_calcium,segment3_ventral_calcium,segment4_ventral_calcium]);



disp('dorsal calcium level (1-4), from anterior to posterior');
disp([segment1_dorsal_calcium,segment2_dorsal_calcium,segment3_dorsal_calcium,segment4_dorsal_calcium]);














