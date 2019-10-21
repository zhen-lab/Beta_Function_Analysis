function plot_position_and_signal(filename)

alldata = xlsread(filename);

position = alldata(:, 36:37); 
position = position-repmat(position(1, :), size(alldata,1), 1);

ratio = alldata(:, 1);
ratio_smooth = smooth(ratio);
max_r = 10; min_r = -10;

cmap = jet(length(position));
lut = floor(ratio_smooth/(max_r)*(length(cmap)-1)/2);
N = length(position); 

hold on;
for i=1:N
    
   if ~isnan(ratio(i)) && (ratio_smooth(i)<=max_r) && (ratio_smooth(i)>=min_r)
    plot(position(i,1),position(i,2),'o','MarkerFaceColor',cmap(floor(abs(lut(i))+(sign(lut(i))>0)*(length(cmap)-1)/2)+1,:),'MarkerEdgeColor','none','MarkerSize',4);
   end
%    disp(cmap(abs(lut(i))+(sign(lut(i))>0)*(length(cmap)-1)/2+1,:));

end

plot(position(:,1), position(:,2), '-k', 'LineWidth', 0.5);
% hold on;

% colorbar;

% hold on; text(position(1,1),position(1,2),'start');
% 
% hold on; text(position(N,1),position(N,2),'end');

%return;

%figure;

%h=quiver(position(1:10:N-1,1),position(1:10:N-1,2),velocity(1:10:N-1,1),velocity(1:10:N-1,2));

%hkid=get(h,'Children');

%X=get(hkid(1),'Xdata');
%Y=get(hkid(1),'Ydata');

%lut_quiver=lut(1:10:N-1)+1;

%figure;

%for ii = 1:3:length(X)-1

%    headWidth = 200 * sqrt((X(ii+1)-X(ii)).^2 + (Y(ii+1)-Y(ii)).^2); % set the headWidth, function of length of arrow
%    ah = annotation('arrow',...
%        'Color', cmap(lut_quiver((ii+2)/3),:),...
%        'headStyle','cback1','HeadLength',5*10^4,'HeadWidth',headWidth);
%    set(ah,'parent',gca);
%    set(ah,'position',[X(ii) Y(ii) X(ii+1)-X(ii) Y(ii+1)-Y(ii)]);
%end

%for j=1:30:N-1
%    h=quiver(position(j,1),position(j,2),velocity(j,1),velocity(j,2),2,'Color',cmap(lut(j)+1,:),'linewidth',2);
%    adjust_quiver_arrowhead_size(h,5);
%    hold on;
%end

%colorbar; caxis([min_r max_r]);

%hold on; text(position(1,1),position(1,2),'start');

%hold on; text(position(N,1),position(N,2),'end');
    
end 


