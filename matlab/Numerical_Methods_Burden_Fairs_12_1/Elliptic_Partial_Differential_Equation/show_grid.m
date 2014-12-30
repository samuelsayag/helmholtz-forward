mygrid = linspace(0.125,0.375,3);
x = repmat(mygrid,3,1);
y = repmat(flip(mygrid'),1,3);

create_grid_figure(x,y)

% % Create figure
% figure1 = figure;
% 
% % Create axes
% axes1 = axes('Parent',figure1,'YTick',[0 0.125 0.25 0.375 0.5],'YGrid','on',...
%     'XTickLabel',{'0','0.125','0.25','0.375','0.5'},...
%     'XTick',[0 0.125 0.25 0.375 0.5],...
%     'XGrid','on');
% 
% hold off
% 
% hold on
% 
% for i=1:3
%     plot(x(i,:),y(i,:),'.', 'MarkerSize', 12)
% end
% 
% hold off


