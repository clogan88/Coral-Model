function coralCoverFigure(fullDir, suffix, XData1, YData1, YData2, YData3, YData4, X1, YMatrix1)
%CREATEFIGURE(YDATA1, XDATA1, YDATA2, YDATA3, YDATA4, X1, YMATRIX1)
%  YDATA1:  patch ydata
%  XDATA1:  patch xdata
%  YDATA2:  patch ydata
%  YDATA3:  patch ydata
%  YDATA4:  patch ydata
%  X1:  vector of x data
%  YMATRIX1:  matrix of y data

%  Auto-generated by MATLAB on 03-Feb-2017 15:39:08

% Create figure
figure1 = figure('Name','Global Coral Cover');

% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

transparency = 0.2;
% Create patch for massive, 10-90
patch('Parent',axes1,'DisplayName','10th - 90th percentile','YData',YData1,...
    'XData',XData1,...
    'FaceAlpha',transparency,...
    'LineStyle','none',...
    'FaceColor',[1 0.8 0.8]);

% Create patch for branching, 10-90
patch('Parent',axes1,'DisplayName','10th - 90th percentile','YData',YData2,...
    'XData',XData1,...
    'FaceAlpha',transparency,...
    'LineStyle','none',...
    'FaceColor',[0.75 0.75 1]);

% Create patch for massive, 25-75
patch('Parent',axes1,'DisplayName','25th - 75th percentile','YData',YData3,...
    'XData',XData1,...
    'FaceAlpha',transparency,...
    'LineStyle','none',...
    'FaceColor',[1 0.5 0.5]);

% Create patch for massive, 25-75
patch('Parent',axes1,'DisplayName','25th - 75th percentile','YData',YData4,...
    'XData',XData1,...
    'FaceAlpha',transparency,...
    'LineStyle','none',...
    'FaceColor',[0.45 0.45 1]);

% Create multiple lines using matrix input to plot
plot1 = plot(X1,YMatrix1);
set(plot1(1),'DisplayName','Massive coral','Color',[1 0 0]);
set(plot1(2),'DisplayName','Branching coral','Color',[0 0 1]);

% Create xlabel
xlabel({'Year'});

% Create title
title('Global Coral Cover','FontWeight','bold');

% Create ylabel
ylabel('Percent of Carrying Capacity, K');

% Uncomment the following line to preserve the X-limits of the axes
% xlim(axes1,[1860 2100]);
box(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',14);
% Create legend
legend1 = legend(axes1,'show', 'Location', 'west');
%set(legend1,...
%    'Position',[0.392272766767302 0.35970390099993 0.134089389747584 0.207646171013633],...
%    'FontSize',13);
set(legend1,'FontSize',13);
print('-dpdf', '-r200', strcat(fullDir, 'GlobalCoralCover', suffix, '.pdf'));
if verLessThan('matlab', '8.2')
    saveas(gcf, strcat(fullDir, 'GlobalCoralCover', suffix), 'fig');
else
    savefig(strcat(fullDir, 'GlobalCoralCover', suffix, '.fig'));
end

