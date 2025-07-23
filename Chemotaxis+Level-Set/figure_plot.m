%% Demo data (replace with your own)
% x = linspace(0,2*pi,200);
% y = sin(x);
%%

P = polyfit((5:1000)*0.1,log(medVals(5:1000)),1);
yfit = polyval(P,(5:1000)*0.1);
% hold on;
% 
%% Create figure and plot
figure('Color','w','Units','centimeters','Position',[4 4 14 10]);  % white background, size in cm
semilogy((1:1000)*0.1,medVals,'LineWidth',2,'Color',[0 0.45 0.74]);                  % thicker line, scientific blue
hold on
plot((5:1000)*0.1,exp(yfit),'r-.','LineWidth',2,'Color',[0.85 0.33 0.10]);
% plot((1:1000)*0.1,1.55891*ones(1000),'LineWidth',2,'Color',[0.85 0.33 0.10]);  
% 
% plot(x,cos(x),'--','LineWidth',1.2,'Color',[0.85 0.33 0.10]);     % example second curve
hold off

%% Axes cosmetics
ax = gca;
ax.FontName   = 'Helvetica';
ax.FontSize   = 9;          % 8–10 pt is typical in figures
ax.LineWidth  = 0.8;        % axes lines
ax.TickDir    = 'out';      % ticks pointing outward
ax.TickLength = [0.015 0.015];
ax.Box        = 'off';      % no box around plot
ax.XMinorTick = 'on';       % minor ticks help scientific reading
ax.YMinorTick = 'on';
grid(ax,'on');              % fine grid
ax.GridAlpha  = 0.15;       % lighter grid lines
ax.GridLineStyle = ':';     % dotted grid

%% Labels & title
xlabel('Time (s)','FontWeight','bold','Interpreter','latex');
ylabel('Mean Spectral Amplitude','FontWeight','bold','Interpreter','latex');
title('Mean Spectral Amplitude of Dominate Frequency','FontWeight','bold','FontSize',10);

%% Legend (optional)
legend({'Simulated','Predicted'},'Location','southwest',...
       'Interpreter','latex','FontSize',8,'Box','off');

%% Export-ready settings
set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 14 10]);   % same size on export
print(gcf,'-dpng','-r600','figure_scientific.png');                % 600 dpi PNG

