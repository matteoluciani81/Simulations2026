
function BL_Histogram(WW,h,xmax)

ntb=length(WW);
grey=[.8 .8 .8];                                                            % grey is for grid

try isempty(h); if isempty(h); h=0.25; end; catch; h=0.25; end              % bin width
try isempty(xmax);      % ------------------------------------------------- % upper bound
    if isempty(xmax); xmax=ML_round(max(abs(ML_minmax(WW))),1,2); end       %
catch                                                                       %
    xmax=ML_round(max(abs(ML_minmax(WW))),1,2);                             %
end                     % ------------------------------------------------- %
xmin=-xmax;                                                                 % lower bound
edges = xmin:h:xmax;
H1 = histc(WW,edges);                                                       % Histogram
% mH=ML_round(max(H1/(ntb*h)),1,2);
Chart=figure;
set(Chart,'visible','off');
axes('Parent',Chart,'FontSize',12); ML_FigureSize,hold on;
bH=bar(edges,H1/(ntb*h),'histc');
set(bH, 'FaceColor', [180,180,250]/255,'EdgeColor',[250,250,250]/255)
f = @(z) normpdf(z,0,1);                                                    % Overlay the distribution
fP=fplot(f,[xmin xmax],'r--','linewidth',2);                                % ------------------------
axis([xmin xmax 0 0.45])
set(gca,'ytick',0:.05:0.45,'yticklabel',num2str(get(gca,'ytick')','%.2f'))
gridxy(get(gca,'xtick'),get(gca,'ytick'),'color',grey,'linewidth',1);
hold off; box on;
set(gca,'Layer','top');