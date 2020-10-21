close all;
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')

a = combinedCompiledProjects_allEnhancers(strcmpi({combinedCompiledProjects_allEnhancers.dataSet}, '1Dg11_2xDl' )  &...
    [combinedCompiledProjects_allEnhancers.cycle]==12);

%%
fig = figure;
ax = axes(fig);
maxBin = max([a.dorsalFluoBin]);
for bin = 1:maxBin
    b = a([a.dorsalFluoBin]==bin);
    fluos = [b.particleFluo3Slice];
    fluos(fluos <= 0) = [];
    fluos(isnan(fluos)) = [];
 if ~isempty(fluos)
    histogram(log(fluos), 'Normalization', 'pdf', 'facealpha', .5)
    pd = fitdist(log(fluos)','Normal');
    x_values = 0:.1:10;
    y = pdf(pd,x_values);
    plot(ax, x_values,y,'-','LineWidth',.5)
    hold on
 end
end
title('log spot fluorescence divided by Dl bin. Normal fit')
ylabel('pdf')


%%
fig2 = figure;
ax2 = axes(fig2);
b = a;
fluos = [];
for k = 1:length(b)
    fluos = [fluos, min(b(k).particleFluo3Slice)];
end  
fluos(fluos <= 0) = [];
fluos(isnan(fluos)) = [];
 if ~isempty(fluos)
    histogram(ax2, fluos, 'Normalization', 'pdf', 'facealpha', .5)
    hold on;
    pd = fitdist(fluos','Lognormal');
    x_values = .1:.1:400;
    y = pdf(pd,x_values);
    plot(ax2, x_values,y,'-','LineWidth',.5)
 end
title('min spot fluorescence. Lognormal fit')
ylabel('pdf')

%%

fig3 = figure;
ax3 = axes(fig3);
    b = a;
    fluos = [b.particleFluo3Slice];
    fluos(fluos <= 0) = [];
    fluos(isnan(fluos)) = [];
    fluos = fluos./max(fluos(:));
 if ~isempty(fluos)
    histogram(ax3, fluos, 'Normalization', 'pdf', 'facealpha', .5)
    pd = fitdist(fluos','Lognormal');
    hold on
    x_values = 0:1E-3:1;
    y = pdf(pd,x_values);
    plot(ax3, x_values,y,'-','LineWidth',1)
    hold on
end
xlabel('spot fluo normalized to max')
ylabel('pdf')

coeffText = getDistributionText(pd);
title(["spot fluorescence fit to lognormal";coeffText'])

paperize(gcf, 4,4)