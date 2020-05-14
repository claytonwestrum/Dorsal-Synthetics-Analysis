function grammPlotDorsalDatabase()

dataTypes = {'1Dg-8D_FFF', '1Dg11_2xDl', '1DgW_2x_Leica',...
    '1DgW_FFF', '1Dg11_og', '1Dg-5_FFF', '1DgVW_FFF', '1Dg11_FFF'};

dataTypes = {'1Dg-8D_FFF', '1DgW_2x_Leica',...
    '1DgW_FFF'};

try
    [~, resultsFolder, ~] = getDorsalPrefixes(dataTypes{1});
    load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')
catch
    dorsalResultsDatabase = createDorsalResultsDatabase(dataTypes);
end

g=gramm('x',dorsalResultsDatabase.dorsalFluoBins,...
    'y',dorsalResultsDatabase.meanFracFluoEmbryo,...
     'subset', dorsalResultsDatabase.nc==12, 'color', dorsalResultsDatabase.mother);

% Subdivide the data in subplots horizontally by region of origin
g.facet_grid(dorsalResultsDatabase.enhancer, [])
g.facet_wrap(dorsalResultsDatabase.enhancer, 'ncols', 2);
% Plot raw data as points
g.geom_point()
% Plot linear fits of the data with associated confidence intervals
g.stat_glm()
% g.stat_summary('type', 'sem');

 %set axis labels and legend labels
 g.set_names('x','Dorsal Concentration (au)','y','Fraction Active', 'row', '', 'column', '');

grammFigurePBoC(g);

g.draw()

end