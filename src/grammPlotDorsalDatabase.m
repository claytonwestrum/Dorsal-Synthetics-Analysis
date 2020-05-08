function grammPlotDorsalDatabase()

dataTypes = {'1Dg_2xDl', '1DgW_2x_Leica',...
    '1DgW_FFF', '1Dg', '1Dg-5_FFF', '1DgVW_FFF'};

try
    [~, resultsFolder, ~] = getDorsalPrefixes(dataTypes{1});
    load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')
catch
    dorsalResultsDatabase = createDorsalResultsDatabase(dataTypes);
end

% color grouping data (number of cylinders) and select a subset of the data
% g=gramm('x',dorsalDatabase.dorsalFluoBins,...
%     'y',dorsalDatabase.meanFracFluoEmbryo,...
%     'color', dorsalDatabase.DataType, 'subset', dorsalDatabase.nc==12);

g=gramm('x',dorsalResultsDatabase.dorsalFluoBins,...
    'y',dorsalResultsDatabase.meanFracFluoEmbryo,...
     'subset', dorsalResultsDatabase.nc==12);


% Subdivide the data in subplots horizontally by region of origin
g.facet_grid([],dorsalResultsDatabase.DataType)
% Plot raw data as points
g.geom_point()
% Plot linear fits of the data with associated confidence intervals
g.stat_glm()
% Set appropriate names for legends
% g.set_names(dorsalDatabase.DataType)
% g.set_title('1Dg_2xDl nc12')

grammFigurePBoC(g);

g.draw()



end