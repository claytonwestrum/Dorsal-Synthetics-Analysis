function dorsalResults = compileAllProjects(DataType)


thisProject = LiveProject(DataType) %#ok<NOPRT>

[~, resultsFolder] = getDorsalFolders;




prefixes = thisProject.includedExperimentNames;

compiledProjects = cell(1, length(prefixes));
% 
% 
for k = 1:length(prefixes)
%     integrateSchnitzFluo(prefixes{k});
%     TrackmRNADynamics(prefixes{k});
    CompileParticles(prefixes{k},  'minBinSize', 0, 'MinParticles', 0,...
        'yToManualAlignmentPrompt');
    alignCompiledParticlesByAnaphase(prefixes{k});
end

%% Validate experiments included in the analysis
hasAllPushed = [thisProject.hasSpots, thisProject.hasParticles,...
    thisProject.hasSchnitzcells, thisProject.hasCompiledParticles,...
    thisProject.anaphaseFramesAnnotated];

% assert( all(hasAllPushed) );
% 
for k = 1:length(prefixes)
    checkSchnitzAssignmentToParticles(prefixes{k})
end
% 
addDVStuffToSchnitzCells(DataType)

binDorsal(DataType, false)


for k = 1:length(prefixes)
    
    compiledProjects{k} = makeCompiledProject(prefixes{k});
    
    if k == 1
        combinedCompiledProjects = compiledProjects{k}; 
    else
        try
            combinedCompiledProjects = [combinedCompiledProjects, compiledProjects{k}]; %#ok<AGROW>
        catch
            [combinedCompiledProjects, compiledProjects{k}] = addFields(combinedCompiledProjects,compiledProjects{k});
            combinedCompiledProjects = [combinedCompiledProjects, compiledProjects{k}]; %#ok<AGROW>
        end
        
    end
    
end

[combinedCompiledProjects.dataSet] = deal(DataType);

save([resultsFolder,filesep,DataType,filesep,'combinedCompiledProjects.mat'], 'combinedCompiledProjects');

% averagedTimeTraces = averageCombinedCompiledProjects(DataType, true);

dorsalResults = createDorsalResults(DataType); 

% plotDorsalResultsLoop(DataType, 'frac', 1:6, 'hill')
% plotDorsalResultsLoop(DataType, activity)

end