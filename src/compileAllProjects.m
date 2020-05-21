function dorsalResults = compileAllProjects(DataType)


thisProject = liveProject(DataType) %#ok<NOPRT>


%% Validate experiments included in the analysis
hasAllPushed = [thisProject.hasSpots, thisProject.hasParticles,...
    thisProject.hasSchnitzcells, thisProject.hasCompiledParticles,...
    thisProject.anaphaseFramesAnnotated];

assert( all(hasAllPushed) ); 

%%

prefixes = thisProject.includedExperimentNames;

compiledProjects = cell(1, length(prefixes));

% 
% for k = 1:length(prefixes)
%     %to prevent issues with stuff in memory
%     clear getMovieMat;
%     clear getHisMat;
% %     fit3DGaussiansToAllSpots(prefixes{k}, 1);
% %     integrateSchnitzFluo(prefixes{k});
%   %  TrackmRNADynamics(prefixes{k});
%     CompileParticles(prefixes{k},  'minBinSize', 0, 'MinParticles', 0,...
%         'yToManualAlignmentPrompt');
%     alignCompiledParticlesByAnaphase(prefixes{k});
% end

addDVStuffToSchnitzCells(DataType)

binDorsal(DataType, false)

for k = 1:length(prefixes)
    compiledProjects{k} = makeCompiledProject(prefixes{k});
end

dorsalResults = createDorsalResults(DataType); 

% plotDorsalResultsLoop(DataType, 'frac', 1:6, 'hill')
% plotDorsalResultsLoop(DataType, activity)

end