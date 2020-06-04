function averagedTimeTraces = averageCombinedCompiledProjects(DataType)

[~, resultsFolder] = getDorsalFolders;

load([resultsFolder,filesep,DataType,filesep,'combinedCompiledProjects.mat'], 'combinedCompiledProjects');
load([resultsFolder,filesep,DataType,filesep,'dlfluobins.mat'], 'dlfluobins');

nBins = length(dlfluobins);

binnedTime_minutes = 0:.1:30; %10s frame interval. should cover nc13 and then some.
nFrames = length(binnedTime_minutes);

combinedCompiledProjects = addTimeBins(binnedTime_minutes, combinedCompiledProjects);

nParticlesDV = zeros(nFrames, nBins);
vectorDV = nan(nFrames, nBins, length(combinedCompiledProjects));


nParticlesDV_aligned = zeros(nFrames, nBins);
vectorDV_aligned = nan(nFrames, nBins, length(combinedCompiledProjects));

for dlFluoBin = 1:nBins
    
    for k = 1:length(combinedCompiledProjects)
        
        if ~isempty(combinedCompiledProjects(k).particleFrames) &&...
                combinedCompiledProjects(k).dorsalFluoBin == dlFluoBin &&...
                combinedCompiledProjects(k).cycle == 12
            
            for frame = 1:length(combinedCompiledProjects(k).binnedTime_index)
                
                t = combinedCompiledProjects(k).binnedTime_index(frame);
                
                nParticlesDV(t, dlFluoBin) =...
                    nParticlesDV(t, dlFluoBin) + 1;
                
                vectorDV(t, dlFluoBin, k) = combinedCompiledProjects(k).particleFluo(frame);
                
                
                %same as above but particles are aligned for averaging
                t_aligned = combinedCompiledProjects(k).binnedTime_indexAligned(frame);
                nParticlesDV_aligned(t_aligned, dlFluoBin) =...
                    nParticlesDV(t_aligned, dlFluoBin) + 1;
                vectorDV_aligned(t_aligned, dlFluoBin, k) = combinedCompiledProjects(k).particleFluo(frame);
                
            end
            
        end
        
    end
    
end

meanVectorDV = nanmean(vectorDV, 3);
sdVectorDV = nanstd(vectorDV, 0, 3) ./ sqrt(nParticlesDV);

meanVectorDV_aligned = nanmean(vectorDV_aligned, 3);
sdVectorDV_aligned = nanstd(vectorDV_aligned, 0, 3) ./ sqrt(nParticlesDV_aligned);

averagedTimeTraces = struct;
averagedTimeTraces.meanVectorDV = meanVectorDV;
averagedTimeTraces.sdVectorDV = sdVectorDV;
averagedTimeTraces.nParticlesDV = nParticlesDV;
averagedTimeTraces.elapsedTime = binnedTime_minutes';
averagedTimeTraces.dlfluobins = dlfluobins;
averagedTimeTraces.unaveragedVectorDV = vectorDV; 

averagedTimeTraces.meanVectorDV_aligned = meanVectorDV_aligned;
averagedTimeTraces.sdVectorDV_aligned = sdVectorDV_aligned;
averagedTimeTraces.unaveragedVectorDV_aligned = vectorDV_aligned; 

save([resultsFolder,filesep,DataType,filesep,'averagedTimeTraces.mat'], 'averagedTimeTraces', '-v6');

figure;
tiledlayout('flow')
for dv = 1:length(averagedTimeTraces.dlfluobins)
    
    nexttile;
    errorbar(averagedTimeTraces.elapsedTime,...
        averagedTimeTraces.meanVectorDV(:,dv),...
        averagedTimeTraces.sdVectorDV(:, dv), 'CapSize', 0);
    
end

saveas(gcf, [resultsFolder, filesep, DataType, filesep, 'trapezoidMontage.fig']);
saveas(gcf, [resultsFolder, filesep, DataType, filesep, 'trapezoidMontage.png']);

figure;
tiledlayout('flow')
for dv = 1:length(averagedTimeTraces.dlfluobins)
    
    nexttile;
    errorbar(averagedTimeTraces.elapsedTime,...
        averagedTimeTraces.meanVectorDV_aligned(:,dv),...
        averagedTimeTraces.sdVectorDV_aligned(:, dv), 'CapSize', 0);
    xlim([0, 5])
    ylim([0, 300])
    
end

saveas(gcf, [resultsFolder, filesep, DataType, filesep, 'trapezoidAlignedMontage.fig']);
saveas(gcf, [resultsFolder, filesep, DataType, filesep, 'trapezoidAlignedMontage.png']);


end

function combinedCompiledProjects = addTimeBins(binnedTime_minutes, combinedCompiledProjects)

alignedStartFrame = 1; %arbitrary start point to align traces

for k = 1:length(combinedCompiledProjects)
    
    if ~isempty(combinedCompiledProjects(k).particleFrames)
        
        t = combinedCompiledProjects(k).particleTimeSinceAnaphase;
        
        for frame = 1:length(t)
            
            
            diffs = binnedTime_minutes - t(frame);
            [~, ind]...
                = min(abs(diffs));
            combinedCompiledProjects(k).binnedTime_index(frame) = ind;
            
            %add a version of the field where we align the traces
            if frame == 1
                firstInd = ind;
                combinedCompiledProjects(k).binnedTime_indexAligned(frame) = alignedStartFrame;
            else
                combinedCompiledProjects(k).binnedTime_indexAligned(frame) = alignedStartFrame + ind - firstInd;
            end
            
        end
        
        
    end
    
end


end


