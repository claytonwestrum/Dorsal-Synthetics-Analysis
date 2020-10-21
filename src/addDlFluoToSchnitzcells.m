function schnitzcells = addDlFluoToSchnitzcells(Prefix)

displayFigures = false;

if displayFigures
%     tileFig = figure(); %commented out because it's probably not useful
    holdFig = figure();
end

liveExperiment = LiveExperiment(Prefix);
schnitzcells = getSchnitzcells(liveExperiment);
schnitzcellsOld = schnitzcells;
ncFrames = [zeros(1, 8), liveExperiment.anaphaseFrames'];
FrameInfo = getFrameInfo(liveExperiment);
imSize = [FrameInfo(1).LinesPerFrame, FrameInfo(1).PixelsPerLine];
resultsFolder = liveExperiment.resultsFolder;

fluoFeatures = [];
for s = 1:length(schnitzcells)
    midFrame = ceil(length(schnitzcells(s).frames)/2);
    dif = double(schnitzcells(s).frames(midFrame)) - ncFrames;
    cycle = find(dif>0, 1, 'last' );
    schnitzcells(s).cycle = uint8(cycle);
end

schnitzcells = addRelativeTimeToSchnitzcells(schnitzcells, FrameInfo, ncFrames);

%do some general QC filtering
try
    schnitzcells = filterSchnitz(schnitzcells, imSize);
catch
    integrateSchnitzFluo(Prefix);
    schnitzcells = filterSchnitz(schnitzcells, imSize);
end

%do some dv specific filtering on the approved schnitzes
schnitzcells = filterSchnitzFurther(...
    schnitzcells);

approvedSchnitzes = find([schnitzcells.Approved]);

for s = approvedSchnitzes
    
    schnitzcells(s).FluoTimeTrace = single(ExtractDlFluo(schnitzcells(s).Fluo, .5));
    
    if schnitzcells(s).cycle ~= 14
        midCycle = floor((ncFrames(schnitzcells(s).cycle) + ncFrames(schnitzcells(s).cycle+1))/2);
    else
        midCycle = 1;
    end
    midCycleFrame = find(schnitzcells(s).frames==midCycle);
    schnitzcells(s).midCycleFrame = find(schnitzcells(s).frames==midCycle);
    
    schnitzcells(s).FluoTimeTraceSmooth = smooth(single(ExtractDlFluo(schnitzcells(s).Fluo, .5)), 5);
    
    if schnitzcells(s).cycle ~= 14
        midCycleSmooth = floor((ncFrames(schnitzcells(s).cycle) + ncFrames(schnitzcells(s).cycle+1))/2);
    else
        midCycleSmooth = 1;
    end
    midCycleFrameSmooth = find(schnitzcells(s).frames==midCycleSmooth);
    
    
    %different characterizations of intensity
    schnitzcells(s).fluoMid= single(schnitzcells(s).FluoTimeTrace(midCycleFrame)); %middle fluorescence
    schnitzcells(s).fluoMidSmooth = single(schnitzcells(s).FluoTimeTraceSmooth(midCycleFrameSmooth)); %middle fluorescence smoothed
    
    schnitzcells(s).FluoFeature =  schnitzcells(s).fluoMid; %the preferred intensity
    
    if isempty(schnitzcells(s).FluoFeature)
        schnitzcells(s).FluoFeature = NaN;
    end
    
    schnitzcellsSizeUnchanged(schnitzcellsOld, schnitzcells);
    
    
    fluoFeatures = [fluoFeatures, schnitzcells(s).FluoFeature]; %#ok<AGROW>
    
    if displayFigures && schnitzcells(s).cycle == 12 && ~isnan(schnitzcells(s).FluoFeature)
        
%         figure(tileFig)
%         nexttile;
%         plot(schnitzcells(s).frames, schnitzcells(s).FluoTimeTrace, '-k');
%         hold on
%         plot(midCycle, schnitzcells(s).FluoFeature, 'ob');
%         hold off
%         xticks([]);
%         yticks([round(min(schnitzcells(s).FluoTimeTrace)), round(max(schnitzcells(s).FluoTimeTrace))]);
        %                 ylim([0, 3000]);
        
        figure(holdFig)
%         plot(1:length(schnitzcells(s).FluoTimeTrace), schnitzcells(s).FluoTimeTrace, '-k');
        plot(schnitzcells(s).timeSinceAnaphase, schnitzcells(s).FluoTimeTrace, '-k');
        hold on
%         plot(midCycleFrame, schnitzcells(s).FluoFeature, 'ob');
        ylim([0, 3000]);
        
    end
    
    
end

save([resultsFolder,filesep,Prefix,'_lin.mat'], 'schnitzcells');

end