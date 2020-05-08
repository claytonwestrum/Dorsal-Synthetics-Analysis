function dorsalResultsDatabase =...
    createDorsalResultsDatabase(dataTypes)

dorsalResultsDatabase = struct;

for i = 1:length(dataTypes)
    
    try
        [~, resultsFolder, ~] = getDorsalPrefixes(dataTypes{i});
        load([resultsFolder,filesep,dataTypes{i},filesep,'dorsalResults.mat'], 'dorsalResults')
    catch
        warning(['skipping: ',dataTypes{i}])
        continue
    end
 
    dorsalResultsClean = addKeysToDorsalResults(dorsalResults);
    
    for ncIndex = 1:length(dorsalResultsClean)
        
        fields = fieldnames(dorsalResultsClean{ncIndex});
        
          if isempty(fieldnames((dorsalResultsDatabase)))
                emptyCell = cell(length(fields),1);
                dorsalResultsDatabase= cell2struct(emptyCell,fields);
          end
        
        for f = 1:length(fields)
                        
             %only concatenate vectors, not matrices
            if isvector(dorsalResultsClean{ncIndex}.(fields{f}))
                dorsalResultsDatabase.(fields{f}) =...
                    [dorsalResultsDatabase.(fields{f}) 
                    dorsalResultsClean{ncIndex}.(fields{f})];
            end
        end
        
    end %nc loop
        
end %dataType loop

save([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], '-v6');