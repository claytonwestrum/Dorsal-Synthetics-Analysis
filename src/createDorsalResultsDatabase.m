function dorsalResultsDatabase =...
    createDorsalResultsDatabase(dataTypes)

dorsalResultsDatabase = struct;

for i = 1:length(dataTypes)
    
    [~, resultsFolder, ~] = getDorsalPrefixes(dataTypes{i});
    
    if exist([resultsFolder,filesep,dataTypes{i},filesep,'dorsalResults.mat'], 'file')
        load([resultsFolder,filesep,dataTypes{i},filesep,'dorsalResults.mat'], 'dorsalResults')
    else
        warning(['skipping: ',dataTypes{i}])
        continue
    end
    
    dorsalResultsClean = addKeysToDorsalResults(dorsalResults);
    
    %skipping nc14
    dorsalResultsClean = dorsalResultsClean(1:2);
    
    for ncIndex = 1:length(dorsalResultsClean)
        
        fields = fieldnames(dorsalResultsClean{ncIndex});
        
        if isempty(fieldnames((dorsalResultsDatabase)))
            emptyCell = cell(length(fields),1);
            dorsalResultsDatabase= cell2struct(emptyCell,fields);
        end
        
        for f = 1:length(fields)
            
            
            if isvector(dorsalResultsClean{ncIndex}.(fields{f}))
                
                dorsalResultsDatabase.(fields{f}) =...
                    [dorsalResultsDatabase.(fields{f})
                    dorsalResultsClean{ncIndex}.(fields{f})];
                
            elseif isnumeric(dorsalResultsClean{ncIndex}.(fields{f}))
                
                %matrices require special care. concacatenate
                %and pad with NaNs since data sets
                %have different numbers of embryos
                dorsalResultsDatabase.(fields{f}) =...
                    padCatMatrices(dorsalResultsDatabase.(fields{f}),...
                    dorsalResultsClean{ncIndex}.(fields{f}) );
                
            end
            
        end
        
    end %nc loop
    
end %dataType loop

save([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], '-v6');