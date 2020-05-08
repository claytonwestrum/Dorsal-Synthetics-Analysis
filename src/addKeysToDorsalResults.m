function dorsalResults = addKeysToDorsalResults(dorsalResults)

 dataType = dorsalResults{1}.DataType;
 
    for ncIndex = 1:3
        
        %this field is problematic. 
         dorsalResults{ncIndex} = rmfield(dorsalResults{ncIndex}, 'allmrnasnomean');

        try
            setLength = length(dorsalResults{ncIndex}.dorsalFluoBins); 
            dorsalResults{ncIndex} = rmfield(dorsalResults{ncIndex}, 'DataType');
            [dorsalResults{ncIndex}.DataType{1:setLength}] = deal(dataType);
            dorsalResults{ncIndex}.DataType = dorsalResults{ncIndex}.DataType';
            dorsalResults{ncIndex}.nc = ones(setLength, 1)*(ncIndex + 11); %nc offset
            
            %transpose to conform with the rest of the data
            if size(dorsalResults{ncIndex}.dorsalFluoBins, 1) == 1
                dorsalResults{ncIndex}.dorsalFluoBins = dorsalResults{ncIndex}.dorsalFluoBins';
            end
            
        catch
            warning(['skipping: ', num2str(ncIndex)])
        end
        
    end
    
end