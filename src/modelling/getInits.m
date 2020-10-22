function [p0, lb, ub] = getInits(expmnt, md, metric, x_max, y_max, nSets, varargin)

minKD = 200;
maxKD = 1E4;
minw = 1E-3; %1E-2
maxw = 1E3; %1E2

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end


if expmnt == "affinities" && md=="simpleweak" && metric=="fraction" 
    %simplebindingweak_fraction- omegaDP kd1..kdn
    p0 = [1; [x_max/2;x_max;maxKD.*ones(1, nSets-2)']];
    lb = [minw; minKD.*ones(1, nSets)'];
    ub = [maxw; maxKD*ones(1, nSets)'];
elseif expmnt == "affinities" && md=="simpleweak" && metric=="fluo"
    %simplebindingweak_fluo- omegaDP kd1..kdn amp off
    p0 = [1; [x_max/2;x_max;maxKD.*ones(1, nSets-2)']; 500; 10];
    lb = [minw; minKD.*ones(1, nSets)'; 10; 10];
    ub = [maxw; maxKD*ones(1, nSets)'; 1E3; 1E3];
elseif expmnt == "phases" && md=="simpleweak" && metric=="fraction"
     %simplebindingweak_fraction- omegaDP1...omegaDPn kd
    p0 = [x_max/2; .02*ones(1, nSets)'];
    lb = [minKD; 1E-3*ones(1, nSets)'];
    ub = [maxKD; 1E3*ones(1, nSets)'];
elseif expmnt == "phases" && md=="simpleweak" && metric=="fluo"
     %simplebindingweak_fraction- omegaDP1...omegaDPn kd amp off
    p0 = [x_max/2; .02*ones(1, nSets)'; 500; 10];
    lb = [minKD; 1E-3*ones(1, nSets)'; 10; 10];
    ub = [maxKD; 1E3*ones(1, nSets)'; 1E3; 1E3];
end


%simple binding- amp kd1 kd2 kd3 p/kdp offset omegadp 
% p0 = [y_max; [x_max/2;x_max;1E5*ones(1, nSets-2)']; 1 ;0; 1];
% lb = [0; 200*ones(1, nSets)'; 0; 0; 1];
% ub = [min(1, y_max*2); Inf*ones(1, nSets)'; Inf; y_max*10; Inf];

%hill- amp kd1..kdn n offset
% p0 = [y_max;  [x_max/2;x_max;1E4*ones(1, nSets-2)'] ; 1 ; 1];
% lb = [0; 200*ones(1, nSets)';1; -y_max];
% ub = [y_max*2; 1E4*ones(1, nSets)'; 8; y_max*10];


%simple binding weak- amp kd1 kd2 kd3 offset omegadp
% p0 = [y_max; [x_max/2;x_max;1E4*ones(1, nSets-2)']; 0; .5];
% lb = [0; 200*ones(1, nSets)'; 0; 0];
% ub = [min(1, y_max*2); Inf*ones(1, nSets)'; y_max*10; Inf];

end