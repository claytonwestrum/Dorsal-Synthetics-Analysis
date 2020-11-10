function fluoMaxNoAverage(varargin)

expmnt = "affinities"; %affinities, phases or phaff
md = "simpleweak"; %simpleweak, simpleweakdimer, repression, tfdriven, artifact, fourstate
metric = "fluo"; %fraction, fluo
lsq = false;
noOff = true;
nSimu = 1E3; %1E3 is bad for real stats but good for debugging. need 1E4-1E6 for good stats
minKD = 200;
maxKD = 1E4;
minw = 1E-2; %1E-2
maxw = 1E1; %1E2
minR = 10;
maxR = 1E3;
displayFigures = true;
wb = true;
fixedKD = NaN; %if this value isn't nan, this value will determine the fixed KD parameter and KD won't be fitted
fixedOffset = NaN;
fixedR = NaN;
fixedw = NaN;
enhancerSubset = {};
scoreSubset = [];
positionSubset = [];
useBatches = false; %fit all the data across embryos and not just means 

%options must be specified as name, value pairs. unpredictable errors will
%occur, otherwise.
for i = 1:2:(numel(varargin)-1)
    if i ~= numel(varargin)
        eval([varargin{i} '=varargin{i+1};']);
    end
end

if expmnt == "phaff"
    lsq = false;
    noOff = true;
end

enhancers_1dg = {'1Dg11'};
enhancers_aff =  {'1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
enhancers_ph = {'1Dg-5', '1Dg-8D'};
scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
positions = [0, -5, -8]';

if strcmpi(expmnt, 'affinities')
    enhancers = [enhancers_1dg, enhancers_aff];
elseif strcmpi(expmnt, 'phases')
    enhancers = [enhancers_1dg, enhancers_ph];
elseif expmnt=="phaff"
    enhancers = [enhancers_1dg, enhancers_aff, enhancers_ph];
end

if ~isempty(enhancerSubset)
    enhancers = enhancerSubset;
    scores= scoreSubset;
    positions = positionSubset;
end


%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = getXRange(enhancers, expmnt);

nSets = length(enhancers);

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'])
% 
a = combinedCompiledProjects_allEnhancers(strcmpi({combinedCompiledProjects_allEnhancers.dataSet}, '1Dg11_2xDl' )  &...
    [combinedCompiledProjects_allEnhancers.cycle]==12);

b = a(cellfun(@any, {a.particleFrames}));

x = [b.dorsalFluoFeature];
y = cellfun(@max, {b.particleFluo3Slice});

% plot(x, y, '.', 'MarkerSize', 10)
% set(gca, 'YScale', 'log')
% ylim([1E2, 1E3])
y_max = nanmax(y);
x_max = max(x);
[p0, lb, ub, names] = getInits(expmnt, md, metric ,x_max, y_max, nSets,...
    'minKD', minKD, 'maxKD', maxKD, 'minw', minw, 'maxw', maxw, 'minR', minR, 'maxR', maxR, 'subset', ~isempty(enhancerSubset));

data.ydata = [x' y'];
k0 = p0;

if noOff
    %ignore the offset in the initial parameters and bounds
    p0(names=="offset") = [];
    lb(names=="offset") = [];
    ub(names=="offset") = [];
    k0(names=="offset") = [];
    names(names=="offset") = [];
end

% put the initial parameters and bounds in a form that the mcmc function
% accepts
params = cell(1, length(k0));
pri_mu = NaN; %default prior gaussian mean
pri_sig = Inf; %default prior gaussian variance
localflag = 0; %is this local to this dataset or shared amongst batches?


for i = 1:length(k0)
    
    targetflag = 1; %is this optimized or not? if this is set to 0, the parameter stays at a constant value equal to the initial value.
    
    if ~isnan(fixedKD) && contains(names(i), "KD")
        k0(i) = fixedKD;
        targetflag = 0;
    end
    
    if ~isnan(fixedOffset) && contains(names(i), "off")
        k0(i) = fixedOffset;
        targetflag = 0;
    end
    
    if ~isnan(fixedR) && contains(names(i), "R")
        k0(i) = fixedR;
        targetflag = 0;
    end
    
    if ~isnan(fixedw) && contains(names(i), "w")
        k0(i) = fixedw;
        targetflag = 0;
    end
    
    params{1, i} = {names(i),k0(i), lb(i), ub(i), pri_mu, pri_sig, targetflag, localflag};
    
end


model = struct;

simpleWeakOptions = struct('noOff', noOff, 'fraction', metric=="fraction",...
    'dimer', contains(md, "dimer"), 'expmnt', expmnt, 'onedsid', true);

mdl = @(x, p) simpleweak(x, p, simpleWeakOptions);


model.modelfun   = mdl;  %use mcmcrun generated ssfun 

if lsq
    model.sigma2 = mse;
end

options.drscale = 5; % a high value (5) is important for multimodal parameter spaces
options.waitbar = wb; %the waitbar is rate limiting sometimes
options.nsimu = nSimu; %should be between 1E3 and 1E6
options.updatesigma = 1; %honestly don't know what this does

rng(1,'twister'); %set the rng seed so we get the same results every run of this function

if lsq
    [results,chain,s2chain] = mcmcrun(model,data,params,options);
else
    %we're gonna run this three times and use the initial results of one
    %run as conditions for the next. this is an alternative when common least
    %squares gives results too poor to initialize with
    results = [];
    [results,~,~,~]=mcmcrun(model,data,params,options,results);
    [results,~,~,~]=mcmcrun(model,data,params,options,results);
    [results,chain,s2chain,~]=mcmcrun(model,data,params,options,results);
end

end