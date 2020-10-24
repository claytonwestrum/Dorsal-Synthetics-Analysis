function [results,chain,s2chain]  = fitstuff_mcmc2glob(varargin)

[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')

expmnt = "affinities"; %affinities, phases or phaff
md = "simpleweak"; %only this
metric = "fraction"; %fraction or fluo
lsq = true;
noOff = false;
nSimu = 1E4;
minKD = 200;
maxKD = 1E4;
minw = 1E-2; %1E-2
maxw = 1E1; %1E2
minR = 10;
maxR = 1E3;
displayFigures = true;
wb = true;
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

if strcmpi(expmnt, 'affinities')
    enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
    scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
elseif strcmpi(expmnt, 'phases')
    enhancers = {'1Dg11', '1Dg-5', '1Dg-8D'};
    scores = [0, -5, -8]';
elseif expmnt=="phaff"
    enhancers_aff =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
    scores = [6.23, 5.81, 5.39, 5.13, 4.80, 4.73, 4.29]';
    enhancers_ph = {'1Dg11', '1Dg-5', '1Dg-8D'};
    positions = [0, -5, -8]';
    enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW', '1Dg-5', '1Dg-8D'};
end


%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = getXRange(enhancers, expmnt);

nSets = length(enhancers);
xo = {};
yo = {};
xs = {};
ys = {};
dsid = [];
T = [];
Y = [];
for k = 1:nSets
    cond = strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k});
    xo{k} = dorsalResultsDatabase.dorsalFluoBins(cond);
    if metric == "fraction"
        yo{k} = dorsalResultsDatabase.meanFracFluoEmbryo(cond);
    elseif metric == "fluo"
        yo{k} = dorsalResultsDatabase.meanAllMaxFluoEmbryo(cond);
    end
    
    [xs{k}, ys{k}]= processVecs(xo{k}, yo{k}, xrange(k, :));
    dsid = [dsid; k*ones(size(xs{k}))];
    
    if isnan(xrange(k, 1))
        xrange(k, 1) = min(xo{k});
    end
    if isnan(xrange(k, 2))
        xrange(k, 2) = max(xo{k});
    end
    
    assert(~isempty(xs{k}));
    
    T = [T; xs{k}];
    Y = [Y; ys{k}];
end

data.ydata = [T, Y];
data.dsid = dsid;
data.X =  [T dsid];

%%


%rate, kd, hill, y offset
y_max = nanmax(Y(:));
x_max = max(T);

[p0, lb, ub] = getInits(expmnt, md, metric ,x_max, y_max, nSets,...
    'minKD', minKD, 'maxKD', maxKD, 'minw', minw, 'maxw', maxw, 'minR', minR, 'maxR', maxR);


if lsq
    % Refine the first guess for the parameters with fminseacrh and calculate residual variance as an estimate of the model error variance.
    [k0, mse] = globfit2('expmnt',expmnt, 'metric', metric, 'md', md, 'maxKD', maxKD, 'displayFigures', displayFigures);
else
    k0 = p0;
end


if noOff && metric=="fluo"
    %ignore the offset in the initial parameters and bounds
    p0 = p0(1:end-1);
    lb = lb(1:end-1);
    ub = ub(1:end-1);
    k0 = k0(1:end-1);
end

% put the initial parameters and bounds in a form that the mcmc function
% accepts
params = cell(1, 3);
for i = 1:length(k0)
    params{1, i} = {['k', num2str(i)],k0(i), lb(i), ub(i)};
end

model = struct;

%%ssfun computes residuals for the mcmc function. mdl is used for computing
%%function values when plotting
mdl = getFitFuns(expmnt, md, metric, noOff);
model.ssfun = @(params, data) sum( (data.ydata(:,2)-mdl(data.X(:,1), params)).^2 );

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
    [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
    [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
    [results,chain,s2chain,sschain]=mcmcrun(model,data,params,options,results);
end

if displayFigures
    
    burnInTime = .25; %let's burn the first 25% of the chain just in case
    chain = chain(round(burnInTime*nSimu):nSimu, :);
    if ~isempty(s2chain)
        s2chain = s2chain(round(.25*nSimu):nSimu, :);
    end
    
    chainfig = figure(); clf
    mcmcplot(chain,[],results,'chainpanel')
    
    % Function chainstats lists some statistics, including the estimated Monte Carlo error of the estimates.
    %geweke is a measure of whether things converged between 0 and 1.
    chainstats(chain,results)
    
    
    
    %ideally, these guys look like ellipses. if certain parameters give weird
    %shapes, it might mean those parameters should be removed from the model if
    %possible
    pairFig = figure; clf
    mcmcplot(chain,[],results,'pairs', .5);
    %
    
    
    figure(4); clf
    til = tiledlayout(1, nSets);
    dsid2 = [];
    
    xx = (0:10:max(data.X(:,1)))';
    
    for k = 1:nSets
        dsid2 = [dsid2; k*ones(length(xx), 1)];
    end
    
    
    %get the prediction intervals for the parameters and function vals
    out = mcmcpred(results,chain,[],repmat(xx, nSets, 1), mdl);
    
    % mcmcpredplot(out);
    nn = (size(out.predlims{1}{1},1) + 1) / 2;
    plimi = out.predlims{1};
    yl = plimi{1}(1,:);
    yu = plimi{1}(2*nn-1,:);
    
    km = mean(chain);
    ks = std(chain);
    X2 = [repmat(xx, nSets, 1), dsid2];
    
    
    yf = plimi{1}(nn,:);
    
    for i = 1:nSets
        
        yy = yf(X2(:, 2)==i);
        yyl = yl(X2(:, 2)==i);
        yyu = yu(X2(:, 2)==i);
        nexttile;
        fillyy(xx,yyl,yyu,[0.9 0.9 0.9]);
        hold on
        plot(xo{i}, yo{i}, 'o-', 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k', 'Color', 'r');
        plot( xo{i}(xo{i} >= xrange(i, 1) & xo{i} <=xrange(i, 2) ), yo{i}(xo{i} >=  xrange(i, 1) & xo{i} <= xrange(i, 2)),'o-','MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'Color', 'r');
        plot(xx,yy,'-k')
        xlim([0,3500])
        
        if expmnt == "affinities"
            vartheta = 'KD = ';
            consttheta = ' \omega'' = ';
        elseif expmnt == "phases"
            consttheta = 'KD = ';
            vartheta = ' \omega'' = ';
        end
        
        
        if expmnt == "phaff"
            naff = 7;
            nph = 3;
            if i <= naff
                
                titleCell = {enhancers{i}, [' \omega'' = ', num2str(round2(km(1))), ' \pm ', num2str(round2(ks(1)))],...
                    [ 'KD = ', num2str(round2(km(nph+i))), ' \pm ', num2str(round2(ks(nph+i)))]};
            elseif i > naff
                titleCell = {enhancers{i}, [' \omega'' = ', num2str(round2(km(i - (naff-1) ))), ' \pm ', num2str(round2(ks(i - (naff-1))))],...
                    [ 'KD = ', num2str(round2(km(nph+1))), ' \pm ', num2str(round2(ks(nph+1)))]};
            end
        else
            titleCell = {enhancers{i}, [vartheta, num2str(round2(km(i+1))), ' \pm ', num2str(round2(ks(i+1)))],...
                [consttheta, num2str(round2(km(1))), ' \pm ', num2str(round2(ks(1)))]};
        end
        
        if metric=="fraction"
            ylim([0, 1])
        elseif metric=="fluo"
            titleCell = [titleCell,...
                [' amp = ', num2str(round2(km(nSets + 1))), ' \pm ', num2str(round2(ks(nSets + 1)))] ];
            if ~noOff
                titleCell = [titleCell, [' off = ', num2str(round2(km(end))), ' \pm ', num2str(round2(ks(end)))] ];
            end
            ylim([0, y_max])
        end
        title(titleCell);
    end
    
    
    figure;
    if expmnt == "phaff"
        
        tilo = tiledlayout('flow');
        nexttile;
        errorbar(scores, km(nph+1:nph+naff), ks(nph+1:nph+naff));
        ylabel('KD (au)')
        xlabel('affinity (Patser score)')
        
        nexttile;
        errorbar(positions, km(1:nph), ks(1:nph));
        ylabel('w'' (au)')
        xlabel('position (bp)')
        xlim([-10, 2])
    else
        
        errorbar(scores, km(2:nSets+1), ks(2:nSets+1));
        if expmnt == "affinities"
            ylabel('KD (au)')
            xlabel('affinity (Patser score)')
        elseif expmnt == "phases"
            ylabel('w'' (au)')
            xlabel('position (bp)')
            xlim([-10, 2])
        end
        
    end
    
    
    covfig = figure;
    cv = @(x, y) sqrt(abs(x)) ./ sqrt((y'*y));
    imagesc(cv(results.cov, results.mean));
    colorbar;
    ylabel('parameter 1')
    xlabel('parameter 2')
    title('Covariance matrix of the parameters');
    
    colormap(viridis);
    
    
    
end

end
%%