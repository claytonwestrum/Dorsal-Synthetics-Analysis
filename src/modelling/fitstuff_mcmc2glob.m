function [results,chain,s2chain]  = fitstuff_mcmc2glob()

close all;
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')
enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};



%we're going to restrict the range of the fits specifically for each
%enhancer. nan values means no restriction
xrange = nan(length(enhancers), 2);
xrange(1,:)  = [1000, 2250];  %1Dg
xrange(2, :) = [500, 1750]; %upper limit for %1DgS
xrange(3, 1) = 500; %lower limit for 1DgW
xrange(4, 2) = 1500; %upper limit for %1DgAW


nSets = length(enhancers);
xo = {};
yo = {};
xs = {};
ys = {};
dsid = [];
T = [];
Y = [];
for k = 1:nSets
    xo{k} = dorsalResultsDatabase.dorsalFluoBins( ...
        strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k}) );
    yo{k} = dorsalResultsDatabase.meanFracFluoEmbryo( ...
        strcmpi(dorsalResultsDatabase.mother,'2x') & strcmpi(dorsalResultsDatabase.enhancer, enhancers{k} ) );
    [xs{k}, ys{k}]= processVecs(xo{k}, yo{k}, xrange(k, :));
    dsid = [dsid; k*ones(size(xs{k}))];
    T = [T; xs{k}];
    Y = [Y; ys{k}];
end

data.ydata = [T, Y];
data.dsid = dsid;
data.X =  [T dsid];

%%
% Refine the first guess for the parameters with fminseacrh and calculate residual variance as an estimate of the model error variance.


%rate, kd, hill, y offset
y_max = nanmax(Y(:));
x_max = max(T);
%hill- amp kd1..kdn n offset
% p0 = [y_max;  [x_max/2;x_max;1E4*ones(1, nSets-2)'] ; 1 ; 1];
% lb = [0; 200*ones(1, nSets)';1; -y_max];
% ub = [y_max*2; Inf*ones(1, nSets)'; 8; y_max*10];

%simplebindingweak_fraction- omegaDP kd1..kdn
% p0 = [.02; [x_max/2;x_max;1E4*ones(1, nSets-2)']];
% lb = [0; 200*ones(1, nSets)'];
% ub = [Inf; Inf*ones(1, nSets)'];


[k0, mse] = globfit2;


%
params = cell(1, 3);
for i = 1:length(k0)
    params{1, i} = {['k', num2str(i)],k0(i), 0};
end

% params{1, 1}{4} = 1;
%set max n of hill function to 10
% params{1, nSets+2}{4} = 10;
for i = 2:nSets+1
    params{1, i}{4} = 1E4;
end

model = struct;
model.ssfun = @funss_simpleweakfraction;
model.sigma2 = mse;

options.nsimu = 10000;
options.updatesigma = 1;
[results,chain,s2chain] = mcmcrun(model,data,params,options);

chainfig = figure(); clf
mcmcplot(chain,[],results,'chainpanel')
% subplot(2,2,4)
% mcmcplot(sqrt(s2chain),[],[],'dens',2)
% title('error std')
% ylabel('pdf')
% xlabel('RMSE')

% Function chainstats lists some statistics, including the estimated Monte Carlo error of the estimates.

chainstats(chain,results)



%ideally, these guys look like ellipses. if certain parameters give weird
%shapes, it might mean those parameters should be removed from the model if
%possible
figure(3); clf
mcmcplot(chain,[],results,'pairs', .5);
% 

figure(4); clf
til = tiledlayout(1, nSets);
dsid2 = [];

xx = (0:10:max(data.X(:,1)))';
for k = 1:nSets
    dsid2 = [dsid2; k*ones(length(xx), 1)];
end

out = mcmcpred(results,chain,[],repmat(xx, nSets, 1), @(x, p)subfun_simplebinding_weak_fraction_std2(x, p));
nn = (size(out.predlims{1}{1},1) + 1) / 2;
plimi = out.predlims{1};
yl = plimi{1}(1,:);
yu = plimi{1}(2*nn-1,:);
  
km = mean(chain);
ks = std(chain);
X2 = [repmat(xx, nSets, 1), dsid2];
yfit2 = subfun_simplebinding_weak_fraction(km, X2);

for i = 1:nSets

    yy = yfit2(X2(:, 2)==i);
    yyl = yl(X2(:, 2)==i);
    yyu = yu(X2(:, 2)==i);
    nexttile;
    fillyy(xx,yyl,yyu,[0.9 0.9 0.9]);
    hold on
    plot(xo{i},yo{i},'o-',xx,yy,'-')

%     xlim([min(xx), max(xx)*1.2])
%     ylim([min(yy), max(yy)*1.2]);
    xlim([0,3500])
    ylim([0, 1])
    title({enhancers{i}, ['KD = ' num2str(round2(km(i+1))), ' \pm ', num2str(round2(ks(i+1)))],...
        [' \omega'' = ', num2str(round2(km(1))), ' \pm ', num2str(round2(ks(1)))]})
end
title(til, 'global fit with mcmc')


end

function ss = funss(k, data)
% sum-of-squares
Aobs = data.ydata(:,2);


x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = k;

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
n = params(nSets + 2);
offset = params(nSets + 3);

Amodel = amplitude.*(((x./KD(dsid)).^n)./(1+((x)./KD(dsid)).^n))+offset;

ss = sum((Aobs-Amodel).^2);
end


function ss = funss_simpleweak(k, data)
% sum-of-squares
Aobs = data.ydata(:,2);


x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = k;

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
omegaDP = params(nSets + 2);
offset = params(nSets + 3);

Amodel = amplitude.*(((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP))+offset;

ss = sum((Aobs-Amodel).^2);
end

function yfit = subfun_simplebinding_weak(params,X)

%simplebinding in the weak promoter limit.

x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
omegaDP = params(nSets + 2);
offset = params(nSets + 3);

yfit = amplitude.*(((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP))+offset;

end


function ss = funss_simpleweakfraction(k, data)
% sum-of-squares
Aobs = data.ydata(:,2);
x = data.X(:,1);        % unpack time from X
dsid = data.X(:,2);     % unpack dataset id from X
nSets = max(dsid);
params = k;
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:nSets+1)';
Amodel = (((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP));
ss = sum((Aobs-Amodel).^2);
end

function yfit = subfun_simplebinding_weak_fraction(params,X)
%simplebinding in the weak promoter limit.
x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:nSets+1)';
yfit = (((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP));
end


function yfit = subfun_simplebinding_weak_fraction_std2(x, params)
%simplebinding in the weak promoter limit.

X = [];
n = 1;
X(1, :) = [x(1), 1];
for k = 2:length(x)
    if x(k) < x(k-1)
        n = n + 1;
    end
    X(k, 1) = x(k);
    X(k, 2) = n;
end

dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

nSets = max(dsid);
params = params(:)'; %need a row vec
omegaDP = params(1);
KD = params(2:nSets+1)';
yfit = (((x./KD(dsid)).*omegaDP)./(1+x./KD(dsid)+(x./KD(dsid)).*omegaDP));

end
