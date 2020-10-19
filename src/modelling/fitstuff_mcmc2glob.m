function fitstuff_mcmc2glob()

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
X = [T dsid];
data.X = X;
%%
% Refine the first guess for the parameters with fminseacrh and calculate residual variance as an estimate of the model error variance.
modd = standardizeModelForGramm(dorsalFitFunction('hill'));
model2 = dorsalFitFunction('hill');
x = xs{1};
y = ys{1};
%rate, kd, hill, y offset
y_max = nanmax(Y(:));
x_max = max(T);
%hill- amp kd1..kdn n offset
p0 = [y_max;  [x_max/2;x_max;1E5*ones(1, nSets-2)'] ; 1 ; 1];
lb = [0; 200*ones(1, nSets)';1; -y_max];
ub = [y_max*2; Inf*ones(1, nSets)'; 8; y_max*10];


[k0, mse] = globfit2;


%
params = cell(1, 3);
for i = 1:length(k0)
    params{1, i} = {['k', num2str(i)],k0(i), 0};
end

params{1, 1}{4} = 1;
%set max n of hill function to 10
params{1, nSets+2}{4} = 10;

model = struct;
model.ssfun = @funss;
model.sigma2 = mse;

options.nsimu = 100000;
options.updatesigma = 1;
[results,chain,s2chain] = mcmcrun(model,data,params,options);

figure(1); clf
mcmcplot(chain,[],results,'chainpanel')
subplot(2,2,4)
mcmcplot(sqrt(s2chain),[],[],'dens',2)
title('error std')

% Function chainstats lists some statistics, including the estimated Monte Carlo error of the estimates.

chainstats(chain,results)



%ideally, these guys look like ellipses. if certain parameters give weird
%shapes, it might mean those parameters should be removed from the model if
%possible
figure(3); clf
mcmcplot(chain,[],results,'pairs', .5);
% 
% figure(2); clf
% xx = (0:1:1E4)';
% yy = model2(mean(chain), xx);
% plot(xo{1},yo{1},'s',xx,yy,'-')
% xlim([min(xx), max(xx)])
% ylim([min(yy), max(yy)]);
% legend({'obs','fit'},'Location','best')
% title('Data and fitted model')



% commented out because these graphs are inside the pairs figure
% figure(4); clf
% mcmcplot(chain,[],results,'denspanel',2);

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