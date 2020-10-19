function [b, mse] = globfit2
% Set up data so that Y is a function of T with a specific functional form,
% but there are multiple groups and one parameter varies across groups.

close all;
[~, resultsFolder] = getDorsalFolders;
load([resultsFolder, filesep, 'dorsalResultsDatabase.mat'], 'dorsalResultsDatabase')
% enhancers =  {'1Dg11', '1DgS2', '1DgW', '1DgAW3', '1DgSVW2', '1DgVVW3', '1DgVW'};
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

gscatter(T,Y,dsid)

%rate, kd, pkd, offset, omegadp
y_max = nanmax(Y);
x_max = max(T);
%simple binding- amp kd1 kd2 kd3 p/kdp offset omegadp 
% p0 = [y_max; [x_max/2;x_max;1E5*ones(1, nSets-2)']; 1 ;0; 1];
% lb = [0; 200*ones(1, nSets)'; 0; 0; 1];
% ub = [min(1, y_max*2); Inf*ones(1, nSets)'; Inf; y_max*10; Inf];

%hill- amp kd1..kdn n offset
p0 = [y_max;  [x_max/2;x_max;1E5*ones(1, nSets-2)'] ; 1 ; 1];
lb = [0; 200*ones(1, nSets)';1; -y_max];
ub = [y_max*2; Inf*ones(1, nSets)'; 8; y_max*10];

%simple binding weak- amp kd1 kd2 kd3 offset omegadp
% p0 = [y_max; [x_max/2;x_max;1E4*ones(1, nSets-2)']; 0; .5];
% lb = [0; 200*ones(1, nSets)'; 0; 0];
% ub = [min(1, y_max*2); Inf*ones(1, nSets)'; y_max*10; Inf];


%%
% Pack up the time and dataset id variables into X for later unpacking
X = [T dsid];

optimoptions = optimset('TolFun',1E-6, 'MaxIter', 1E6);
[b, resnorm, res] = lsqcurvefit(@subfun_hill,p0,X,Y, lb, ub, optimoptions);

mse = mean(res.^2);

tiledlayout(1, nSets);
dsid2 = [];

xx = (0:1:max(X(:,1)))';
for k = 1:nSets
    dsid2 = [dsid2; k*ones(length(xx), 1)];
end

X2 = [repmat(xx, nSets, 1), dsid2];
yfit2 = subfun_hill(b, X2);

for k = 1:nSets
    
    y = yfit2(X2(:, 2)==k);
    
    nexttile;
    plot(xo{k}, yo{k}, xx, y);
    ylim([0, max(yfit2)*1.1]);
    xlim([0, max(xx)]);
    title(enhancers{k});
    
end

end



function yfit = subfun_simplebinding(params,X)

x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
n = params(nSets + 2);
omegaDP = params(nSets + 3);
offset = params(nSets + 4);

yfit = amplitude.*((n+(x./KD(dsid)).*n.*omegaDP)./(1+x./KD(dsid)+n+(x./KD(dsid)).*n.*omegaDP))+offset;

end


function yfit = subfun_hill(params,X)

x = X(:,1);        % unpack time from X
dsid = X(:,2);     % unpack dataset id from X
nSets = max(dsid);

params = params(:)'; %need a row vec

amplitude = params(1);
KD = params(2:nSets+1)';
n = params(nSets + 2);
offset = params(nSets + 3);

yfit = amplitude.*(((x./KD(dsid)).^n)./(1+((x)./KD(dsid)).^n))+offset;

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

% %% example function if i want to use different functional forms for
% different datasets
% function yfit = subfun(params,X)
% T = X(:,1);        % unpack time from X
% dsid = X(:,2);     % unpack dataset id from X
% A0 = params(1);    % same A0 for all datasets
% A1 = params(2:4)'; % different A1 for each dataset
% tau = params(5);   % same tau
% yfit = zeros(size(T));
% % Find separate groups and call the F function for each group
% idx = (dsid==1);
% yfit(idx) = F1(X(idx),[A0 A1(1) tau]);
% idx = (dsid==2);
% yfit(idx) = F2(X(idx),[A0 A1(2) tau]);
% idx = (dsid==3);
% yfit(idx) = F3(X(idx),[A0 A1(3) tau]);
% end
% % Below are the three functions for the three groups
% function y = F1(x,b)
% y = b(1) + b(2)*exp(-x/b(3));
% end
% function y = F2(x,b)
% y = b(1) + b(2)*sin(x/b(3));
% end
% function y = F3(x,b)
% y = b(1) + b(2)*cos(x/(2*b(3)));
% end