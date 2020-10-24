close all force;

[results, chain, s2chain] =fitstuff_mcmc2glob('nSimu', 1E3, 'expmnt', 'affinities', 'metric', 'fluo', 'minKD', 0,...
    'maxKD', 9000, 'noOff', true, 'minw', 1E-4, 'maxw', 1E4, 'minR', 1E-3, 'maxR', 1E4, 'displayFigures', false);

kd_1dg = results.mean(2)

[results, chain, s2chain] =fitstuff_mcmc2glob('nSimu', 1E3, 'expmnt', 'phases', 'metric', 'fluo',...
    'noOff', true, 'minw', 1E-4, 'maxw', 1E4, 'minR', 1E-3, 'maxR', 1E4, 'fixedKD', kd_1dg);
