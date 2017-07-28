%% Problem 2: Weighted averaging

means = [50, 55, 40, 36, 51];
stds = [8, 10, 8, 10, 12];
weights = 1./(stds.^2);
numer = sum(weights.*means);
denom = sum(weights);
wam = numer./denom
wam_uncertainty = numer