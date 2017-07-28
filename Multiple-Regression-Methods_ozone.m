fname = load('Atlanta_O3.txt');
days = fname(:,2);
days(32:62)=days(32:62)+31;
obs = fname(:,3);
sim = fname(:,4);
sim(obs<0) = []; % this order is important
days(obs<0) = [];
obs(obs<0) = [];
endval = length(sim);

% Introducing an outlier: 
% sim(endval) = 30;
% obs(endval) = 70;


%% Regression methods
% x = observed; y = simulated

% LS REGRESSION
plot(obs, sim, 'g*'); hold on;
p = polyfit(obs, sim,1);
fit_sims = polyval(p,obs);
plot(obs, fit_sims, 'k-');

    % Why do people ever use normal linear regression? Why don't we always
    % use reduced-axis regression and account for errors in the x- and y-?
    % Or do more scientists use reduced-axis regression (when even x-
    % variable is not guaranteed) than I realize? 

% Reduced major axis
beta1 = std(sim)./std(obs);
beta0 = mean(sim)-(beta1.*mean(obs));
reduc_sims = polyval([beta1 beta0],obs);
plot(obs, reduc_sims, 'r')

% Principle component 
means_out = [obs-mean(obs) sim-mean(sim)];
[u,s,v] = svd(means_out, 'econ');
v=v';
obs_pc = u(:,1)*s(1,1)*v(1,1)+mean(obs);
sim_pc = u(:,1)*s(1,1)*v(1,2)+mean(sim);
coeffs_pc = polyfit(obs_pc, sim_pc, 1);
pc_sims = polyval(coeffs_pc,obs);
plot(obs, pc_sims, 'b');

xlabel('Observed ozone (ppbv)')
ylabel('Simulated ozone (ppbv)')
legend('Data', 'Linear-Regression', 'Reduced major axis', 'Principle component')
hold off

% slopes and intercepts
sprintf('Lin reg slope is %d and intercept is %d', p(1), p(2))
% sl = 0.367; y-int = 43.56
sprintf('Red maj slope is %d and intercept is %d', beta1, beta0)
% sl = 0.891; y-int = 24.17
sprintf('Principle comp. slope is %d and intercept is %d', coeffs_pc(1), coeffs_pc(2))
% sl = 0.757; y-int = 29.11

%% Fitting residuals

% linear regression
lin_res = sim-fit_sims;
stem(obs, lin_res, 'fill','--', 'MarkerFaceColor','r')
title('linear regression residuals')
ylabel('Ozone (ppbv)')
%hold on

% reduced major axis
figure; 
red_res = sim-reduc_sims;
stem(obs, red_res, 'fill','--', 'MarkerFaceColor','g')
title('reduced major axis regression residuals') 
ylabel('Ozone (ppbv)')

% principle component
figure; 
prin_res = sim-pc_sims;
stem(obs, prin_res, 'fill','--', 'MarkerFaceColor','m')
title('principle component regression residuals') 
ylabel('Ozone (ppbv)')
% legend('linear regression', 'reduced major axis', 'principle component')
% title('regression residuals')

% residuals as function of time
plot(days, lin_res, 'r-');
hold on
plot(days, red_res, 'g-');
plot(days, prin_res, 'b-');
plot(days, days.*0, 'k--');
legend('linear regression', 'reduced major axis', 'principle component')
title('Residuals of various regressions as a function of time')
xlabel('Days since June 30')
ylabel('Ozone residuals (ppbv)')
hold off

% chi-squared short cut way 
[lin_reject_that_its_normal prob] = chi2gof(lin_res)
% yes normal, p = 0.5418
[red_reject_that_its_normal prob] = chi2gof(red_res)
% yes normal, p = 0.3938
[prin_reject_that_its_normal prob] = chi2gof(prin_res)
% yes normal, p = 0.2163

% chi-square long way

% If chi2_obs < chi2_mark, the hypothesis cannot be rejected (Yes, probably normal!)

%linear 
lin_bins = 12;
[lin1 centers_linreg] = hist(lin_res, lin_bins); 
lin_created = pdf('norm', centers_linreg, mean(lin_res), std(lin_res));
lin_created_new = lin_created.*(sum(lin1)./sum(lin_created));
lin_num = (lin1 - lin_created_new).^2;
chi2_lin = sum(lin_num./lin_created_new)
chi2_lin_mark = chi2inv(0.95,lin_bins-2-1)
% chi2_obs = 14.1933
% chi2_obs_mark = 16.919

%reduced axis
red_bins = 12;
[red1 centers_red] = hist(red_res, red_bins); 
red_created = pdf('norm', centers_red, mean(red_res), std(red_res));
red_created_new = red_created.*(sum(red1)./sum(red_created));
red_num = (red1 - red_created_new).^2;
chi2_red = sum(red_num./red_created_new)
chi2_red_mark = chi2inv(0.95,red_bins-2-1)
% chi2_obs = 10.3697
% chi2_obs_mark = 16.919


%princ
prin_bins = 12;
[prin1 centers_prinreg] = hist(prin_res, prin_bins); 
prin_created = pdf('norm', centers_prinreg, mean(prin_res), std(prin_res));
prin_created_new = prin_created.*(sum(prin1)./sum(prin_created));
prin_num = (prin1 - prin_created_new).^2;
chi2_prin = sum(prin_num./prin_created_new)
chi2_prin_mark = chi2inv(0.95,prin_bins-2-1)
% chi2_obs = 12.5015
% chi2_obs_mark = 16.919
