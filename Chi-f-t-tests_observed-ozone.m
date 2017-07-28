
%% PART A: Chi-square test

% Null Hyp is that data sampled comes from a normal distribution
% Small p-value = reject the null!! (so, you want a large p-value in this
% case)

file = load('../Atlanta_O3.txt'); % the function "open" doesn't work here; it doesn't actually assign the data to the variable "file"
obs = file(:,3);
sim = file(:,4);
obs(obs<0) = [];
% hist(obs); figure; hist(sim) %looks roughly normal!

bins = 20;
[obs_exp centers_exp] = hist(obs, bins); 
hist(obs, bins); title('Distribution of observed ozone mixing ratios without outlier');
xlabel('Mixing ratios (ppbv, parts per billion by volume)')
ylabel('Occurrences');
obs_created = pdf('norm', centers_exp, mean(obs), std(obs));
obs_created_new = obs_created.*(sum(obs_exp)./sum(obs_created));
obs1_num = (obs_exp - obs_created_new).^2;
chi2_obs = sum(obs1_num./obs_created_new)
chi2_obs_mark = chi2inv(0.95,bins-2-1)
% chi2_obs = 24.6576
% chi2_obs_mark = 27.5871


[sim_exp centers_exp2] = hist(sim, bins); % assumes 10 bins
figure; hist(sim, bins); title('Distribution of simulated ozone mixing ratios without outlier');
xlabel('Mixing ratios (ppbv, parts per billion by volume)')
ylabel('Occurrences')
sim_created = pdf('norm', centers_exp2, mean(sim), std(sim));
sim_created_new = sim_created.*(sum(sim_exp)./sum(sim_created));
sim_num = (sim_exp-sim_created_new).^2;
chi2_sim = sum(sim_num./sim_created_new)
chi2_sim_mark = chi2inv(0.95,bins-2-1)
% chi2_sim = 10.9941
% chi2_sim_mark = 27.5871

% Plot part
figure;
plot(centers_exp,obs_created, 'b');
hold on
plot(centers_exp2, sim_created, 'r')
xlabel('Mixing ratios (ppbv, parts per billion by volume)')
ylabel('Probability distribution')
title('PDF of simulated and observed ozone values in Atlanta') 
legend('observed', 'simulated')
hold off

% If chi2_sim/obs < chi2_mark, the hypothesis cannot be rejected (Yes, probably normal!)
% Thus: obs ozone values ARE likely normal and sim ozone values ARE also
% likely normal 

% Double-check via built-in chi-squared goodness-of-fit test: 
[reject prob] = chi2gof(obs)
% 0 0.1110 (greater than 0.05, normal!)
[reject prob] = chi2gof(sim)
% 0 0.4091 (greater than 0.05, normal!) 


%% PART B: F-test

% Null hyp is that variances are equal 

obs_std = std(obs);
obs_dof = length(obs)-1;
sim_std = std(sim);
sim_dof = length(sim)-1;
f_val = (obs_std^2)./(sim_std^2)
f_crit = finv(0.95,obs_dof, sim_dof)

% If observed (f-val) is less than the ideal value (f_crit), you CANNOT
% reject the null hypothesis that variances are equal (ie. variances are
% probably equal!)

% f-value = 1.2378 and f_crit = 1.5331
% --> Variances are about equal! 




%Student's T-test
[h,p,ci]=ttest2(obs,sim)
% h = 1 -> null hypothesis must be rejected


%% PART C: Introducing the outlier
obs2 = obs;
lastval = length(obs2);
obs2(lastval) = 430;

% Still a normal distribution? 

bins = 20;
[obs2_exp centers_exp] = hist(obs2, bins); figure;
hist(obs2, bins); title('Distribution of observed ozone mixing ratios WITH outlier');
xlabel('Mixing ratios (ppbv, parts per billion by volume)')
ylabel('Occurrences') 
obs2_created = pdf('norm', centers_exp, mean(obs2), std(obs2));
obs2_created_new = obs2_created.*(sum(obs2_exp)./sum(obs2_created));
obs2_num = (obs2_exp - obs2_created_new).^2;
chi2_obs2 = sum(obs2_num./obs2_created_new)
chi2_obs2_mark = chi2inv(0.95,bins-2-1)
% chi2_obs2 = 1.6339e10
% chi2_obs_mark = 27.5871


% Variances still about equal? 

obs2_std = std(obs2);
obs2_dof = length(obs2)-1;
sim_std = std(sim);
sim_dof = length(sim)-1;
f_val2 = (obs2_std^2)./(sim_std^2)
f_crit2 = finv(0.95,obs2_dof, sim_dof)

% f-value = 23.7838 and f_crit = 1.5331
% f-val is now more than f-crit, so variances are not equal anymore! 

%Student's T-test
[h2,p2,ci2]=ttest2(obs2,sim)


% Plot part
figure;
plot(centers_exp,obs2_created, 'b');
hold on
plot(centers_exp2, sim_created, 'r')
xlabel('Mixing ratios (ppbv, parts per billion by volume)')
ylabel('Probability distribution')
title('PDF of simulated and observed ozone values in Atlanta INCLUDING outlier') 
legend('observed', 'simulated')
hold off




%% TRASH CODE 

% e_obs = ones(length(obs),1).*mean(obs)
% chi2_obs = sum(((e_obs-obs).^2)./e_obs)
% chi2_obs_mark = chi2inv(0.95,length(obs)-1)
% 
% e_sim = ones(length(sim),1).*mean(sim)
% chi2_sim = sum(((e_sim-sim).^2)./e_sim)
% chi2_sim_mark = chi2inv(0.95,length(sim)-1)

% Class notes: Chi2 test
% org = load('organicmatter_one.txt');
% [n_exp,v]=hist(org,8);
% n_syn = pdf('norm', v, mean(org), std(org));
% n_syn = n_syn.*sum(n_exp)/sum(n_syn);
% figure
% bar(v,n_syn,'b')
% figure
% bar(v,n_exp,'r')
% chi2 = sum((n_exp-n_syn).^2 ./n_syn);
% % dof = # of bins - # of parameters - 1
% dof = 8-2-1;
% chi2inv(0.95,dof)

%% Jarque-Bera test

%{

The Jarque-Bera test is a goodness-of-fit test of whether sample data have
the skewness and kurtosis matching a normal distribution. In the given 
equation, n is the number of observations (or degrees of freedom in 
general); S is the sample skewness, and K is the sample kurtosis

Q: So is n number of observations or degrees of freedom? Those are 2
different numbers 


file = load('Atlanta_O3.txt'); % the function "open" doesn't work here; it doesn't actually assign the data to the variable "file"
obs = file(:,3);
sim = file(:,4);
obs(obs<0) = [];

obs_sk = skewness(obs);
obs_kur = kurtosis(obs);
obs_n = length(obs);
obs_JB = (obs_n/6).*((obs_sk^2)+(0.25.*((obs_kur-3)^2)))

sim_sk = skewness(sim);
sim_kur = kurtosis(sim);
sim_n = length(sim);
sim_JB = (sim_n/6).*((sim_sk^2)+(0.25.*((sim_kur-3)^2)))

% small value = normally distributed <-- is this correct? What should the
% cut-off be?? Use a reference site and cite in PDF/paper? 

obs_JB =

    5.1616


sim_JB =

    0.0805

%}
