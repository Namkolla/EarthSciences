fname = load('O3_CO.txt');
mos = fname(:,1);
dys = fname(:,2);
o3 = fname(:,3); %O3, ppbv
co = fname(:,4); %CO, ppbv
% both measured in Atlanta

%% PROBLEM 1.1 Jackknife method for mean and std dev of corr coeff

for i = 1 : length(o3)
    % Define two temporary variables for o3 and age
    x = co; 
    y = o3;
    % Eliminate the i-th data point
    x(i) = [];
    y(i) = [];
    % Compute regression line from the n-1 data points
    [r p rlow rhigh] = corrcoef(x,y);
    corrs(i)=r(1,2);
    rlows(i)=rlow(1,2);
    rhighs(i)=rhigh(1,2);
end

mean(corrs) % 0.2631
std(corrs) % 0.0180
figure; hist(corrs,20); title('Resampled correlation coefficient by Jackknife method')

    %is this the right way to get the confidence interval? 
sorted_corrs = sort(corrs); 
y1=prctile(sorted_corrs, 2.5) %0.2255
y2=prctile(sorted_corrs, 97.5) %0.3019

    % Had to assume that the 2.5 and 97.5 percentiles of the sorted vector
    % after the jack-knife loop could give the interval for 95% 

[cov_corr,lags_corr] = xcov(corrs,30,'unbiased');
stem(lags_corr, cov_corr); title('Check for autocorrelation of correlation coefficient using Jackknife method') 
    
%% PROBLEM 1.2 Bootstrap method for mean and std dev of LS slope

p_bootstrp = bootstrp(1000, 'polyfit', co,o3,1);
mean(p_bootstrp(:,1)) % 0.0192
std(p_bootstrp(:,1)) % 0.0095
figure;
hist(p_bootstrp(:,1)); 
title('Resampled regression slopes values by Bootstrap method');
ci = bootci(1000,@polyfit,co,o3,1);
ci(:,1) % conf interval is [-0.0036 0.0347]

[cov_slops,lags_slops] = xcov(p_bootstrp(:,1),20,'unbiased');
stem(lags_slops, cov_slops); title('Check for autocorrelation of least-squares slopes using Bootstrap method') 
    
%% PROBLEM 1.3 TRY: Jackknife method for mean and std dev of LS slope
for i = 1 : length(o3)
    % Define two temporary variables for o3 and age
    y = o3; 
    x = co;
    % Eliminate the i-th data point
    x(i) = [];
    y(i) = [];
    % Compute regression line from the n-1 data points
    [p s] = polyfit(x,y,1);
%     [fittedvals delta] = polyval(p,x,s);
%     upper(i) = mean(fittedvals+2.*delta);
%     lower(i) = mean(fittedvals-2.*delta);
    slopes(i)=p(1,1);
%     deltas(i)=delta;
end

mean(slopes) % 0.0205
std(slopes) % 0.0012
figure; hist(slopes); title('Resampled regression slopes by Jackknife method')

    % calculating 95% conf. interval
sorted_slopes = sort(slopes);
y3=prctile(sorted_slopes, 2.5) % 0.0184
y4=prctile(sorted_slopes, 97.5) % 0.0233

[cov_slops,lags_slops] = xcov(slopes,20,'unbiased');
stem(lags_slops, cov_slops); title('Check for autocorrelation of least-squares slopes using Jackknife method') 

%% PROBLEM 1.4 TRY: Bootstrap method for mean and std dev of corr coeff

p_bootstrp = bootstrp(1000, 'corr', co,o3);
mean(p_bootstrp(:,1)) % 0.2600
std(p_bootstrp(:,1)) % 0.1321
figure;
hist(p_bootstrp(:,1)); 
title('Resampled correlation coefficients by Bootstrap method');
ci = bootci(1000,@corr,co,o3) % conf interval is [0.0206 0.5217]

[cov_slops,lags_slops] = xcov(p_bootstrp(:,1),20,'unbiased');
stem(lags_slops, cov_slops); title('Check for autocorrelation of correlation coefficents using Bootstrap method') 

%% PROBLEM 2.1 - Fitting error for first-order model 

    % "ozone as a function of CO" = ozone is y, CO is x

[p1 s1] = polyfit(co,o3,1);
[fitted1 delta1] = polyval(p1, co, s1); 
%plot(co,o3,'ro')
resids1 = o3-fitted1;
figure; stem(resids1); title('Residuals of 1st order fit of ozone as a function of CO')

    % (1) test for normality of residuals
% bins = round(sqrt(length(data))
bins = 8;
[lin1 centers] = hist(resids1, bins); % lin1 = actual values 
figure; hist(resids1,bins); title('Histogram of residuals of 1st order fit')
created1 = pdf('norm', centers, mean(resids1), std(resids1));
created_new = created1.*(sum(lin1)./sum(created1)); % created1 = synthetic values
lin_num = (lin1 - created_new).^2;
chi2_lin = sum(lin_num./created_new)
chi2_lin_mark = chi2inv(0.95,bins-2-1)
    % chi2_lin = 17.8021
    % chi2_lin_mark = 11.0705
    % Chi2_lin < chi2_lin_mark, so hypothesis cannot be rejected 
    % ie. Yes, probably normal!
subplot(2,1,1); bar(centers, lin1); title('actual values of first order')
subplot(2,1,2); bar(centers, created1); title('theoretical values of first order')

    % (2) test for autocorrelation
[cov1,lags1] = xcov(resids1,30,'unbiased');
stem(lags1, cov1)

    % (3) bootstrap resampling of polynomial coefficients
res1_output=bootstrp(10.^4,'polyfit',co,o3,1);
figure; subplot(2,1,1); hist(res1_output(:,1),round(sqrt(10.^4))); title('distribution of first coefficient (slope), first-order')
subplot(2,1,2); hist(res1_output(:,2),round(sqrt(10.^4))); title('distribution of second coefficient (y-intercept), first-order')


%% PROBLEM 2.2 - Fitting error for third-order model 
[p2 s2] = polyfit(co,o3,3);
[fitted2 delta2] = polyval(p2, co, s2);
resids2 = o3-fitted2; 
figure; stem(resids2); title('Residuals of 3rd order fit of ozone as a function of CO')

    % (1) test for normality of residuals
bins2 = 8;
[lin2 centers2] = hist(resids2, bins2);
figure; hist(resids2,bins2); title('Histogram of residuals of 1st order fit')
created2 = pdf('norm', centers2, mean(resids2), std(resids2));
created_new2 = created2.*(sum(lin2)./sum(created2));
lin_num2 = (lin2 - created_new2).^2;
chi2_lin2 = sum(lin_num2./created_new2)
chi2_lin_mark2 = chi2inv(0.95,bins2-2-1)
    % chi2_lin2 = 8.2996
    % chi2_lin_mark2 = 11.0705
    % Chi2_lin2 < chi2_lin_mark2, so hypothesis cannot be rejected 
    % ie. Yes, probably normal!
subplot(2,1,1); bar(centers, lin2); title('actual values of third order')
subplot(2,1,2); bar(centers, created2); title('theoretical values of third order')
    
    % (2) test for autocorrelation
[cov2,lags2] = xcov(resids2,30,'unbiased');
stem(lags2, cov2)

    % (3) bootstrap resampling of polynomial coefficients
res2_output=bootstrp(10.^4,'polyfit',co,o3,3);
figure; subplot(4,1,1); hist(res2_output(:,1),round(sqrt(10.^4))); title('distribution of first coefficient, third-order')
subplot(4,1,2); hist(res2_output(:,2),round(sqrt(10.^4))); title('distribution of second coefficient, third-order')
subplot(4,1,3); hist(res2_output(:,3),round(sqrt(10.^4))); title('distribution of third coefficient, third-order')
subplot(4,1,4); hist(res2_output(:,4),round(sqrt(10.^4))); title('distribution of fourth coefficient (y-intercept), third-order')