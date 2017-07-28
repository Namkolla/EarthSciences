% o3 = load(./'Atlanta_O3.txt');
o3_obs = o3(:,3)
o3_sim = o3(:,4)
ind1 = find(o3_obs<0);
o3_sim(ind1) = [];
o3_obs(ind1) = [];
ind2 = find(o3_sim<0);
o3_sim(ind2) = [];
o3_obs(ind2) = [];





%% Part A: Least-squares regression and correlation coefficient

% (1)

% observed data = x-variable
% modeled data as predictor variable = y-variable
[p s] = polyfit(o3_obs,o3_sim,1);
sprintf('The slope is %d',p(1))
sprintf('The y-intercept is %d',p(2))
[o3_sim_fit delta] = polyval(p,o3_obs,s);
err_var = sum(((o3_sim-o3_sim_fit).^2)/(length(o3_obs)-2));
slope_std = sqrt(err_var./var(o3_obs));
tval = abs(tinv(0.025,length(o3_obs)-2));
low = p(1)-(tval.*slope_std);
high = p(1)+(tval.*slope_std);
sprintf('The 95%% confidence interval of the least-squares regression slope is %d to %d', low, high)


% (2) 

[r p rlow rhigh] = corrcoef(o3_obs,o3_sim)
sprintf('The 95%% confidence interval of the correlation coefficient is %d to %d', rlow(1,2), rhigh(1,2))
p(find(p<0.05))

% (3) 

plot(o3_obs, o3_sim, 'go')
hold on
plot(o3_obs, o3_sim_fit, 'b')
plot(o3_obs, o3_sim_fit+(2.*delta), 'r--')
plot(o3_obs, o3_sim_fit-(2.*delta), 'r--')
xlabel('Observed Ozone (ppbv)')
ylabel('Simulated Ozone (ppbv)')
title('Relationship between observed and simulated ozone measurements in Atlanta')
legend('data', 'LS regression line', 'error bounds')




%% Part B: Through-the-origin least-squares regression

% Attempt 1
% Function must be in the form of y = ax+0
% Use the '\' operator

%A = o3_obs(:)\o3_sim_fit(:);
%plot(o3_obs, A*o3_obs,'m');

% Attempt 2
%plot(o3_obs, o3_sim_fit-p(2), 'm')

% Attempt 3 - produces the same thing as attempt 1! 
A = sum(o3_obs.*o3_sim_fit)./sum(o3_obs.^2);
plot(o3_obs, A*o3_obs,'m');

legend('Data', 'Least-squares regression line', 'Upper bound', 'Lower bound', 'Least-squares regression line forced through origin')
axis([0 80 0 110]);
hold off

% Can't figure out how to do (4) ! 





%% Part C: Resampled statistics

[o3_obs_n ind] = sort(o3_obs);
o3_sim_n = o3_sim(ind);
figure

% Subset 1

o3_obs1 = o3_obs_n(1:10);
o3_sim1 = o3_sim_n(1:10);

[p1 s1] = polyfit(o3_obs1,o3_sim1,1); % LS stuff
ls_slope1 = p1(1,1)
[o3_sim_fit1 delta1] = polyval(p1,o3_obs1,s1);
errvar1 = sum((o3_sim1-o3_sim_fit1).^2)/(length(o3_obs1)-2);
xvar1 = var(o3_obs1);
slopestd1 = sqrt(errvar1./xvar1);
tvalue1 = abs(tinv(0.025,length(o3_obs1)-2));
ls_slope_low1 = ls_slope1-tvalue1*slopestd1
ls_slope_high1 = ls_slope1+tvalue1*slopestd1

subplot(2,3,1)
plot(o3_obs1, o3_sim1, 'go') % data points
hold on
plot(o3_obs1, o3_sim_fit1, 'b') % linear regression
plot(o3_obs1, o3_sim_fit1 + 2*delta1, 'b--') % error bars
plot(o3_obs1, o3_sim_fit1 - 2*delta1, 'b--') % error bars
A1 = sum(o3_obs1.*o3_sim1)./sum(o3_obs1.^2)
plot(o3_obs1, A1*o3_obs1,'m'); % through origin
title('Ozone in Atlanta, Subset 1')
xlabel('Observed Ozone (ppbv)')
ylabel('Simulated Ozone (ppbv)')
legend('Data', 'LS', 'Lower Error Bound', 'Upper Error Bound', 'Through-the-origin LS')
axis square
hold off



% Subset 2

o3_obs2 = o3_obs_n(11:20);
o3_sim2 = o3_sim_n(11:20);

[p2 s2] = polyfit(o3_obs2,o3_sim2,1);
ls_slope2 = p2(1,1)
[o3_sim_fit2 delta2] = polyval(p2,o3_obs2,s2)
errvar2 = sum((o3_sim2-o3_sim_fit2).^2)/(length(o3_obs2)-2);
xvar2 = var(o3_obs2);
slopestd2 = sqrt(errvar2./xvar2);
tvalue2 = abs(tinv(0.025,length(o3_obs2)-2));
ls_slope_low2 = ls_slope2-tvalue2*slopestd2
ls_slope_high2 = ls_slope2+tvalue2*slopestd2

subplot(2,3,2)
plot(o3_obs2, o3_sim2, 'go') % data points
hold on
plot(o3_obs2, o3_sim_fit2, 'b') % linear regression
plot(o3_obs2, o3_sim_fit2 + 2*delta2, 'b--') % error bars
plot(o3_obs2, o3_sim_fit2 - 2*delta2, 'b--') % error bars
A2 = sum(o3_obs2.*o3_sim2)./sum(o3_obs2.^2)
plot(o3_obs2, A2*o3_obs2,'m'); % through origin
title('Ozone in Atlanta, Subset 2')
xlabel('Observed Ozone (ppbv)')
ylabel('Simulated Ozone (ppbv)')
axis square
hold off

% Subset 3

o3_obs3 = o3_obs_n(21:30);
o3_sim3 = o3_sim_n(21:30);

[p3 s3] = polyfit(o3_obs3,o3_sim3,1);
ls_slope3 = p3(1,1)
[o3_sim_fit3 delta3] = polyval(p3,o3_obs3,s3)
errvar3 = sum((o3_sim3-o3_sim_fit3).^2)/(length(o3_obs3)-2);
xvar3 = var(o3_obs3);
slopestd3 = sqrt(errvar3./xvar3);
tvalue3 = abs(tinv(0.025,length(o3_obs3)-2));
ls_slope_low3 = ls_slope3-tvalue3*slopestd3
ls_slope_high3 = ls_slope3+tvalue3*slopestd3

subplot(2,3,3)
plot(o3_obs3, o3_sim3, 'go') % data points
hold on
plot(o3_obs3, o3_sim_fit3, 'b') % linear regression
plot(o3_obs3, o3_sim_fit3 + 2*delta3, 'b--') % error bars
plot(o3_obs3, o3_sim_fit3 - 2*delta3, 'b--') % error bars
A3 = sum(o3_obs3.*o3_sim3)./sum(o3_obs3.^2)
plot(o3_obs3, A3*o3_obs3,'m'); % through origin
title('Ozone in Atlanta, Subset 3')
xlabel('Observed Ozone (ppbv)')
ylabel('Simulated Ozone (ppbv)')
axis square
hold off

% Subset 4

o3_obs4 = o3_obs_n(31:40);
o3_sim4 = o3_sim_n(31:40);

[p4 s4] = polyfit(o3_obs4,o3_sim4,1);
ls_slope4 = p4(1,1)
[o3_sim_fit4 delta4] = polyval(p4,o3_obs4,s4)
errvar4 = sum((o3_sim4-o3_sim_fit4).^2)/(length(o3_obs4)-2);
xvar4 = var(o3_obs4);
slopestd4 = sqrt(errvar4./xvar4);
tvalue4 = abs(tinv(0.025,length(o3_obs4)-2));
ls_slope_low4 = ls_slope4-tvalue4*slopestd4
ls_slope_high4 = ls_slope4+tvalue4*slopestd4


subplot(2,3,4)
plot(o3_obs4, o3_sim4, 'go') % data points
hold on
plot(o3_obs4, o3_sim_fit4, 'b') % linear regression
plot(o3_obs4, o3_sim_fit4 + 2*delta4, 'b--') % error bars
plot(o3_obs4, o3_sim_fit4 - 2*delta4, 'b--') % error bars
A4 = sum(o3_obs4.*o3_sim4)./sum(o3_obs4.^2)
plot(o3_obs4, A4*o3_obs4,'m'); % through origin
title('Ozone in Atlanta, Subset 4')
xlabel('Observed Ozone (ppbv)')
ylabel('Simulated Ozone (ppbv)')
axis square
hold off

% Subset 5

o3_obs5 = o3_obs_n(41:50);
o3_sim5 = o3_sim_n(41:50);

[p5 s5] = polyfit(o3_obs5,o3_sim5,1);
ls_slope5 = p5(1,1)
[o3_sim_fit5 delta5] = polyval(p5,o3_obs5,s5)
errvar5 = sum((o3_sim5-o3_sim_fit5).^2)/(length(o3_obs5)-2);
xvar5 = var(o3_obs5);
slopestd5 = sqrt(errvar5./xvar5);
tvalue5 = abs(tinv(0.025,length(o3_obs5)-2));
ls_slope_low5 = ls_slope5-tvalue5*slopestd5
ls_slope_high5 = ls_slope5+tvalue5*slopestd5

subplot(2,3,5)
plot(o3_obs5, o3_sim5, 'go') % data points
hold on
plot(o3_obs5, o3_sim_fit5, 'b') % linear regression
plot(o3_obs5, o3_sim_fit5 + 2*delta5, 'b--') % error bars
plot(o3_obs5, o3_sim_fit5 - 2*delta5, 'b--') % error bars
A5 = sum(o3_obs5.*o3_sim5)./sum(o3_obs5.^2)
plot(o3_obs5, A5*o3_obs5,'m'); % through origin
title('Ozone in Atlanta, Subset 5')
xlabel('Observed Ozone (ppbv)')
ylabel('Simulated Ozone (ppbv)')
axis square
hold off

% Subset 6

o3_obs6 = o3_obs_n(51:60);
o3_sim6 = o3_sim_n(51:60);

[p6 s6] = polyfit(o3_obs6,o3_sim6,1);
ls_slope6 = p6(1,1)
[o3_sim_fit6 delta6] = polyval(p6,o3_obs6,s6)
errvar6 = sum((o3_sim6-o3_sim_fit6).^2)/(length(o3_obs6)-2);
xvar6 = var(o3_obs6);
slopestd6 = sqrt(errvar6./xvar6);
tvalue6 = abs(tinv(0.025,length(o3_obs6)-2));
ls_slope_low6 = ls_slope6-tvalue6*slopestd6
ls_slope_high6 = ls_slope6+tvalue6*slopestd6

subplot(2,3,6)
plot(o3_obs6, o3_sim6, 'go') % data points
hold on
plot(o3_obs6, o3_sim_fit6, 'b') % linear regression
plot(o3_obs6, o3_sim_fit6 + 2*delta6, 'b--') % error bars
plot(o3_obs6, o3_sim_fit6 - 2*delta6, 'b--') % error bars
A6 = sum(o3_obs6.*o3_sim6)./sum(o3_obs6.^2)
plot(o3_obs6, A6*o3_obs6,'m'); % through origin
title('Ozone in Atlanta, Subset 6')
xlabel('Observed Ozone (ppbv)')
ylabel('Simulated Ozone (ppbv)')
axis square
hold off

