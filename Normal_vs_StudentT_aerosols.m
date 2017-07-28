
ocarbs = [4.5, 5.9, 6.0, 5.8, 6.4, 5.8, 6.3, 3.5, 4.5, 4.4];
std_ocarbs = std(ocarbs);
mean_ocarbs = mean(ocarbs);
len_ocarbs = length(ocarbs);

% student-t distribution

%X=tinv(P,V) returns the inverse of Student's T cdf with V degrees 
%    of freedom, at the values in P.
t_val1 = abs(tinv(0.025,len_ocarbs-1));
low1 = mean_ocarbs-((t_val1.*std_ocarbs)/sqrt(len_ocarbs-1));
high1 = mean_ocarbs+((t_val1.*std_ocarbs)/sqrt(len_ocarbs-1));
sprintf('The 95%% confidence interval of the measurement using Student''s t-distribution is %d to %d', low1, high1)

% normal distribution

% X = norminv(P,MU,SIGMA) returns the inverse cdf for the normal
%    distribution with mean MU and standard deviation SIGMA, evaluated at
%    the values in P.
z_val1 = norminv(0.025,mean_ocarbs,std_ocarbs);
low2 = mean_ocarbs-((z_val1.*std_ocarbs)/sqrt(len_ocarbs));
high2 = mean_ocarbs+((z_val1.*std_ocarbs)/sqrt(len_ocarbs));
sprintf('The 95%% confidence interval of the measurement using a normal distribution is %d to %d', low2, high2)

%should I be doing this or giving it its own xaxis values like below? 
% xaxis = linspace(-5,5,100); 
% len_xaxis = length(xaxis);

% plot comparing 2 PDFs
xvals = linspace(-10,15,100);
studentt_pdf1 = pdf('t',xvals,len_ocarbs-1);
normal_pdf1 = pdf('norm',xvals,mean_ocarbs,std_ocarbs);
plot(xvals+mean_ocarbs,studentt_pdf1,'b');
hold on
% title('pdf of the student t-distribution')
% xlabel('aerosol organic C (ug m^(-3))')
% ylabel('probability density')
plot(xvals,normal_pdf1, 'r');
% title('pdf of the normal distribution')
xlabel('aerosol organic C (ug m^-^3)')
ylabel('probability density')
legend('student t', 'normal')
axis([-2,12,0,0.45])
title('Normal and Student-t PDFs of aerosol organic carbon emissions in Atlanta')
