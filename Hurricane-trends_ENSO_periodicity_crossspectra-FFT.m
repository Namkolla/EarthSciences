%% Hurricane, regression (1) 

hurr_data=open('hurricane_atlantic.dat');
all_data = hurr_data.data;
all_data(:,1)=[];
all_data = all_data';
[x y] = size(all_data);
endval = x.*y;
all_data = reshape(all_data,1,endval);
months = 1948+((1:endval)/12);
plot(months, all_data, 'ro'); hold on
pvals = polyfit(months, all_data,1); 
newyvals = polyval(pvals, months);
plot(months, newyvals, 'k-')
title('Number of hurricanes from 1948 to 2004')
ylabel('Number of hurricanes')
xlabel('Year')
hold off

%% Hurricane, bootstrap (2) and confidence interval (3)  

p_bootstrp = bootstrp(1000, 'polyfit', months,all_data,1);
hurr_mean = mean(p_bootstrp(:,1)); %
hurr_std = std(p_bootstrp(:,1)); % 
figure;
hist(p_bootstrp(:,1),40); 
ylabel('frequency')
title('Resampled regression slopes values by Bootstrap method');
ci(1) = hurr_mean-(2.*hurr_std);
ci(2) = hurr_mean+(2.*hurr_std);
ci % conf interval is [-0.0002 0.0011]

%% Hurricane, detrended periodogram (4) and (5) 
de_all_data = detrend(all_data);
[Pxx3,f3]=periodogram(de_all_data,[],length(de_all_data),12); %sampling frequency is 12
semilogx(1./f3(2:end),Pxx3(2:end)) 
xlabel('period')
ylabel('PSD') 
title('Detrended periodogram')

%% ENSO, FFT (1) 

load('nino12.txt'); 
nino12(:,1) = [];
nino12 = nino12';
[x2 y2] = size(nino12);
endval2 = x2.*y2;
nino12_all = reshape(nino12,1,endval2);

load('nino34.txt'); 
nino34(:,1) = [];
nino34 = nino34';
[x3 y3] = size(nino34);
endval3 = x3.*y3;
nino34_all = reshape(nino34,1,endval3);

nino12_all = detrend(nino12_all);
nino34_all = detrend(nino34_all);
[Pxy4,F4]=cpsd(nino12_all, nino34_all,[],[],length(nino12_all), 12);
semilogx(1./F4(2:end), abs(Pxy4(2:end))); 
xlabel('year');
ylabel('PSD');
title('Cross-spectrum FFT of two data series')

%% ENSO, phase values (2) 

phase = angle(Pxy4)/(2*pi)*360;
lags = interp1(F4,phase,1);
lags=lags/360; % lags = 0.1687

plot(F4,phase,'r-')
xlabel('frequency (years)')
ylabel('phase values')
title('phase values of the cross spectra')

%% ENSO, coherence (3) 

[Cxy,f2] = mscohere(nino12_all,nino34_all,[],[],length(F4),12);
plot(f2,Cxy)
xlabel('frequency (years)')
ylabel('coherance')
title('coherance of cross spectra')
coh = interp1(f2,Cxy,1); % coh = 0.9696
