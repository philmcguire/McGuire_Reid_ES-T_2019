clc;
clear all;
close all;

%Input Data
t = xlsread('GB1_br_data.csv','A2:A28');
t = t';
tsec = t.*60;
concdata = xlsread('GB1_br_data.csv','B2:B28');
concdata = concdata';

%Make pretty scatter plot
figure;
hold on
scatter(tsec,concdata,'filled')
xlabel('Time (sec)');
ylabel('C/C0');
title('Glass Beads 1 - Bromide Breakthrough Curve');

%Input Model Parameters
c0 = 1; %mg/L
k = 0; %1/s
Dx = 3.8967e-06; %1.5*10^-6.1; %m2/s 3*10^-7  1.5*10^-6.1
v = 6.4933e-05; %m/s 5*10^-5
x = 0.3048; %m

%Minimization
Dxv = [Dx, v];
fun = @(Dxv, tsec) (c0/2)*((erfc((x-Dxv(2).*tsec)./(2.*sqrt(Dxv(1).*tsec))))+exp(Dxv(2)*x/Dxv(1)).*erfc((x+Dxv(2).*tsec)./(2.*sqrt(Dxv(1).*tsec))));
% time = tsec;
% newy = fun(Dxv, time);
% scatter (time, newy)
Dxv0 = [3.8967e-06, 6.4933e-05];%[1.5*10^-6.1, 5.2*10^-5];
Dxv = lsqcurvefit(fun, Dxv0, tsec, concdata);

%New Model Parameters
Dx = Dxv(1)
v = Dxv(2)
%Model
modelt = [0:1:45000];
cmodel = (c0/2)*((erfc((x-v.*modelt)./(2.*sqrt(Dx.*modelt))))+exp(v*x/Dx).*erfc((x+v.*modelt)./(2.*sqrt(Dx.*modelt))));


%Plot Model
modeltmin = modelt./60;
modelthr = modeltmin./60;
plot(modelt,cmodel);
 ylim([0 1.05]);

%Find C/C0 = 0.5
% tsec2 = 6257;
% cmodel = (c0/2)*erfc((x-v*tsec2)/(2*sqrt(Dx*tsec2)))+exp((v*x)/Dx)*erfc((x+v*tsec2)/(2*sqrt(Dx*tsec2)));
% CoverC0 = cmodel/c0;
% brkthrtime = tsec2/60/60;
% brktime = ['Breakthrough time is ', num2str(brkthrtime) ' hours'];
% disp(brktime)