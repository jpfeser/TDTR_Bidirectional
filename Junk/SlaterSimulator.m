clear  all
close all

absorbance =0.1; %Aluminum
r_pump=0.290e-6; %probe 1/e^2 radius, m
P_laser = 3.5e-3; %laser power (Watts)
Tcrit = 212; %Critical temperature for desorption

%layers below the laser spot (layer 1 is the laser absorbing layer)
lambda_down=[200 1]; %W/m-K
C_down=[1.6 1.6]*1e6; %J/m^3-K
h_down=[10 1e6]*1e-9; %m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)
eta_down=ones(1,numel(lambda_down)); %isotropic layers, eta=kx/ky;

%layers above the laser spot (layer 1 is the one closest to the laser absorbing layer)
lambda_up=[0.02]; %W/m-K
C_up=[0.05]*1e6; %J/m^3-K
h_up=[1e6]*1e-9; %m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)
eta_up=ones(1,numel(lambda_down)); %isotropic layers, eta=kx/ky;

%% temperture profiles for the "base case"
r = linspace(0,r_pump*4,200);
T_profile = Bi_SS_Heating(r,lambda_down,C_down,h_down,eta_down,...
    lambda_up,C_up,h_up,eta_up,...
    r_pump,absorbance,P_laser);

figure(4)
plot(r*1e9,T_profile+25,r*1e9,Tcrit*ones(size(r)));
xlabel('\r (nm)')
ylabel('T (C)')
legend('Surface Temperature','Desorption Temperature')
figure(gcf)
        
%% simulation of Tmax, rcrit vs k_metal, h_metal

X = logspace(0,3,100);
Y = logspace(-9,-6,100);

XX = zeros(length(X),length(Y));
YY = zeros(length(X),length(Y));
Tss_max = zeros(size(XX));
rcrit = zeros(size(XX));

nx = length(X);
ny = length(Y);
parfor i = 1:nx
    %lambda_down(1)=X(i);
    for j = 1:ny
        %h_down(1)=Y(j);
        XX(i,j) = X(i);
        YY(i,j) = Y(j);
        Tss_max(i,j) = Bi_SS_Heating(0,[X(i) lambda_down(2:end)],C_down,[Y(j) h_down(2:end)],eta_down,...
            lambda_up,C_up,h_up,eta_up,...
            r_pump,absorbance,P_laser);
        if Tss_max(i,j)>Tcrit
            rcrit(i,j)=fzero(@(r) (Tcrit-25) - Bi_SS_Heating(r,[X(i) lambda_down(2:end)],C_down,[Y(j) h_down(2:end)],eta_down,...
                lambda_up,C_up,h_up,eta_up,...
                r_pump,absorbance,P_laser),...
                [0 10*r_pump]);
        else
            rcrit(i,j) = -1e-12;
        end
    end
end

%% Tmax contour plot
figure(1)
[cs,h]=contour(XX,YY*1e9,Tss_max);
xlabel('\kappa_{Au} W/m-K')
ylabel('h (nm)')
clabel(cs,h)
set(gca,'Xscale','log')
set(gca,'Yscale','log')

%% rcrit contour plot
figure(2)
[cs,h]=contour(XX,YY*1e9,rcrit*1e9);
xlabel('\kappa_{Au} W/m-K')
ylabel('h (nm)')
set(gca,'Xscale','log')
set(gca,'Yscale','log')


%% rcrit contour plot
figure(3)
[cs,h]=contour(XX,YY*1e9,rcrit*1e9);
xlabel('\kappa (W/m-K)')
ylabel('h (nm)')
axis([0 300 0 20])
labels = 0:100:1000;
clabel(cs,h,labels)