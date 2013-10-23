function dT_SS = Bi_SS_Heating(r,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_pump,absorbance,A_tot_powermeter)

if nargin==0 % only used if not called as an external function
    absorbance =0.5; %Aluminum
    %abslayer =40;
    lambda_down=[200 1]; %W/m-K
    C_down=[1.6 1.6]*1e6; %J/m^3-K
    h_down=[10 1e6]*1e-9; %m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)
    eta_down=ones(1,numel(lambda_down)); %isotropic layers, eta=kx/ky;
    
    lambda_up=[0.02]; %W/m-K
    C_up=[0.01]*1e6; %J/m^3-K
    h_up=[1e6]*1e-9; %m ("Al" thickness estimated at 84nm=81nm Al+3nm oxide, from picosecond acoustic)
    eta_up=ones(1,numel(lambda_down)); %isotropic layers, eta=kx/ky;
    
    r_pump=0.4e-6; %probe 1/e^2 radius, m
    r=linspace(0,2*r_pump,20)';
    A_tot_powermeter = 3.5e-3; %laser power (Watts) . . . (assumes the chopper isn't on!)
    t_rep = 12.5e-9;
%    suggested_powers = [A_tot_powermeter/2,A_tot_powermeter/4];
%    p_pulse_heating = (A_tot_powermeter/2)*t_rep/(C(1)*h(1)*2*pi*r^2)
end

%-----------DO NOT MODIFY BELOW HERE------------------------------
f=0; %laser Modulation frequency, Hz
A_abs = absorbance*A_tot_powermeter;
kmin=1/(10000*max(r_pump));
kmax=1/sqrt(r_pump^2)*1.5;
dT_SS  =rombint_multi(@(kvect) BiTDTR_SS_TsKernal(kvect,r,f,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_pump,A_abs),kmin,kmax,length(r));
%-----------------------------------------------------------------