for ii=1:length(lambda_down) %the "down" loop
    %-------------Specific heat-----------------
    C_temp=C_down;
    C_temp(ii)=C_down(ii)*1.01;
    [deltaR_data_C(:,ii),ratio_data_C(:,ii)]=BiTDTR_REFL(tdelay,TCR,tau_rep,f,lambda_down,C_temp,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_pump,r_probe,A_pump);
    delta_C(:,ii)=ratio_data_C(:,ii)-ratio_data;
    Num=log(ratio_data_C(:,ii))-log(ratio_data);
    Denom=log(C_temp(ii))-log(C_down(ii));
    S_C_down(:,ii)=Num/Denom;
    %-------------Thermal Conductivity (ky)----------
    lambda_temp=lambda_down;
    lambda_temp(ii)=lambda_down(ii)*1.01;
    eta_temp = eta_down.*lambda_down./lambda_temp;
    [deltaR_data_L(:,ii),ratio_data_L(:,ii)]=BiTDTR_REFL(tdelay,TCR,tau_rep,f,lambda_temp,C_down,h_down,eta_temp,lambda_up,C_up,h_up,eta_up,r_pump,r_probe,A_pump);%TDTR_REFL(tdelay,TCR,tau_rep,f,lambda_temp,C,h,eta_temp,r_pump,r_probe,A_pump);
    delta_L(:,ii)=ratio_data_L(:,ii)-ratio_data;
    Num=log(ratio_data_L(:,ii))-log(ratio_data);
    Denom=log(lambda_temp(ii))-log(lambda_down(ii));
    S_L_down(:,ii)=Num/Denom;
    %-------------Layer Thickess---------------
    h_temp=h_down;
    h_temp(ii)=h_down(ii)*1.01;
    [deltaR_data_h(:,ii),ratio_data_h(:,ii)]=BiTDTR_REFL(tdelay,TCR,tau_rep,f,lambda_down,C_down,h_temp,eta_down,lambda_up,C_up,h_up,eta_up,r_pump,r_probe,A_pump);
    delta_h(:,ii)=ratio_data_h(:,ii)-ratio_data;
    Num=log(ratio_data_h(:,ii))-log(ratio_data);
    Denom=log(h_temp(ii))-log(h_down(ii));
    S_h_down(:,ii)=Num/Denom;
    %--------------------------------------------
    %-------------Anisotropy---------------
 
    eta_temp(ii) = eta_down(ii)*1.01;
    %eta_temp=eta;
    %eta_temp(ii)=eta_temp(ii)*1.01;
    [deltaR_data_eta(:,ii),ratio_data_eta(:,ii)]=BiTDTR_REFL(tdelay,TCR,tau_rep,f,lambda_down,C_down,h_down,eta_temp,lambda_up,C_up,h_up,eta_up,r_pump,r_probe,A_pump);
    delta_eta(:,ii)=ratio_data_eta(:,ii)-ratio_data;
    Num=log(ratio_data_eta(:,ii))-log(ratio_data);
    Denom=log(eta_temp(ii))-log(eta_down(ii));
    S_eta_down(:,ii)=Num/Denom;
    %--------------------------------------------
end

clear ratio_data_eta ratio_data_h ratio_data_L ratio_data_C

for ii=1:length(lambda_up) %the "up" loop
    %-------------Specific heat-----------------
    C_temp=C_up;
    C_temp(ii)=C_up(ii)*1.01;
    [deltaR_data_C(:,ii),ratio_data_C(:,ii)]=BiTDTR_REFL(tdelay,TCR,tau_rep,f,lambda_down,C_down,h_down,eta_down,lambda_up,C_temp,h_up,eta_up,r_pump,r_probe,A_pump);
    delta_C(:,ii)=ratio_data_C(:,ii)-ratio_data;
    Num=log(ratio_data_C(:,ii))-log(ratio_data);
    Denom=log(C_temp(ii))-log(C_up(ii));
    S_C_up(:,ii)=Num/Denom;
    %-------------Thermal Conductivity (ky)----------
    lambda_temp=lambda_up;
    lambda_temp(ii)=lambda_up(ii)*1.01;
    eta_temp = eta_up.*lambda_up./lambda_temp;
    [deltaR_data_L(:,ii),ratio_data_L(:,ii)]=BiTDTR_REFL(tdelay,TCR,tau_rep,f,lambda_down,C_down,h_down,eta_down,lambda_temp,C_up,h_up,eta_temp,r_pump,r_probe,A_pump);%TDTR_REFL(tdelay,TCR,tau_rep,f,lambda_temp,C,h,eta_temp,r_pump,r_probe,A_pump);
    delta_L(:,ii)=ratio_data_L(:,ii)-ratio_data;
    Num=log(ratio_data_L(:,ii))-log(ratio_data);
    Denom=log(lambda_temp(ii))-log(lambda_up(ii));
    S_L_up(:,ii)=Num/Denom;
    %-------------Layer Thickess---------------
    h_temp=h_up;
    h_temp(ii)=h_up(ii)*1.01;
    [deltaR_data_h(:,ii),ratio_data_h(:,ii)]=BiTDTR_REFL(tdelay,TCR,tau_rep,f,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_temp,eta_up,r_pump,r_probe,A_pump);
    delta_h(:,ii)=ratio_data_h(:,ii)-ratio_data;
    Num=log(ratio_data_h(:,ii))-log(ratio_data);
    Denom=log(h_temp(ii))-log(h_up(ii));
    S_h_up(:,ii)=Num/Denom;
    %--------------------------------------------
    %-------------Anisotropy---------------
 
    eta_temp(ii) = eta_up(ii)*1.01;
    %eta_temp=eta;
    %eta_temp(ii)=eta_temp(ii)*1.01;
    [deltaR_data_eta(:,ii),ratio_data_eta(:,ii)]=BiTDTR_REFL(tdelay,TCR,tau_rep,f,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_temp,r_pump,r_probe,A_pump);
    delta_eta(:,ii)=ratio_data_eta(:,ii)-ratio_data;
    Num=log(ratio_data_eta(:,ii))-log(ratio_data);
    Denom=log(eta_temp(ii))-log(eta_up(ii));
    S_eta_up(:,ii)=Num/Denom;
    %--------------------------------------------
end

r_pumptemp=r_pump*1.01;
r_probetemp=r_probe*1.01
[deltaR_data_r_pump,ratio_data_r_pump]=BiTDTR_REFL(tdelay,TCR,tau_rep,f,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_pumptemp,r_probetemp,A_pump);
delta_r_pump(:,1)=ratio_data_r_pump-ratio_data;
Num=log(ratio_data_r_pump)-log(ratio_data);
Denom=log(r_pumptemp)-log(r_pump);
S_r_pump=Num/Denom;

close all
figure(1)
semilogx(tdelay,[S_C_down],'*','MarkerSize',8)
hold on
semilogx(tdelay,[S_C_up],'-*','MarkerSize',8)
semilogx(tdelay,[S_L_down],'o','MarkerSize',8)
semilogx(tdelay,[S_L_up],'-o','MarkerSize',8)
semilogx(tdelay,[S_h_down],'x','MarkerSize',8)
semilogx(tdelay,[S_h_up],'-x','MarkerSize',8)
semilogx(tdelay,[S_r_pump],'-','LineWidth',2)
semilogx(tdelay,[S_eta_down],'+','MarkerSize',8)
semilogx(tdelay,[S_eta_up],'-+','MarkerSize',8)
Cplotlab_down=strcat('C_D',int2str((1:length(lambda_down))'));
Lplotlab_down=strcat('kzD',int2str((1:length(lambda_down))'));
tplotlab_down=strcat('h_D',int2str((1:length(lambda_down))'));
etaplotlab_down=strcat('kxD',int2str((1:length(lambda_down))'));
Cplotlab_up=strcat('C_U',int2str((1:length(lambda_up))'));
Lplotlab_up=strcat('kzU',int2str((1:length(lambda_up))'));
tplotlab_up=strcat('h_U',int2str((1:length(lambda_up))'));
etaplotlab_up=strcat('kxU',int2str((1:length(lambda_up))'));
%fplotlab='f__';
legend([Cplotlab_down;Cplotlab_up;Lplotlab_down;Lplotlab_up;tplotlab_down;tplotlab_up;'R   ';etaplotlab_down;etaplotlab_up])
set(gca,'FontSize',16)
xlabel('td (ps)','Fontsize',16)
ylabel('Ratio Sensitivity (ps)','FontSize',16)

%dlmwrite('sensitivity.txt',[tdelay,S_C,S_L,S_t,S_r_pump,S_r_probe,S_eta]);


% figure(2)
% semilogx(tdelay,ratio_data,'o')
% axis([100e-12 10e-9 0 max(ratio_data)*1.3])
% title('-Vin/Vout')
% 
% figure(3)
% semilogx(tdelay,real(deltaR_data)/TCR,'b-',tdelay,imag(deltaR_data)/TCR,'g-')
% ampl=sqrt(conj(deltaR_data).*deltaR_data);
% axis([100e-12 10e-9 -max(ampl)*1.3/TCR max(ampl)*1.3/TCR])
% title('Vin')
toc
fprintf('Sensitivities calculated\n')