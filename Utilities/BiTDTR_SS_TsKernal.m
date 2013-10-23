function Ts_Kernal = BiTDTR_SS_TsKernal(kvectin,r,f,lambda_down,C_down,h_down,eta_down,lambda_up,C_up,h_up,eta_up,r_pump,A_pump)
r = r(:); %ensures r is a column vector
kvectin = kvectin(:); %ensures r is a column vector

Nr=length(r);
Nk=length(kvectin);
kmatrix= kvectin*ones(1,Nr); %Nk x Nr matrix
rmatrix = ones(Nk,1)*r';
f = zeros(size(r));

[Integrand_down,G_down]=TDTR_TEMP(kvectin(:)',f(:)',lambda_down,C_down,h_down,eta_down,r_pump,r_pump,A_pump);
[Integrand_up,G_up]=TDTR_TEMP(kvectin(:)',f(:)',lambda_up,C_up,h_up,eta_up,r_pump,r_pump,A_pump);

G_Bi = (G_down.*G_up)./(G_down + G_up);
P = A_pump*exp(-pi^2*kmatrix.^2*r_pump^2/2);
Ts_Kernal = 2*pi*P.*G_Bi.*besselj(0,2*pi*kmatrix.*rmatrix).*kmatrix;