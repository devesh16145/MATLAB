clc; clear variables; close all;

N = 10^5;
n=6;
L=10;
U=2;
alpha = 2;
Rd = 3;
d_user1 = 1000; d_user2 = 500;    %Distances 
beta_i_1 = 0.75; beta_i_2 = 0.25;   %Power allocation factors
eta = 4;                %Path loss exponent
err = [0 10^-3 10^-2 10^-1]; %SIC error
oma_alpha = 0.5 ;
%rayleigh fading coefficient for both users
h1 = sqrt(d_user1^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
h2 = sqrt(d_user2^-eta)*(randn(1,N)+1i*randn(1,N))/sqrt(2);

h1_mod = (abs(h1));
h2_mod = (abs(h2));

Pt = 0:1:20;                %Transmit power in dBm
pt = (10^-3)*10.^(Pt/10);   %Transmit power in linear scale
BW = 10^6;                  %bandwidth
npower = noisepow(BW,1,300);%Noise power
npowerdb = 10*log10(npower);

% snr_db = 10*log10((pt/npower));
p = length(Pt);
p1 = zeros(1,length(Pt));
p2 = zeros(1,length(Pt));
p1_av = zeros(1,length(Pt));
p2_av = zeros(1,length(Pt));
npower_db = log10(npower*10^15);
p1_oma = zeros(1,length(Pt));
p2_oma = zeros(1,length(Pt));

sum_noma = zeros(1,length(Pt));
sum_oma = zeros(1,length(Pt));
%AWGN
w1 = sqrt(npower)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
w2 = sqrt(npower)*(randn(1,N)+1i*randn(1,N))/sqrt(2);
snr_db = Pt/npower_db ;
%Generate random binary data for two users
data1 = randi([0 1],1,N);  
data2 = randi([0 1],1,N); 

S1 = 2*data1 - 1;
S2 = 2*data2 - 1;
rate1 = 0.2; rate2 = 0.5;

for u = 1:p
        
    %Calculate SNRs
    snr_1 = beta_i_1*pt(u)*(h1_mod.^2)./(beta_i_2*pt(u)*(h1_mod.^2)+npower);
    snr_12 = beta_i_1*pt(u)*(h2_mod.^2)./(beta_i_2*pt(u)*(h2_mod.^2)+npower);
    snr_2 = beta_i_2*pt(u)*(h2_mod.^2)/npower;
    snr_1_oma = beta_i_1*pt(u)*(h1_mod.^2)./npower;
   
    %Calculate achievable rates
    R1 = log2(1+snr_1);
    R12 = log2(1+snr_12);
    R2 = log2(1+snr_2);
    sum_noma(u) = mean(R1+R2);
    
    R1_av(u) = mean(R1);
    R12_av(u) = mean(R12);
    R2_av(u) = mean(R2);
    
    R1_oma = oma_alpha*log2(1+snr_1_oma);
    R2_oma = (1-oma_alpha)*log2(1+snr_2);
    sum_oma(u) = mean(R1_oma+R2_oma);
    R1_oma_av(u) = mean(R1_oma);
    R2_oma_av(u) = mean(R2_oma);

    theta_l = cos(((2*n - 1)/2*L)*3.14);
    beta_sum = 0 ;
    for l = 1:L 
        beta_l = 3.14/l * sqrt((1-theta_l^2)) * ((Rd/2)*theta_l + Rd/2) * (1+ ((Rd/2)*theta_l + Rd/2)^alpha); 
        beta_sum = beta_sum + beta_l;
    end
    eta_i = (1/Rd) * beta_sum ;
    tau_1 = factorial(U)/factorial(1-1) * factorial(U - 1);
    tau_2 = factorial(U)/factorial(2-1) * factorial(U - 2);
    p1(u) = tau_1 * eta_i * max(snr_1) ;
    p2(u) = (tau_2/2) * (eta_i^2) * max(snr_2)^2;
    p1_av(u) = mean(p1);
    p2_av(u) = mean(p2);
    
    %Check for outage
    for k = 1:N
        if R1_oma(k) < rate1
            p1_oma(u) = p1_oma(u)+1;
        end
        if (R2_oma(k) < rate2)
            p2_oma(u) = p2_oma(u)+1;
        end
    end
end
pout1_oma = p1_oma/N; 
pout2_oma = p2_oma/N;
for ep = 1:length(err)
    for u= 1:p
         snr_err_2 = beta_i_2*pt(u)*(h2_mod.^2)./(err(ep)*beta_i_1*pt(u)*(h2_mod.^2)+npower);
         R2n = log2(1+snr_err_2);
         R2n_av(u) = mean(R2n) ;
    end
    plot(snr_db, R2n_av, 'linewidth', 2); hold on;
end
grid on;
xlabel('SNR (dB)');
ylabel('spectral efficiency (bps/Hz)');
xlim([0 28])
title('imperfect SIC');
legend('\epsilon = 0 (perfect SIC)', '\epsilon = 10^{-3}','\epsilon = 10^{-2}','\epsilon = 10^{-1}','Location','northwest')

figure;
semilogy(snr_db, sort(p1_av,'descend'),'-r*','linewidth', 2); hold on; grid on;
semilogy(snr_db, sort(p2_av,'descend'),'-go','linewidth', 2);
xlim([0 28])
xlabel('SNR (dB)');
ylabel('Outage probability');
legend('User 1 NOMA(far user)','User 2 NOMA(near user)');

figure;
semilogy(snr_db, pout1_oma,'-b*', 'linewidth', 1.5);hold on; grid on;
semilogy(snr_db, pout2_oma,'-yo', 'linewidth', 1.5);
xlim([0 28])
xlabel('SNR (dB)');
ylabel('Outage probability');
legend('User 1 OMA(far user)','User 2 OMA(near user)');

figure;
plot(snr_db, R1_av,'-r','linewidth', 2); hold on; grid on;
plot(snr_db, R2_av,'-g','linewidth', 2);
plot(snr_db, R1_oma_av,'-b','linewidth', 2);
plot(snr_db, R2_oma_av,'-y','linewidth', 2);
xlim([0 28])
xlabel('SNR (dB)');
ylabel('spectral efficiency (bps/Hz)');
legend('R_1 noma (far)','R_2 noma (near)','R_1 oma (far)','R_2 oma (near)','Location','northwest')

figure;
plot(snr_db, sum_noma,'-r','linewidth', 2); hold on; grid on;
plot(snr_db, sum_oma,'-g','linewidth', 2);
xlim([0 28])
xlabel('SNR (dB)');
ylabel('Sum spectral efficiency(bps/Hz)');
legend('NOMA','OMA')



% figure;
% scatter(R1,R2); hold on; grid on;
% scatter(R1_oma,R2_oma);
% xlabel('user1');
% ylabel('user2');
% legend('NOMA','OMA')