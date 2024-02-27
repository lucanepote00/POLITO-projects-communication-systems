clc
clear all
close all

%% Exercise 3

N=16;
SNR=-10; %in dB

p=randi(2,1,N)-1;
p_p=2*p-1;
Ns=100000;
Ne=20;
sigma=1/sqrt((10^(SNR/10)));
T0=zeros(1,Ns); T1=zeros(1,N);

for i=1:Ns
    noise=randn(1,N)*sigma;
    r0=noise;
    r1=p_p+noise;
    T0(i)=r0*p_p';
    T1(i)=r1*p_p';
end

T0_m=mean(T0)*ones(1,Ns);
T1_m=mean(T1)*ones(1,Ns);
figure, plot(1:Ns,T0,'r'), hold on, plot(1:Ns,T0_m, 'w');
ylim([-70 70]);
xlabel('simulation number')
grid on
h=ylabel('$$T_0$$')
set(h,'Interpreter','Latex')
h=title({'soft correlation: simulation outcomes';'SNR=-10.00 dB'})
set(h,'Interpreter','Latex')
saveas(figure(1),'T0','epsc')
figure, plot(1:Ns,T1,'b'), hold on, plot(1:Ns,T1_m, 'w');
ylim([-70 70]);
xlabel('simulation number')
grid on
h=ylabel('$$T_1$$')
set(h,'Interpreter','Latex')
h=title({'soft correlation: simulation outcomes';'SNR=-10.00 dB'})
set(h,'Interpreter','Latex')
saveas(figure(2),'T1','epsc')

figure, h=histogram(T0,'Normalization','probability'), h.FaceColor=[1 0 0];
hold on, h=histogram(T1,'Normalization','probability'),h.FaceColor=[0 0 1];
h=title({'soft correlation: pdf of T under $$H_0$$ and $$H_1$$';'SNR=-10.00 dB'})
set(h,'Interpreter','Latex')
saveas(figure(3),'PDF','epsc')

%% False Alarm and missed detection probabilities 

t=[-50:70];
l=length(t);
P_fa=zeros(1,l);
P_md=zeros(1,l);
t0=zeros(1,l);
t1=zeros(1,l);

for j=1:l
    N_fa=0; %number of false alarms
    N_md=0; %number of missed detections
    for i=1:Ns
            if T0(i)>t(j)
                N_fa=N_fa+1;
            end
    end
    if N_fa>=Ne
       P_fa(j)=N_fa/Ns;
    end
    
    for i=1:Ns
        if T1(i)<t(j)
            N_md=N_md+1;
        end
    end
    if N_md>=Ne
       P_md(j)=N_md/Ns;
    end   
end

figure, semilogy(t,P_fa,'or')
hold on 
plot(t,P_md,'ob')
grid on
xlim([-50 70])
ylim([1e-8 1e0])
xlabel('threshold t')
h=ylabel('P$$_{fa}$$, P$$_{md}$$'), set(h,'Interpreter','Latex')
h=title({'soft correlation: P$$_{fa}$$ and P$$_{md}$$ vs. threshold';'SNR=-10.00 dB'}),set(h,'Interpreter','Latex')
h=legend('P$$_{fa,s}$$','P$$_{md,s}$$'), set(h,'Interpreter','Latex','Location','southeast')
saveas(figure(4),'SoftCorrelation1','epsc')

%% ROC curve

P_d=(1-P_md);
figure, h=loglog(P_fa,P_d,'ok')
set(h,'LineWidth',1.2)
xlim([1e-8 1e0])
h=xlabel('P$$_{fa}$$')
set(h,'Interpreter','Latex')
ylim([1e-8 1e0])
h=ylabel('P$$_{d}$$')
set(h,'Interpreter','Latex')
legend('simulation','Location','southeast')
grid on
title({'ROC curve';'SNR=-10.00 dB'})
saveas(figure(5),'ROC','epsc')

%% Analytical curves

P_fa_a=1/2*erfc(t./sqrt((2*N*sigma^2)));
P_md_a=1-1/2*erfc((t-N)./sqrt((2*N*sigma^2))); 

figure
semilogy(t,P_fa,'or')
hold on, plot(t,P_md,'ob')
grid on
xlim([-50 70])
ylim([1e-8 1e0])
xlabel('threshold t')
plot(t,P_fa_a,'-m'), plot(t,P_md_a,'-k')
h=ylabel('$$P_{fa}$$, $$P_{md}$$'), set(h,'Interpreter','Latex')
h=title({'soft correlation: $$P_{fa}$$ and $$P_{md}$$ vs. threshold';'SNR=-10.00 dB'}),set(h,'Interpreter','Latex')
h=legend('$$P_{fa,s}$$','$$P_{md,s}$$','$$P_{fa,a}$$','$$P_{md,a}$$'), set(h,'Interpreter','Latex')
hold off
saveas(figure(6),'SoftCorrelation2','epsc')

%% Analytic ROC curve

P_d_a=1-P_md_a;

figure, h=loglog(P_fa,P_d,'ok')
set(h,'LineWidth',1.2)
xlim([1e-8 1e0])
h=xlabel('P$$_{fa}$$')
set(h,'Interpreter','Latex')
ylim([1e-8 1e0])
h=ylabel('P$$_{d}$$')
set(h,'Interpreter','Latex')
hold on, h=loglog(P_fa_a,P_d_a,'--k')
set(h,'LineWidth',1.2)
legend('simulation','analytic','Location','southeast')
grid on
hold off
title({'ROC curve';'SNR=-10.00 dB'})
saveas(figure(7),'ROC2','epsc')

%% P_md vs SNR

N=16;
SNR_2=-10:8; % in dB
P_fa_A=[1e-4,1e-2];
sigma=1./sqrt((10.^(SNR_2/10)));
l=length(sigma);
P_md_A=zeros(2,l);
thr=zeros(2,l);

for j=1:2
    thr(j,:)=(erfcinv(2*P_fa_A(j)))*sqrt(2*N*sigma.^2);
    P_md_A(j,:)=1-1/2*erfc((thr(j,:)-N)./sqrt(2*N*sigma.^2));
end

figure, semilogy(SNR_2,P_md_A(1,:),'-ob')
hold on, semilogy(SNR_2,P_md_A(2,:),'-or')
grid on
xlim([-10 8]), ylim([1e-10 1e0])
xlabel('SNR [dB]'), h=ylabel('P$$_{md}$$')
set(h,'Interpreter','Latex')
h=legend('P$$_{fa}$$=1e-4','P$$_{fa}$$=1e-2','Location','southwest'), set(h,'Interpreter','Latex')
h=title('P$$_{md}$$ for P$$_{fa}$$=1e-4 and P$$_{fa}$$=1e-2'), set(h,'Interpreter','Latex')
saveas(figure(8),'PFA_PMD','epsc')