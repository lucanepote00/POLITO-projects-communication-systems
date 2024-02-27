clc 
clear all
close all

%% Exercise 2: first scenario
N=8;
Eb_No=0:20;
l=length(Eb_No);
Eb_No_l=10.^(Eb_No./10);
k=randi(2^(N-1));
c_1=de2bi(k,N);
c_1p=2*(c_1)-1;
sigma=sqrt(8./(2*Eb_No_l));
var=sigma.^2;
BER=zeros(1,8);

for index=1:8
    Ntc=0;
    Ncw=0;
    while Ncw<500  
            Ntc=Ntc+1;
            r=zeros(1,N);
            v1=randi(2)-1;
            v1_p=2*v1-1;
            s1=v1_p*c_1p;
            noise=randn(1,N)*sigma(index);
            r=s1+noise;
            dotProduct = dot(r, c_1p);
            test =dotProduct / norm(c_1p)^2;
            if test<0
               proj=0;
            elseif test>0
                proj=1;
            end
            if proj~=v1
               Ncw=Ncw+1;
            end
    end
    BER(index)=Ncw/Ntc;    
end

P_e=0.5*erfc(sqrt(Eb_No_l));

figure, semilogy(0:7,BER,'-ok')
hold on, plot(0:20,P_e,'-b'), hold off
ylim([1e-10 1e0]), xlim([0 20]), grid on 
title('CDMA, single user, N=8')
h=xlabel('$$E_b/N_0$$ [dB]'), set(h,'Interpreter','Latex')
h=ylabel('$$P_b(e)$$'), set(h,'Interpreter','Latex')
legend('simulation','analytical')
saveas(figure(1),'BER1_ex2','epsc')

%% second scenario

k2=randi(2^(N-1));
c_2=de2bi(k2,N);
c_2p=2*(c_2)-1;
BER2=zeros(1,8);

for index=1:8
    Ntc=0;
    Ncw=0;
    while Ncw<500  
            Ntc=Ntc+1;
            r=zeros(1,N);
            v1=randi(2)-1;
            v1_p=2*v1-1;
            s1=v1_p*c_1p;
            v2=randi(2)-1;
            v2_p=2*v2-1;
            s2=v2_p*c_2p;
            noise=randn(1,N)*sigma(index);
            r1=s1+s2+noise;
            dotProduct = dot(r1,c_1p);
            test =dotProduct / norm(c_1p)^2;
            if test<0
               proj=0;
            elseif test>0
                proj=1;
            end
            if proj~=v1
               Ncw=Ncw+1;
            end
    end
    BER2(index)=Ncw/Ntc;   
end

p=abs(dot(c_1p,c_2p));
P_e1=0.5.*erfc(sqrt((8+p)^2./(2*8*var)));
P_e2=0.5.*erfc(sqrt((8-p)^2./(2*8*var)));
an_c=(P_e1+P_e2)/2; %analytic curve

figure, s=semilogy(0:20,P_e,'--r'); s.LineWidth=2;
hold on, plot(0:7,BER,'-or'),plot(0:20,an_c,'-b'),plot(0:7,BER2,'-ok'), plot(0:20,P_e1,'-g',0:20,P_e2,'-m')
ylim([1e-10 1e0]), xlim([0 20]), grid on, xticks(0:20)
h=xlabel('$$E_b/N_0$$ [dB]'); set(h, 'Interpreter', 'Latex')
h=ylabel('$$P_b(e)$$'); set(h,'Interpreter','Latex')
title(['CDMA, 2 users, N = 8, p = ', num2str(p)])
h=legend('2-PAM','sim. single user','analytic','simulation','$$P_{e1}$$','$$P_{e2}$$'); set(h,'Interpreter','Latex')
saveas(figure(2),'BER2_4_ex2','epsc')

figure, subplot(2,1,1)
stairs(c_1p, 'LineWidth',2),ylim([-2 2]), grid on, title('s1')
subplot(2,1,2)
stairs(c_2p, 'LineWidth',2),ylim([-2 2]), grid on, title('s2')
saveas(figure(3),'S_s_4_ex2','epsc')
%% part 4

BER3=zeros(1,8);
an_c_2p=1/8*an_c;
 
for index=1:8
    Ntc=0;
    Ncw=0;
    while Ncw<500  
            Ntc=Ntc+1;
            s=randi(8,1)-1; %random shift
            c_2p_s=[c_2p(N-s+1:N),c_2p(1:N-s)];
            r=zeros(1,N);
            v1=randi(2)-1;
            v1_p=2*v1-1;
            s1=v1_p*c_1p;
            v2=randi(2)-1;
            v2_p=2*v2-1;
            s2=v2_p*c_2p_s;
            noise=randn(1,N)*sigma(index);
            r1=s1+s2+noise;
            dotProduct = dot(r1,c_1p);
            test =dotProduct / norm(c_1p)^2;
            if test<0
               proj=0;
            elseif test>0
                proj=1;
            end
            if proj~=v1
               Ncw=Ncw+1;
            end
    end
    BER3(index)=Ncw/Ntc;   
end

for i=1:7
  c_2p_s=[c_2p(N-i+1:N),c_2p(1:N-i)];
  projection=abs(dot(c_1p,c_2p_s));
  P_e1_2=0.5.*erfc(sqrt((8+projection)^2./(16.*var)));
  P_e2_2=0.5.*erfc(sqrt((8-projection)^2./(16.*var)));
  an_c_2p= an_c_2p+1/8*((P_e1_2+P_e2_2)/2);%analytic curve
end

figure, s=semilogy(0:20,P_e,'--r'); s.LineWidth=2;
hold on, plot(0:20,an_c_2p,'-b'),plot(0:7,BER3,'-ok')
ylim([1e-10 1e0]), xlim([0 20]), grid on, xticks(0:20)
h=xlabel('$$E_b/N_0$$ [dB]'); set(h, 'Interpreter', 'Latex')
h=ylabel('$$P_b(e)$$'); set(h,'Interpreter','Latex')
title('CDMA, 2 users, N = 8 with random cyclic shift')
h=legend('2-PAM','analytic','simulation'); set(h,'Interpreter','Latex')
saveas(figure(4),'BER_4_cyclic','epsc')