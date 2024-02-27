clc
clear all
close all

%% Hartogs - Hughes algorithm

delta=1;
N=128;
D1=1;
D2=1.01;
D3=1.015;
D4=1.02;
a_1=1;
a_2=0.5;
a_3=0.9;
a_4=0.3;
phi1=2*pi*rand;
phi2=2*pi*rand;
phi3=2*pi*rand;
phi4=2*pi*rand;

for j=1:N
    f=j*delta;
    H_f = a_1*exp(1i*phi1)*exp(-2*pi*1i*f*D1)+a_2*exp(1i*phi2)*exp(-2*pi*1i*f*D2)+a_3*exp(1i*phi3)*exp(-2*pi*1i*f*D3)+a_4*exp(1i*phi4)*exp(-2*pi*1i*f*D4);
    A(j)=1/(abs(H_f))^2;
end

f=figure; subplot(4,1,1)
f.Position = [100 100 1000 700];
plot(1:128,10*log10(A),'-o'), grid on, xlim([0 N]), xlabel('f'), ylabel('att[dB]')
title('Attenuation')
P_TX=median(A)*N;
n_bits=zeros(1,N);
n_tones=0;
flag=0;
A1=A;

while flag==0
    [M,i]=min(A1);
    if M<P_TX
        P_TX=P_TX-M;
        A1(i)=2*A1(i);
        n_bits(i)=n_bits(i)+1;
    else 
        flag=1;
    end
end

for j=1:N
    if n_bits(j)~=0
        n_tones=n_tones+1;
    end
end

N_BITS=sum(n_bits); % total number of transmitted bits
subplot(4,1,2)
plot(1:N,n_bits,'squarer'), grid on
xlim([0 N]), ylim([0 3]), xlabel('f'), ylabel('bit')
title(['Number of active tones:',num2str(n_tones),'; number of transmitted bits:',num2str(N_BITS)])

P_TX=median(A)*N;
ptx=P_TX/n_tones*ones(1,N);
n_bits_2=zeros(1,N);
A2=A;

for i=1:N
    if n_bits(i)~=0
        while A2(i)<ptx(i)
        ptx(i)=ptx(i)-A2(i);
        A2(i)=2*A2(i);
        n_bits_2(i)=n_bits_2(i)+1;
        end
    end
end

N_BITS_2=sum(n_bits_2); % total number of transmitted bits
subplot(4,1,3)
plot(1:N,n_bits_2,'og'), grid on
xlim([0 N]), ylim([0 3]), xlabel('f'), ylabel('bit')
title(['Number of active tones:',num2str(n_tones),'; number of transmitted bits:',num2str(N_BITS_2)])


P_TX=median(A)*N;
ptx=P_TX/N*ones(1,N);
n_bits_3=zeros(1,N);
n_tones_3=0;
A3=A;

for i=1:N
    while A3(i)<ptx(i)
    ptx(i)=ptx(i)-A3(i);
    A3(i)=2*A3(i);
    n_bits_3(i)=n_bits_3(i)+1;
    end
    if n_bits_3(i)~=0
       n_tones_3=n_tones_3+1;
    end
end

N_BITS_3=sum(n_bits_3); % total number of transmitted bits
subplot(4,1,4)
plot(1:N,n_bits_3,'*c'), grid on
xlim([0 N]), ylim([0 3]), xlabel('f'), ylabel('bit')
title(['Number of active tones:',num2str(n_tones_3),'; number of transmitted bits:',num2str(N_BITS_3)])
saveas(figure(1),'case3_Ex3','epsc')


