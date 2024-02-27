%% Exercise 3

clc;
clear all;
close all

k=8;
N = 2^k - 1;
V=[];
C=[];

i=1;
for i=1:N
    vector=de2bi(i,k);
    v=fliplr(vector); % it inverts the order of the bits
    V=[V;v];
    vect=v';
    crc_encoder = comm.CRCGenerator('Polynomial','z^8 + z^7 + z^4  + 1');
    codeword=crc_encoder(vect);
    C=[C,codeword];
end 

HW=zeros(1,N);

for j=1:N
    for i=1:length(codeword)
        if C(i,j)== 1 
            HW(j)=HW(j)+1;
        end 
    end
end 

A=zeros(1,length(codeword));
for i=1:length(codeword)
    for j=1:N
        if HW(j)== i
            A(i)= A(i) +1;
        end
    end
end
figure
bar(1:length(A),A, 'b')
hold on 
bar(0,1,'b')
grid on
xlim([0, length(A)])
xlabel('WH')
ylabel('multiplicity')
title("Distance Profile")
saveas (figure(1),'DistProf3','epsc')

n=length(codeword);
p=[4*1e-1, 3*1e-1, 2*1e-1, 1e-1, 5*1e-2, 1e-2, 1e-3, 1e-4, 1e-5];
P_UE=zeros(1,length(p));

for j=1:length(p)
    for i=1:n
        s_p = A(i)*(p(j)^i)*(1-p(j))^(n-i); % single probability
        P_UE(j)= P_UE(j) + s_p;
    end
end 

figure;
loglog(p,P_UE,'-o');
xlim([1e-5 1e-0])
ylim([1e-20 1e0])
ax=gca;
ax.XDir= 'Reverse';
grid on;
xlabel('p')
ylabel('P(UE)')
title('CRC(16,8) Undetected Error Probability')
saveas(figure(2),'CRC(16,8)_Undetected_Error_Probability','epsc');

%% Simulation part

p2=[4e-1, 3e-1, 2e-1, 1e-1];
crc_detector = comm.CRCDetector('Polynomial','z^8 + z^7 + z^4  + 1');

for j=1:length(p2)
    n_UE=0;
    ntr=0;
    while n_UE<100
        for i=1:N
            ntr=ntr+1;
            y=bsc(C(:,i), p2(j)); 
            [~,syndrome]=crc_detector(y);  
            if isequal(syndrome,zeros(1,length(syndrome)))
                if ~isequal(C(:,i), y)
                    n_UE=n_UE+1;
                end
            end
        end 
    end
    P_UE=n_UE/ntr;
    hold on
    loglog(p2(j),P_UE,'r -square');
end 

title('CRC(16,8) Undetected Error Probability with Simulation')
legend('Analytic Curve', 'Simulated Curve')
saveas(figure(2),'CRC(16,8)_Undetected_Error_Probability_With_Simulation','epsc');

