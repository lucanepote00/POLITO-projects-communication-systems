clc;
clear all;
close;

k=4;
n=8;
V=[]; % total matrix

G1 = [1 0 0 0 0 1 1 1; 0 1 0 0 1 0 1 1; 0 0 1 0 1 1 0 1; 0 0 0 1 1 1 1 0];
G2 = [1 0 0 0 0 1 1 0; 0 1 0 0 1 0 1 1; 0 0 1 0 1 1 0 1; 0 0 0 1 1 1 1 0];

N = 2^k - 1;

i=1;
for i=1:N
    vector=de2bi(i,k);
    v=fliplr(vector); % it inverts the order of the bits
    V=[V;v]; %we save all the vectors into a matrix
end 

%%% Calculation of the C matrix

C1=mod(V*G1,2); % binary multiplication of matrices

%% Hamming weight

HW1=zeros(1,N);

for i=1:N
    for j=1:n
        if C1(i,j)== 1 
            HW1(i)=HW1(i)+1;
        end 
    end
end 

%%% Ai

A1=zeros(1,n);
for i=1:8
    for j=1:N
        if HW1(j)== i
            A1(i)= A1(i) +1;
        end
    end
end

%% P(UE)

p=[1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];
P_UE1=zeros(1,length(p));

for j=1:length(p)
    for i=1:n
        s_p = A1(i)*(p(j)^i)*(1-p(j))^(n-i); % single probability
        P_UE1(j)= P_UE1(j) + s_p;
    end
end 

x=logspace(-1,-6, 6);

figure;
loglog(x,P_UE1,'-o');

xlim([10e-7 10e-2])
ylim([10e-21 10e0])
ax=gca;
ax.XDir= 'Reverse';
grid on;

%% Repetition for the second matrix G2

C2=mod(V*G2,2); % binary multiplication of matrices

%% Hamming weight

HW2=zeros(1,N);

for i=1:N
    for j=1:n
        if C2(i,j)== 1 
            HW2(i)=HW2(i)+1;
        end 
    end
end 

%% Ai

A2=zeros(1,n);
for i=1:8
    for j=1:N
        if HW2(j)== i
            A2(i)= A2(i) +1;
        end
    end
end

P_UE2=zeros(1,6);

for j=1:6
    for i=1:n
        s_p = A2(i)*(p(j)^i)*(1-p(j))^(n-i); % single probability
        P_UE2(j)= P_UE2(j) + s_p;
    end
end 

hold on
loglog(x,P_UE2,'-o');

xlabel('p');
ylabel('P(UE)');
title('Undetected Error Probability');

legend('Code1 Analytic', 'Code2 Analytic');
saveas(figure(1),'Undetected_Error_Probability','epsc');

%% Second Part

H1=[G1(1:4, 5:8);eye(4)];
H2=[G2(1:4, 5:8);eye(4)];

p_s = [0.1, 0.09, 0.08, 0.07, 0.06];
P_UE_1=zeros(1,length(p_s));

% H1

for i=1:length(p_s)
    
    n_UE=0;
    n_tr=0; % number of transmitted codewords
    while n_UE < 100
        
        vect=randi(2,1,k)-1;
        n_tr=n_tr+1;
        codeword=mod(vect*G1,2);
        test=1;
        y=bsc(codeword, p_s(i));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
        s=mod(y*H1,2);
        o=zeros(1,4);
        tf=isequal(s,o); %check if the syndrome is zero
        if tf==1
           e=mod(y-codeword,2);
           test=isequal(e,zeros(1,8));
           if test==0
              n_UE = n_UE +1;
           end
        end
    end

    P_UE_1(i)= n_UE/n_tr;
      
end
 hold on
 loglog(p_s,P_UE_1,'m -square');

% H1
P_UE_2=zeros(1,length(p_s));

for i=1:length(p_s)
    
    n_UE=0;
    n_tr=0; % number of transmitted codewords
    while n_UE < 100

        vect=randi(2,1,k)-1;
        n_tr=n_tr+1;
        codeword=mod(vect*G2,2);

        y=bsc(codeword, p_s(i));                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
        s=mod(y*H2,2);
        o=zeros(1,4);
        tf=isequal(s,o); %check if the syndrome is zero
        if tf==1
           e=mod(y-codeword,2);
           test=isequal(e,zeros(1,8));
           if test==0
              n_UE = n_UE +1;
           end
        end
    end

    P_UE_2(i)= n_UE/n_tr;
       
end
hold on
loglog(p_s,P_UE_2,'k -square');

legend('Code1 Analytic', 'Code2 Analytic', 'Code1 Simulation', 'Code2 Simulation');
saveas(figure(1),'Undetected_Error_Probability_With_Simulation','epsc');

figure
subplot(1,2,1)
bar(1:length(A1),A1, 'b')
grid on
xlabel('WH')
ylabel('multiplicity')
title('G1')
subplot(1,2,2)
bar(1:length(A2),A2, 'b')
grid on
xlabel('WH')
ylabel('multiplicity')
title('G2')
sgtitle("Distance Profile")
saveas (figure(2),'DistProf','epsc')
