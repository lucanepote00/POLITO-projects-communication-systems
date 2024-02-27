clc
clear all
close all

%% Hard Correction Codeword Error Rate analytic curve

k=4;
n=8;
M=2^k-1;
G = [1 0 0 0 0 1 1 0; 0 1 0 0 1 0 1 1; 0 0 1 0 1 1 0 1; 0 0 0 1 1 1 1 0];
Eb_N0=0:10;
Eb_N0_l=10.^(Eb_N0./10);
l=length(Eb_N0);

Pb=0.5.*erfc(sqrt(k/n*Eb_N0_l));
Pwe=ones(1,l);

for i=0:1
    Pwe = Pwe-nchoosek(n,i).*Pb.^(i).*(1-Pb).^(n-i)
end

semilogy(Eb_N0,Pwe,'b')
grid on
xlim([0 10]), h=xlabel('$$E_{b}$$/$$N_0$$ [dB]'), set(h,'Interpreter','Latex')
ylim([1e-6 1e0]), h=ylabel('$$P_{w}(e)$$'), set(h,'Interpreter','Latex')
title('Codeword Error Rate')

%% LUT

r=n-k;
N=2^r;
e=[zeros(1,n);eye(n)];
s=zeros(N,k);
H=[G(:,5:8);G(:,1:4)];
sigma=sqrt(1./((2*k/n).*Eb_N0_l));
V=[];

s=mod(e*H,2); % add the all zero vector
LUT=[s,e]; % the first 4 columns correspond to the syndrome, the others to the error vectors

%% Codeword Error Rate 

CER=zeros(1,n);

for index=1:7
    Ntc=0;
    Ncw=0;
    while Ncw<500  
            y=zeros(1,n);
            C_R=zeros(1,n);
            vector=randi(2,1,k)-1;
            Ntc=Ntc+1;
            CT=mod(vector*G,2);
            CT_p=(CT*2)-1;
            noise=randn(1,n)*sigma(index); %generate the noise vector
            r_vect=CT_p+noise; %add the noise to the bipolar vector
            for j=1:n
                    if r_vect(j)<= 0
                       y(j)=0;
                    elseif r_vect(j) > 0
                        y(j)=1;
                    end
            end

            s_test=mod(y*H,2); %calculate the syndrome
            flag=0;
            error=zeros(1,n);
            for i=1:9
                    if isequal(s_test,LUT(i,1:4))
                        error=LUT(i,5:12);
                        flag=1;
                    end
            end
            
            if flag==1
                C_R=mod(y+error,2);
                if ~isequal(C_R,CT)
                    Ncw=Ncw+1;
                end
            elseif flag==0
                Ncw=Ncw+1;
            end
    end
    CER(index)=Ncw/Ntc;   
end

hold on 
plot(Eb_N0(1:8),CER,'-ok')

%% Hard Correction vs Soft Correction
for i=1:M
    vector=de2bi(i,k);
    V=[V;vector]; %we save all the vectors into a matrix
end 

C=mod(V*G,2);
C_p=(C*2)-1;
CER=zeros(1,n);

Eb_N0=0:10;
Eb_N0_l=10.^(Eb_N0./10);

HW=zeros(1,M);

for i=1:M
    for j=1:n
        if C(i,j)== 1 
            HW(i)=HW(i)+1;
        end 
    end
end 
dmin=min(HW); % in this case we know to be 3

A=zeros(1,n);
for i=1:n
    for j=1:M
        if HW(j)== i
            A(i)= A(i)+1;
        end
    end      
end

i=1;
while A(i)==0
    i=i+1;
end
Amin=A(i);

Pwe_s_a=0.5*Amin*erfc(sqrt(k/n*dmin*Eb_N0_l)); %error probability for soft correction, asymptotically
hold on, plot(Eb_N0,Pwe_s_a,'--m')

Pwe_s_ub=0; %union bound
for i=1:n
    Pwe_s_ub=Pwe_s_ub + 0.5*A(i)*erfc(sqrt(k/n*i*Eb_N0_l));  
end
hold on, plot(Eb_N0,Pwe_s_ub,'-m')
            
%% Simulation part

Eb_N0=0:7;
Eb_N0_l=10.^(Eb_N0./10);
l=length(Eb_N0);
CER_s=zeros(1,n);
C=[zeros(1,n);C]; %add the all zero vector to the codebook
C_p=(2*C)-1;

for index=1:7
    Ntc=0;
    Ncw=0;
    
    while Ncw<500  
            CR=zeros(1,n);
            dE_2=zeros(1,M+1); %euclidian distance
            vector=randi(2,1,k)-1;
            Ntc=Ntc+1;
            CT=mod(vector*G,2);
            CT_p=(CT*2)-1;
            noise=randn(1,n)*sigma(index); %generate the noise vector
            r_vect=CT_p+noise; %add the noise to the bipolar vector
            for j=1:M+1
                for i=1:n
                    dE_2(j)=dE_2(j)+(r_vect(i)-C_p(j,i))^2;
                end
            end
            dE_2_min=min(dE_2);
            for j=1:M+1
                if dE_2(j)==dE_2_min
                   i_k=j; %correct index
                end
            end

            dE=sqrt(dE_2_min); 
            CR=C(i_k,:);

            if ~isequal(CR,CT)
                Ncw=Ncw+1;
            else 
                Ncw=Ncw;
            end
    end
    CER_s(index)=Ncw/Ntc;   
end

hold on 
plot(Eb_N0,CER_s,'-or')
legend('CER HARD analytic','CER HARD simulation','CER SOFT union bound','CER SOFT asympt. app.', 'CER SOFT simulation', 'Location', 'southwest')
saveas(figure(1),'Correction','epsc')
%% Complete LUT

vector2=zeros(n-1,n);
LUT2=LUT;
z=n;
for ind=0:n-2
    for i=1:n-1
        vector2(i,:)=de2bi(2^ind+2^i,n)
        s2_v=mod(vector2(i,:)*H,2);
        flag=0;
        for j=1:z+1
            if isequal(s2_v,LUT2(j,1:4))
               flag=1;
            end
        end
        if flag==0
           LUT2=[LUT2;s2_v,vector2(i,:)];
           z=z+1;
        end
    end
end

CER=zeros(1,n);

for index=1:7
    Ntc=0;
    Ncw=0;
    while Ncw<500  
            y=zeros(1,n);
            C_R=zeros(1,n);
            vector=randi(2,1,k)-1;
            Ntc=Ntc+1;
            CT=mod(vector*G,2);
            CT_p=(CT*2)-1;
            noise=randn(1,n)*sigma(index); %generate the noise vector
            r_vect=CT_p+noise; %add the noise to the bipolar vector
            for j=1:n
                    if r_vect(j)<= 0
                       y(j)=0;
                    elseif r_vect(j) > 0
                        y(j)=1;
                    end
            end

            s_test=mod(y*H,2); %calculate the syndrome
            flag=0;
            error=zeros(1,n);
            for i=1:M+1
                    if isequal(s_test,LUT2(i,1:4))
                        error=LUT2(i,5:12);
                        flag=1;
                    end
            end
            
            if flag==1
                C_R=mod(y+error,2);
                if ~isequal(C_R,CT)
                    Ncw=Ncw+1;
                end
            elseif flag==0
                Ncw=Ncw+1;
            end
    end
    CER(index)=Ncw/Ntc;   
end

hold on 
plot(Eb_N0(1:8),CER,'Marker','o','Color','#4DBEEE' )
legend('CER HARD analytic','CER HARD simulation','CER SOFT union bound','CER SOFT asympt. app.', 'CER SOFT simulation', 'CER HARD sim. entire LUT', 'Location', 'southwest')
saveas(figure(1),'CorrectionLUT2','epsc')


