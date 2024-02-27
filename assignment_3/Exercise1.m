clc
clear all
close all

%% Exercise 1

m=7;
L=2^m-1;
goldseq = comm.GoldSequence('FirstPolynomial','x^7+x^4+1','SecondPolynomial','x^7+x+1','FirstInitialConditions',[0 0 0 0 0 0 1],'SecondInitialConditions',[0 0 0 0 0 0 1],'Index',1,'SamplesPerFrame',L);
x1 = goldseq()';
x1b=1-2*x1;
R=ifft(fft(x1b).*conj(fft(x1b)));
r=length(R);
R=[R(65:r),R(1:64)];
f=figure, plot(-(r-1)/2:1:(r-1)/2,R);
xlim([-80 80]), ylim([-20 140])
grid on
h1=xlabel('$$\tau$$');
h2=ylabel('$$R(\tau)$$');
h=title('Gold Sequence, $$i$$ positive, 127 symbols, Cyclic Autocorrelation');
set(h1,'Interpreter','Latex'),set(h2,'Interpreter','Latex'),set(h,'Interpreter','Latex');
f.Position = [100 100 1000 500];
saveas(figure(1),'Autocorrelation1_ex1','epsc')

goldseq2 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1','SecondPolynomial','x^7+x+1','FirstInitialConditions',[0 0 0 0 0 0 1],'SecondInitialConditions',[0 0 0 0 0 0 1],'Index',-1,'SamplesPerFrame',L);
x2 = goldseq2()';
x2b=1-2*x2;
R=ifft(fft(x2b).*conj(fft(x2b)));
r=length(R);
R=[R(65:r),R(1:64)];
f=figure, plot(-(r-1)/2:1:(r-1)/2,R);
xlim([-80 80]), ylim([-20 140])
grid on
h1=xlabel('$$\tau$$');
h2=ylabel('$$R(\tau)$$');
h=title('Gold Sequence, $$i$$ negative, 127 symbols, Cyclic Autocorrelation');
set(h1,'Interpreter','Latex'),set(h2,'Interpreter','Latex'),set(h,'Interpreter','Latex');
f.Position = [100 100 1000 500];
saveas(figure(2),'Autocorrelation2_ex1','epsc')
%% cross correlation 

% case 1

i1=randi(L,1);
goldseq1 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1','SecondPolynomial','x^7+x+1','FirstInitialConditions',[0 0 0 0 0 0 1],'SecondInitialConditions',[0 0 0 0 0 0 1],'Index',i1,'SamplesPerFrame',L);
x1 = goldseq1()';
x1b=1-2*x1;
i2=randi(L,1);
goldseq2 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1','SecondPolynomial','x^7+x+1','FirstInitialConditions',[0 0 0 0 0 0 1],'SecondInitialConditions',[0 0 0 0 0 0 1],'Index',i2,'SamplesPerFrame',L);
x2 = goldseq2()';
x2b=1-2*x2;
R=ifft(fft(x1b).*conj(fft(x2b)));
f=figure, plot(R);
xlim([0 130]), ylim([-20 20])
grid on
h1=xlabel('$$\tau$$');
h2=ylabel('$$R(\tau)$$');
h=title('cross-correlation of two Gold Sequences with N=127, i>0');
set(h1,'Interpreter','Latex'),set(h2,'Interpreter','Latex');
saveas(figure(3),'CrossCor1_ex1','epsc')
% case 2

goldseq1 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1','SecondPolynomial','x^7+x+1','FirstInitialConditions',[0 0 0 0 0 0 1],'SecondInitialConditions',[0 0 0 0 0 0 1],'Index',i1,'SamplesPerFrame',L);
x1 = goldseq1()';
x1b=1-2*x1;
goldseq2 = comm.GoldSequence('FirstPolynomial','x^7+x^4+1','SecondPolynomial','x^7+x+1','FirstInitialConditions',[0 0 0 0 0 0 1],'SecondInitialConditions',[0 0 0 0 0 0 1],'Index',-2,'SamplesPerFrame',L);
x2 = goldseq2()';
x2b=1-2*x2;
R=ifft(fft(x1b).*conj(fft(x2b)));
f=figure, plot(R);
xlim([0 130]), ylim([-20 20])
grid on
h1=xlabel('$$\tau$$');
h2=ylabel('$$R(\tau)$$');
h=title('cross-correlation of two Gold Sequences with N=127, i<0');
set(h1,'Interpreter','Latex'),set(h2,'Interpreter','Latex');
saveas(figure(4),'CrossCor2_ex1','epsc')
%% welch bound

N=L;
K=1:200;
WB_1=zeros(1,200);
SB=zeros(1,200);

for k=1:200
    for s=1:10
        arg(k)=1/(k*N-1)*(k*N/nchoosek(N+s-1,s)-1);
   
        if arg(k)<0
            arg(k)=0;

        else
        wb=ceil(N*nthroot(arg(k),2*s));
            if wb>WB_1(k)
               WB_1(k)=wb;
            end
        end
    end
end

for k=1:200
    for s=0:(2*N/5-1)
        sb=ceil(sqrt((2*s+1)*(N-s)+s*(s+1)/2-2^s*N^(2*s+1)/(k*factorial(2*s)*nchoosek(N,s))));
        
        if sb>SB(k)
           SB(k)=sb;
        end
    end
end

WB_2=ceil(N*sqrt((K-1)./(K.*N-1)));
WB_3=ceil(sqrt(N)).*ones(1,200);

f=figure 
plot(K,SB,'-r')
hold on,  plot(K,WB_2,'-c'), plot(K,WB_3,'-k'),plot(K,WB_1,'-ob')
hold off
xlim([0 200]), ylim([0 20])
grid on
h1=xlabel('K');
h2=ylabel('bound');
h=title('N=127');
f.Position = [100 100 1000 500];

%entire set of gold sequences
GC=[];
R_A=zeros(129,127);
for i=-2:L-1
    goldseq = comm.GoldSequence('FirstPolynomial','x^7+x^4+1','SecondPolynomial','x^7+x+1','FirstInitialConditions',[0 0 0 0 0 0 1],'SecondInitialConditions',[0 0 0 0 0 0 1],'Index',i,'SamplesPerFrame',L);
    x1 = goldseq()';
    x1b=1-2*x1;
    GC=[GC;x1b]; 
    R_A(i+3,:)=ifft(fft(x1b).*conj(fft(x1b)));
end
R_A_1=R_A(:,2:127);
r_A=max(max(abs(R_A_1)));
R_C_1=[];

for i=1:129
    for j=(i+1):129
         Rc1=ifft(fft(GC(i,:)).*conj(fft(GC(j,:))));
         R_C_1=[R_C_1;Rc1];
    end
end
r_C=max(max(abs(R_C_1)));
r_M_1=max(r_A,r_C);
hold on, plot(129,r_M_1,'*'), hold off, legend('Sidelnikov binary','Welch original','Welch formula','Welch sqrt(N)','max. cyc. cor.','Location','south')
saveas(figure(5),'Bound1_ex1','epsc')
%% last part

K_1=1:200;
N=126;
SB_1=zeros(1,200);
WB_1_2=zeros(1,200);

for k=1:200
    for s=1:10
        arg(k)=1/(k*N-1)*(k*N/nchoosek(N+s-1,s)-1);
   
        if arg(k)<0
            arg(k)=0;

        else
        wb=ceil(N*nthroot(arg(k),2*s));
            if wb>WB_1_2(k)
               WB_1_2(k)=wb;
            end
        end
    end
end

for k=1:200
    for s=0:(2*N/5-1)
        sb_1=ceil(sqrt((2*s+1)*(N-s)+s*(s+1)/2-2^s*N^(2*s+1)/(k*factorial(2*s)*nchoosek(N,s))));
        
        if sb_1>SB_1(k)
           SB_1(k)=sb_1;
        end
    end
end

WB_2_2=ceil(N*sqrt((K-1)./(K.*N-1)));
WB_3_2=ceil(sqrt(N)).*ones(1,200);

GC_1=GC(:,1:126);
R_A_2=zeros(129,126);

for i=1:129
    R_A_2(i,:)=ifft(fft(GC_1(i,:)).*conj(fft(GC_1(i,:))));
end
R_A_1_2=R_A_2(:,2:126);
r_A_2=max(max(abs(R_A_1_2)));
R_C_1_2=[];

for i=1:129
    for j=(i+1):129
         Rc1=ifft(fft(GC_1(i,:)).*conj(fft(GC_1(j,:))));
         R_C_1_2=[R_C_1_2;Rc1];
    end
end
r_C_2=max(max(abs(R_C_1_2)));
r_M_2=max(r_A_2,r_C_2);

f=figure 
plot(K_1,SB_1,'-r')
hold on, plot(K_1,WB_1_2,'ob'), plot(K_1,WB_2_2,'-c'), plot(K_1,WB_3_2,'-k')
hold off
xlim([0 200])
grid on
h1=xlabel('K');
h2=ylabel('bound');
h=title('N=126');
f.Position = [100 100 1000 500];
hold on, plot(129,r_M_2,'*'), hold off, legend('Sidelnikov binary','Welch original','Welch formula','Welch sqrt(N)','max. cyc. cor.','Location','south')
saveas(figure(6),'Bound2_ex1','epsc')