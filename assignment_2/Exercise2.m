clc
close all
clear all

%% Exercise 2

m=7;
Nb=2^m-1;
N1=0;
N0=0;
Nt=0;
Nnt=0;
numR0=0;
numR1=0;
pnSequence=comm.PNSequence('Polynomial',[7 4 0], 'SamplesPerFrame',Nb,'InitialConditions',[1 1 1 0 1 1 0]);
x1=pnSequence()';

for i=1:Nb
    if x1(i)==0
        v1(i)=-1;
        N0=N0+1;
    elseif x1(i)==1
        v1(i)=1;
        N1=N1+1;
    end        
end
f=figure
h=stairs(0:126,v1)
xlim([0 130]),xticks(0:10:130)
ylim([-2 2]),yticks(-1:1)
set(h,'LineWidth',1.2)
h=title({'5G PSS m-sequence, 127 symbols';'$$N_1=64$$, $$N_0$$=63, $$N_T=64$$, $$N_{NT}=63$$'}), set(h,'Interpreter','Latex')
grid on 
f.Position = [100 100 1000 500];
saveas (figure(1),'signal5G','epsc')

%% Number of transition/no transition

for i=2:Nb+1 %to confront 2 different values I need to start from the second position
    if i==Nb+1
       y=x1(1); %first bit of the sequence because it is cyclic
       if y == x1(i-1)
        Nnt=Nnt+1;
        else
        Nt=Nt+1;
       end
       
    elseif x1(i)== x1(i-1)
        Nnt=Nnt+1;
    else
        Nt=Nt+1;
    end
end

table(N0,N1,Nt,Nnt)

numR0=0; 
numR1=0;
pos0=0;
pos1=0;
n=1;
f=figure, hold on
for i=2:Nb+1
    flag=0;
    if i==Nb+1
       q=1;
       t=Nb;
       flag=1;
       
    else
        q=i;
        t=i-1;
    end
    
    if x1(q)==0
       if x1(q)== x1(t)
           n=n+1;
       else 
           numR0=[numR0;n];
           p=q-n;
           if flag==1
               p=p+t;
           end
           pos0=[pos0,p];
           line([p p],[0 n],'Color','r')
           n=1;
           p=0;
       end
        elseif x1(q)==1
        if x1(q)==x1(t)
            n=n+1;
        else
            numR1=[numR1;n];
            p=q-n;
            if flag==1
               p=p+t;
            end
            pos1=[pos1,p];
            line([p p],[0 n],'Color','b')
            n=1;
            p=0;
        end
         end
end
   
plot(pos0, numR0,'.r'), plot(pos1,numR1,'.b')
hold off
grid on
x=linspace(0,130,14), xticks(x)
xlim([0 130]), xlabel('starting position in the sequence')
ylim([0 8]), ylabel('run length')
title('5G PSS m-sequence, 127 symbols')
f.Position = [100 100 1000 500];
saveas (figure(2),'5G_PSS_m_sequence','epsc')

numR0=numR0';
numR1=numR1';
maximum=max(max(numR0),max(numR0));
NR0=zeros(1,maximum);
NR1=zeros(1,maximum);
index=[1:maximum]';

for i=1:length(numR0)
    for j=1:maximum
        if numR0(i)==j;
            NR0(j)=NR0(j)+1;
        end
    end
end

for i=1:length(numR1)
    for j=1:maximum
        if numR1(i)==j;
            NR1(j)=NR1(j)+1;
        end
    end
end   
NR0=NR0';
NR1=NR1';
i=index;
table(i,NR0,NR1)
 

%% Autocorrelation

x1b=2*x1-1; % bipolar version 0 → -1 1 → +1
R=ifft(fft(x1b).*conj(fft(x1b))); % non-normalized periodic autocorr.
r=length(R);
R=[R(65:r),R(1:64)];
f=figure, plot(-(r-1)/2:1:(r-1)/2,R);
xlim([-80 80]), ylim([-20 140])
grid on
h1=xlabel('$$\tau$$');
h2=ylabel('$$R(\tau)$$');
h=title({'5G PSS m-sequence, 127 symbols, Periodic Autocorrelation';' Max Peak Side Lobe = 1'});
set(h1,'Interpreter','Latex'),set(h2,'Interpreter','Latex');
f.Position = [100 100 1000 500];
saveas (figure(3),'Autocorrelation','epsc');

%% Truncated m sequences

x1_tr=x1(1:117);
x1b_tr=2*x1_tr-1; % bipolar version 0 → -1 1 → +1
R2=ifft(fft(x1b_tr).*conj(fft(x1b_tr))); % non-normalized periodic autocorr.
r2=length(R2);
R2=[R2(60:r2),R2(1:59)];
f=figure, plot(-(r2-1)/2:1:(r2-1)/2,R2);
xlim([-80 80]),ylim([-20 130])
grid on
h1=xlabel('$$\tau$$');
h2=ylabel('$$R(\tau)$$');
h=title({'5G PSS m-sequence, 117 symbols, Periodic Autocorrelation';' Max Peak Side Lobe = 1, starting seed'});
set(h1,'Interpreter','Latex'),set(h2,'Interpreter','Latex');
f.Position = [100 100 1000 500];
saveas (figure(4),'Autocorrelation2','epsc');

R2_t=[R2(60:r2),R2(1:58)];
MPSL=max(abs(R2_t))

%% changing the seed

pnSequence=comm.PNSequence('Polynomial',[7 4 0], 'SamplesPerFrame',Nb,'InitialConditions',[0 1 0 1 0 0 0]);
x2=pnSequence()';

x2_tr=x2(1:117);
x2b_tr=2*x2_tr-1; % bipolar version 0 → -1 1 → +1
R2=ifft(fft(x2b_tr).*conj(fft(x2b_tr))); % non-normalized periodic autocorr.
r2=length(R2);
R2=[R2(60:r2),R2(1:59)];
f=figure, plot(-(r2-1)/2:1:(r2-1)/2,R2);
xlim([-80 80]),ylim([-20 130])
grid on
h1=xlabel('$$\tau$$');
h2=ylabel('$$R(\tau)$$');
h=title({'5G PSS second m-sequence, 117 symbols, Periodic Autocorrelation';' Max Peak Side Lobe = 1, changed seed'});
set(h1,'Interpreter','Latex'),set(h2,'Interpreter','Latex');
f.Position = [100 100 1000 500];
saveas (figure(5),'Autocorrelation3','epsc');
R2_t=[R2(60:r2),R2(1:58)];
MPSL=max(abs(R2_t))