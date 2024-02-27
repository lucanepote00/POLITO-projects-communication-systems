clc;
clear all
close all

%% exercise 2 new
figure;

for i=2:8
    A=GenPolynomials(i);
    subplot(3,3,i-1)
    bar(1:length(A),A, 'b')
    hold on 
    bar(0,1,'b')
    grid on
    xlim([0, length(A)])
    xlabel('WH')
    ylabel('multiplicity')
    title(['k =', num2str(i)])
end
sgtitle('Distance Profiles')
saveas (figure(1),'WH','epsc')
function [A] = GenPolynomials (k)

N = 2^k;
V=[];
C=[];

for i=1:N
    vector=de2bi(i-1,k);
    v=fliplr(vector)'; % it inverts the order of the bits
    V=[V,v];
    crc_encoder = comm.CRCGenerator('Polynomial','z^4 + z^3 + z^2 + 1');
    codeword=crc_encoder(v);
    C=[C, codeword];
end 

C=C';
HW=zeros(1,N);

for i=1:N
    for j=1:length(codeword)
        if C(i,j)== 1 
            HW(i)=HW(i)+1;
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
end