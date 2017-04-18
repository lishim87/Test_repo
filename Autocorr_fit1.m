%Finite monomer pool with two flagella, with IFT mechanism. This code also
%calculates the autocorrelation and compares it with the theoretical result
%XYZ.

clear all;
close all;

Ntot=500;k=1;
Nift1 = 5;
Nift2 = 5;
v = 2.5;
alpha = 0.01;
r1=((Nift1*v*alpha)/2);
r2=((Nift2*v*alpha)/2);
gamma=0.1;
p1 = zeros(1,Ntot);
p2 = zeros(1,Ntot);
ptot = zeros(1,Ntot);
MaxT=1000000;
t1=1000000;
Lstar=r1*Ntot/(2*r1+gamma);
% eta1=0.002;
% eta2=0.004;


for j=1:1
    j
    
m1=zeros(1,MaxT);
m1(1)=140;
m2=zeros(1,MaxT);
m2(1)=1;
mtot=zeros(1,MaxT);
mtot(1)=m1(1)+m2(1);
monomers=Ntot;
T=zeros(1,MaxT);
T(1)=0;

for i=1:MaxT
    
    k1=r1*(monomers-mtot(1))/m1(i);
    k3=r2*(monomers-mtot(1))/m2(i);
%     k1=r1*(monomers-mtot(1));
%     k3=r2*(monomers-mtot(1));
    
    if (m1(i)<=1)
        k2=0;
    else 
        k2=gamma;
    end
    
    if (m2(i)<=1)
        k4=0;
    else 
        k4=gamma;
    end
        

    
    k0=k1+k2+k3+k4;
    
    % Determine time spent 
    
    CoinFlip1=rand;
    tau(i)=(1/k0)*log(1/CoinFlip1);
    T(i+1)= T(i)+tau(i);
    % Determine reaction
    
    CoinFlip=rand;
    if CoinFlip<=k1/k0
        m1(i+1)=m1(i)+1;
        m2(i+1)=m2(i);
        monomers=monomers-1;
    else
            if (CoinFlip>=k1/k0) && (CoinFlip<=(k1+k2)/k0)
                m1(i+1)=m1(i)-1;
                m2(i+1)=m2(i);
                monomers=monomers+1;
            else
                if (CoinFlip>=(k1+k2)/k0) && (CoinFlip<=(k1+k2+k3)/k0)
                    m1(i+1)=m1(i);
                    m2(i+1)=m2(i)+1;
                    monomers=monomers-1;
                else
                    m1(i+1)=m1(i);
                    m2(i+1)=m2(i)-1;
                    monomers=monomers+1;
                end
            end
    end
    
    mtot(i+1)=m1(i+1)+m2(i+1);
    
  if T(i+1)>=t1
    break;
  end
end
 
 

% figure(1)
% plot(T,m1,'b.')
% hold on
% plot(T,m2,'r.')
% plot(T,mtot,'g.')
% %hold off
% xlabel('time')
% ylabel('filament length')

p1(m1(i+1))= p1(m1(i+1))+1;
p2(m2(i+1))= p2(m2(i+1))+1;
ptot(mtot(i+1))= ptot(mtot(i+1))+1;
 

end

p1 = p1/sum(p1);
p2 = p2/sum(p2);
ptot = ptot/sum(ptot);
x = 1:1:Ntot;

Avg= sum(x.*p1);
variance= sum((x.^2).*p1)-(sum(x.*p1))^2;


% figure(1)
figure;
plot(T(1:i),m1(1:i),'b')
hold on
plot(T(1:i),m2(1:i),'r')
plot(T(1:i),mtot(1:i),'g')
hold off
xlabel('time')
ylabel('filament length')
legend('L1','L2','L1+L2','location', 'SE');

% 
% figure(2)
% plot(T,m1-m2,'b.')
% hold on
% plot(T,m1+m2,'r.')
% hold off
% xlabel('time')
% ylabel('filament difference length')

figure;
t1=T(10000:i);
sig1=m1(10000:i)-m2(10000:i);
sig0=(m1(10000:i)+m2(10000:i))- 2*Lstar;
plot(t1,sig1,'b.',t1,sig0,'r.')
title('Input signals');
xlabel('time');
ylabel('filament difference and sum lengths');
legend('L1-L2','L1+L2-2Lss','location', 'SE');

xq=t1(1):1:t1(end);
sig=interp1(t1,sig1,xq);
sig2=interp1(t1,sig0,xq);
figure;
plot(t1,sig1,'o',xq,sig,':.');
%plot(sig,'b.')
 title('Input signals');
 xlabel('time');
 ylabel('filament difference length');
 legend('signal','Interpolated','location', 'SE');

figure;
plot(t1,sig0,'o',xq,sig2,':.');
%plot(sig,'b.')
 title('Input signals');
 xlabel('time');
 ylabel('filament sum length');
 legend('signal','Interpolated','location', 'SE');

%Arbitrary lags; How do I input time information here?
lagS=size(xq,2)-1;
lagS2=size(xq,2)-1;
%lagS2=1000;

figure;
 autocorr(sig,lagS);
% hold on;
 [h,lag]=autocorr(sig,lagS);
%j1=(r1/Lstar)*(Ntot/Lstar-2);
% %k1=(exp(-j1.*lag))*eta1/j1;
% %k1=(exp(-j1.*lag))/j1;
% %k1=(exp(-j1.*lag));


j1=(r1/Lstar)*(Ntot/Lstar-2);
k1=(exp(-j1.*lag))/j1;
eta1= 1/k1(1);

k1=(exp(-j1.*lag))*eta1/j1;
figure
plot(lag,h,'b');
hold on;
plot(lag,k1,'g');
hold off;
title('auto-correlation l1-l2');
xlabel('time');
ylabel('autocorr');
xlim([0,2000]);
legend('auto(computed)','auto(theory)','location', 'NE')


Eqn = @(a,x)(a/j1)*exp(-j1*x);
startPoints = 1;
[f1,gof,output] = fit(lag',h',Eqn,'Start', startPoints);
figure;
plot(f1,lag',h','predfunc');
xlabel('time');
ylabel('Autocorrelation');
title('Fit with Theorectical function a*exp(-b*x)');
xlim([0,2000]);
eta3=coeffvalues(f1);


figure;
semilogy(lag(1:1000),h(1:1000),lag(1:1000),k1(1:1000))
title('log auto-correlation l1-l2');
xlabel('lag');
ylabel('log autocorr');
legend('auto(computed)','auto(theory)','location', 'SW')
%xlim([0,2000]);

% 
 figure;
 autocorr(sig2,lagS2);
% xlim([0,1500]);
% hold on;
 [h2,lag2]=autocorr(sig2,lagS2);
 j2=(r1*Ntot)/(Lstar^2);
% k2=exp(-j2.*lag2);
% k2=exp(-j2.*lag2)*eta2/j2;
% plot(lag2,k2,'b');
% title('auto-correlation l1+l2');
% hold off;
% xlabel('time');
% ylabel('autocorr');
% 

figure;
k2=exp(-j2.*lag2)/j2;
eta2= 1/k2(1);
k2=exp(-j2.*lag2)*eta2/j2;
plot(lag2,h2,'b');
hold on;
plot(lag2,k2,'g');
hold off;
title('auto-correlation l1+l2');
xlabel('time');
ylabel('autocorr');
xlim([0,2000]);


Eqn = @(a,x)(a/j2)*exp(-j2*x);
startPoints = 1;
[f1,gof,output] = fit(lag2',h2',Eqn,'Start', startPoints);
figure;
plot(f1,lag2',h2','predfunc');
xlabel('time');
ylabel('Autocorrelation');
title('Fit with Theorectical function a*exp(-b*x)');
xlim([0,2000]);
eta4=coeffvalues(f1);

figure;
semilogy(lag2(1:2000),h2(1:2000),lag2(1:2000),k2(1:2000))
title('log auto-correlation l1+l2');
xlabel('lag');
ylabel('log autocorr');
xlim([0,1000]);

