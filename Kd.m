close all
clear all

figure()
pause(0.001)
T = linspace(700+273,900+273,50);
x = 3;
for i = 1:50
Dz(i)    = log(DiffusionCoefficient(T(i),x)/1e4);
Dif(i,:) = TraceDiff(T(i),x);
% kdd(i,:)=Kdd(T(i));
% kz(i)=log(490000/ZrSaturation(T(i)));
end
% plot(T,kdd(:,1),T,kdd(:,2),T,kdd(:,3),T,kdd(:,4),T,kdd(:,5),T,kz,'LineWidth',5);
% legend 'Hf' 'Y' 'U' 'Th' 'P' 'Zr'
%csvwrite('KD.csv',[T', (1000./T)',kz',kdd(:,1:8)]);
csvwrite('D.csv',[T', (1000./T)',Dz',Dif(:,1:8)]);
stop

x1000T=A(:,1);
KD=A(:,6);
k=1;
for i=1:6
    if isnan( KD(i))
    else
        TT(k)=x1000T(i);
        kd(k)=KD(i);
        k=k+1;
    end
end
plot(TT,log(kd),'+r');
hold on
fit1=fit(TT',log(kd)','poly1')
plot(fit1,'--b')
stop


T=linspace(700+273,900+273,50);
x=linspace(1,6,50);
for i=1:50
    for j=1:50
        Dz=log(DiffusionCoefficient(T(i),x(j))/1e4);
        Dif(:)=TraceDiff(T(i),x(j));
         DHF(i,j)=exp(Dif(2)) ;%/Dz; %exp(Dif(1))/exp(Dz);,J)
         DZR(i,j)=exp(Dz);
    end
end
surf(x,1000./T,log10(DHF),'EdgeColor','red');
hold on
surf(x,1000./T,log10(DZR),'EdgeColor','blue');
stop





for i=1:200 %,	950,	1000,	1050];
T(i)=1000/0.9+i*4;
Dz(i)=log(DiffusionCoefficient(T(i),x)/1e4);
Dif(i,:)=TraceDiff(T(i),x);
DDD(i)=exp(Dif(i,1))/exp(Dz(i));
%DHf_z(i,:)=exp(Dif(i,:))*1e4/Dz(i);
end
T1=1000./T;
csvwrite('diffZr.csv',[T1',DDD']);
% semilogy(T-273,DHf_z)
% legend 'Hf' 'Y' 'U' 'Th' 'P'
stop


T1=(T+273);
logaSiO2=0;
logaTiO2=log(1);

logTi = 5.711   - 4800./T1 - logaSiO2 + logaTiO2;
Ti=exp(logTi);
plot(T1-273,Ti)

% Hf=[3476,1829,1100,698,	492];
% figure()
% % pause(0.01)
% plot(T1',log(Hf)','+r');
% hold on
% % fit1=fit(T1',Hf','poly1')
% % plot(fit1,'--b')
% plot(T1,log(20729.61963 * T1 - 15843.31),'--g')
% % KD_Hf=490000/1.34./ZrSaturation(T1);
% % plot(T1,KD_Hf,'--y')
% % plot(T1,89.362659 * T1 - 73.283)

function Csat=ZrSaturation(T)
%Csat=4.414e7/exp(13352/T)/2; %Watson 96, Eq 1, in ppm Zr (divide by 1),or mol Zr (divide by 2) 
%Mfactor = 0.0000048*(T)^2-0.0083626*(T)+4.8484463; % empirical relations from magma
%differentiation calc (file  M_factorsforOleg.xlsx
Mfactor=1.62;
Csat=490000/exp(10108/T+1.16*(Mfactor-1)-1.48); % Boehnkeetal2013ChemGeol351,324 assuming 490000ppm in Zircon
%Csat=490000/(exp(12900/T-0.85*(Mfactor-1)-3.80));
end


function Dif=DiffusionCoefficient(T,x) %X-wt%water
theta=1000/T;
lnD=-(11.4*x+3.13)/(0.84*x+1)-(21.4*x+47)/(1.06*x+1)*theta;%best fit of Zr Diff coefficients (several workers) and WH dependence on XH2O
Dif=exp(lnD)*1e4;% in cm2/s
%Dif=0.1*exp(-235980/8.31/T); %Watson 96, Eq 2 in cm2/s; switch off for our
end

function Dif=TraceDiff(T,x)
% lnD_Hf=(-8.620340372*T-42705.17449-.318918919*x*T+4049.500765*x)/T;
lnD_Hf=log(DiffusionCoefficient(T,x)/3e4);
lnD_Y=-9.925340370-(35587.24428-2615.214492*x)/T;
P=0.01; %Pressure in GPa
lnD_U=-6.37-2.65*sqrt(x) - (44729-1093*P-8944*sqrt(x))/T;
lnD_Th=-7.02-1.30*sqrt(x) - (44682-1370*P-7281*sqrt(x))/T;
lnD_P=((0.9469e1 * x - 0.1815e3) / (x + 0.2374e1) * ((1000 / T) - 0.54e0) - 0.267e2);
Dif=[lnD_Hf,lnD_Y,lnD_U,lnD_Th,lnD_Y,lnD_Y,lnD_Y,lnD_P];
end
 

function kdd=Kdd(T)
X=1000./T;
KD_Hf=11.29 .* X - 2.275;
KD_Y= 19.47 .* X - 13.04;
KD_U= 15.32 .* X - 9.17;
KD_Th=13.02 .* X - 8.354;
Csat=ZrSaturation(T);
KD_Sm=log(13.338*Csat^(-0.622));
KD_Dy=log(2460.0*Csat^(-0.867));
KD_Yb=log(33460.*Csat^(-1.040));
KD_P= 7.646 .* X - 5.047;
kdd=[KD_Hf,KD_Y,KD_U,KD_Th,KD_Sm,KD_Dy,KD_Yb,KD_P];
end
