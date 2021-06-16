function  Zircon_Hf_Dif_movBC_sinT( varargin )
%close all

fname='c:\Users\1\Documents\MATLAB\try\Zircon_Mel\zirc_marker.54680.magma.csv';

global  n R A B D F alpha1 beta Xs UCR ZrPl ratio Kd n_trace d_Dif T0
global Dplag Temp tyear XH2O Tsolidus Csupsat V Dscale tscale length Rd W t tfin

[filepath,name,ext]       = fileparts(fname);
fsave                     = ['zirc\zirc_',name,ext];
[time, Temp_arr,CrFrac,xc,yc] = TemperatureHistory(fname);
nt                        = round(max(size(time,1),size(time,2))/1);
Temp_arr =  Temp_arr + 273;

fig_setup;

%Parameters_____________________________________________________
ZirconRadius = 5e-4; %zircon radius in cm, can be changed
tfin         = tyear*100000; 
Temp         =  Temp_arr(1)
XH2O         = 2.; %water content in melt
U            = -0; %arbitrary undersaturation set up by user, - growth, + dissol
Tsolidus     = 500+273;
T0           = 750+273;
Csupsat      = 3; %ppm supersaturation for nucleation upon cooling
%total time in years
%nt           = 500;  %number of timesteps
Dplag        = 1; %Distribution coefficient for main minerals
UCR          = 20000; %Critical concentration for microZircon nucleation
Csat         = ZrSaturation(Temp);
Cmelt        = Csat-U; %melt Zr concentration
CbulkZr      = 100; %bulk rock ppm, set up by user
length       = 100*ZirconRadius/(CbulkZr-Cmelt)^0.333; % computed (spherical geometry only) sphere of influence
length       = 0.2; %initial length of the domain, cm

%END:Parameters_____________________________________________________

%SCALING-----------------
Dscale = DiffusionCoefficient(Temp,XH2O)

n      = 2000;%round(length*2500*8);
tscale = length^2/Dscale;

tfin_d = tfin/tscale;
dt     = tfin_d/(nt-1);
t      = 0;
Xs     = ZirconRadius/length; %initial crystal boundary coordinate
Rd     = 1;
W      = 0;
V      = 0; %V - dissolution rate in cm/s
%END:SCALING-----------------

matrixes; % storage matrixes, separate file


%Initial Conditions 
C0(1:n)    = Cmelt; %U - Zr undersaturation in ppm
n_trace    = 8;
CH0(1:n,1) = 4.7;%Hf  %FCT from Bachmann etal JPet 2002.
CH0(1:n,1) = 6;%Hf  %FCT from Bachmann etal JPet 2002.
CH0(1:n,2) = 21.9; %Y
CH0(1:n,3) = 3.9;%U
CH0(1:n,4) = 12.4; %Th
CH0(1:n,5) = 4.;%Sm;%P
CH0(1:n,6) = 4.;%Dy
CH0(1:n,7) = 4;%Yb
CH0(1:n,8) = 1000;%P
Cint(1,1:n_trace) = 0;
index      = 1;

pause(1e-5)
%MAIn LOOP in time _______________________
for i = 1:nt
  C       = progonka(C0,dt,i,Temp_arr);
  C0      = C;
  CH      = HF_progonka(CH0,dt);
  CH0     = CH;
  CC(i,:) = C(:);
  rr  = R*(Rd-Xs)+Xs;
  set(f,'CurrentAxes',ha1)
  hold off
  plot(rr*1e4*length,C(),'r-') %Linear R scale in um
  %plot (rr*1e4*length,CH(:,1),'r-') %Linear R scale in um
  pause(1e-10);
  t      = t+dt;
  Xs     = max(0.,Xs-V*dt);
  Rd     = max(Xs,Rd+W*dt);
  XXs(i) = Xs*1e4*length; %zircon radius in um
  RRd(i) = Rd*1e4*length; %melt shell radius in um
  VV(i)  = V*length/tscale; %*length/tscale % array of dissolution rate  
  tt(i)  = t/tyear*tscale;
  UU(i)  = C(1); %-C(n)
  Kd     = Kdd(Temp);
  UH(i,1:n_trace) = CH(1,1:n_trace).*Kd(1:n_trace); %-C(n)
  UH(i,3:7)       = UH(i,3:7)*10;
  Tsave(i)        = Temp-273;
  ZrPls(i)        = ZrPl;
  if(i>1), Cint(i,:) = Cint(i-1,:)+4*pi*(UH(i,:)*XXs(i)^2+UH(i-1,:)*XXs(i-1)^2)*(XXs(i)-XXs(i-1))/2; end;
  Temp = Temp_arr(i);
  for p = 1:1:1
      if(V<=0)
        index = index + 1;
        XXs_in(index,1) = Xs*1e4*length;
      else
          index_2 = 2;
          while XXs_in(index_2,1) < Xs*1e4*length 
              index_2 = index_2 +1;
          end
          index_3 = index_2;
          while(index_3 <= n)
          XXs_in(index_3,1)  = 0; 
          for pp = 2:1:n_trace+1
            XXs_in(index_3,pp) = 0;
          end
          index_3 = index_3 + 1;
          end
      XXs_in(index_2,1) = Xs*1e4*length;
      index = index_2;
      nn = 14;
      nnn = 4;
      for pppp = nnn:nn
          C0(pppp) = C0(pppp) + 490000*abs(V*dt/((R(nn)-R(nnn))/(Rd-Xs)))/(nn-nnn+1);
      end
      nn = 2;
      nnn = 1;
      for ppp = 1:1:n_trace
          for pppp = nnn:nn
              CH0(pppp,ppp) = CH0(pppp,ppp) + XXs_in(index-1, ppp+1)*abs(V*dt/((R(nn)-R(nnn))/(Rd-Xs)))/(nn-nnn+1);
          end
      end
      end
      for pp = 2:1:n_trace+1
            XXs_in(index,pp) = UH(i,pp-1);
      end
  end 
  
%  if (Temp<T0), break; end;
end
Vol  = 4*pi*XXs(nt)^3/3;
%csvwrite('integ.csv',[XXs,Cint/Vol]);
isav = 1;
%END:MAIn LOOP in time _______________________
%Saving results
if isav==1
Colheader={'Time','Velocity','Zr Radius','Cell Radius','Temperature','Saturation','Plag content','Zr in Plag','Hf','Y','U','Th','Sm','Dy','Yb','P'};
csvwrite_with_headers('time-evolution.csv',[tt,-VV,XXs,RRd,Tsave,UU,Tsave+273,ZrPls,UH],Colheader)
% fid=fopen('cons.bin','w+');
% fwrite(fid,n,'int');
% fwrite(fid,nt,'int');
% fwrite(fid,R,'single');
% fwrite(fid,tt,'single');
% fwrite(fid,CC,'single');
% fclose(fid);
end
%END:Saving results

%Plotting data ____________________________

set(f,'CurrentAxes',ha)
%figure
hold on
xlim([0 max(RRd)])
ylim([0 tfin/tyear])
zlim([0 max(max(CC))])
for i=1:5:nt
    XXX(1:n) = R(1:n)*(RRd(i)-XXs(i))+XXs(i);
    YYY(1:n) = i*dt*tscale/tyear;
    hold on
    plot3(XXX,YYY,CC(i,:),'b-');
end
plot3(XXs,tt,CC(:,1),'g-');
aaa = zeros(nt,1);
plot3(RRd,tt,aaa,'r-');
plot3(XXs,tt,aaa,'r-');
plot3(RRd,tt,CC(:,n),'g-');
grid on
view(-45,30)
hold off

ylabel('time in years')
xlabel('radial distance, um')
zlabel('Zr concentration')
Colheader={'r','Hf','Y','U','Th'};
fileID = fopen('trace_profiles.txt','w');
fprintf(fileID,'%6s  %12s %12s %12s %12s %12s %12s %12s\n','r','Z','Hf','Y','U','Th','Dy','P');
for k=1:n
    CH(k,:) = CH(k,:)./CH(n,:);
    fprintf(fileID,'%6.2f %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e %12.8e\n',XXX(k),C(k)/C(n),CH(k,1),CH(k,2),CH(k,3),CH(k,4),CH(k,6),CH(k,8));
end
fclose(fileID);


%stop
set(f,'CurrentAxes',ha2)
hold on
plot(tt(1:i-1), XXs(1:i-1),'b-')
set(f,'CurrentAxes',ha3)
hold on
plot(tt(1:i-1), -VV(1:i-1)*dt*tscale,'b-')
set(f,'CurrentAxes',ha4)
hold on
plot(tt(1:i-1), Tsave(1:i-1),'b-')
set(f,'CurrentAxes',ha5)
hold on
kk   = 1;
smax = XXs(i-1);
for k = i-2:-1:2
    if(XXs_in(k,1) ~=0 )
        UHs(kk,:) = XXs_in(k,2:9);
        Xss(kk)   = XXs_in(k,1);
        kk        = kk+1;
        smax      = XXs(k);
    end
end
UHmax    = max(UHs);
sz       = size(UHs);
ind      = 1:n_trace;
ave(ind) = 0;
for k = 2:sz(1)
    ave(ind) = ave(ind)-(UHs(k,ind).*Xss(k)^2+UHs(k-1,ind).*Xss(k-1)^2)*(Xss(k)-Xss(k-1))/2;
end
ave=3*ave/Xss(1)^3;
disp('Hf Y U Th P Sm Dy Yb');
disp(num2str(ave,'%12.0f'));
% ave=(ave*0.014-0.08*[13,2683,84,2.6,1619])/(0.014+0.08);
% disp(num2str(ave,'%6.0f'));
%plot(XXs(1:i-1), UH(1:i-1,1)/UHmax(1),'b-')
for k=1:n_trace-1
plot(Xss(:), UHs(:,k)/UHmax(k),'LineWidth',2)
hold on
end
%plot(Xss(:), UHs(:,2)/UHmax(2),'r-')
%plot(XXs(1:i-1), UH(1:i-1,2)/UHmax(2),'r-')
%xlabel(ha5,['Hf Y U Th P=',num2str(UHmax)])
legend(['Hf= ', num2str(UHmax(1),'%6.0f')], ['Y= ', num2str(UHmax(2),'%6.0f')], ['U= ', num2str(UHmax(3),'%6.0f')], ...
    ['Th= ', num2str(UHmax(4),'%6.0f')],['Sm= ', num2str(UHmax(5),'%6.0f')],...
    ['Dy= ', num2str(UHmax(6),'%6.0f')],['Yb= ', num2str(UHmax(7),'%6.0f')],'Location','eastoutside') %['P= ', num2str(UHmax(5),'%6.0f')],
%END:Plotting data ____________________________


% figure()
% fid=fopen('cons.bin','r');
% R=fread(fid,n+1,'single');
% tt=fread(fid,nt,'single');
% CC=reshape(fread(fid,n*nt,'single'),nt,n);
% fclose(fid);
% mesh(R(1:n),tt,CC)
end
function C=progonka(C0,dt,i,Temp_arr)
global  n R  Rd Dplag ZrPl
global A B D F alpha1 beta Xs XH2O Tsolidus V W Csupsat Dscale UCR 

Temp  = Temp_arr(i);
Dif   = DiffusionCoefficient(Temp,XH2O)/Dscale; %see below Diff Coeff dependednt on water and T in cm2/s
Csat  = ZrSaturation(Temp);
Dflux = Dif*(C0(2) - C0(1))/(R(2)-R(1))/(Rd-Xs); 
dxp   = (Xplag(Temp+1e-5)-Xplag(Temp-1e-5))/2e-5;
%dTdt = (Temp_arr(i+1)-Temp_arr(i))/dt;
dTdt  = 0;
W     = -dxp*dTdt/3/Rd^2;
if W > 0
    a = 0
end
V = Dflux/(C0(1)-490000.); %497871 ppm Zr in pure Zr
if(Xs < 1e-5) && (C0(1)-Csat > Csupsat)
    Xs    = 1.005e-4;
    C0(1) = Csat;
    V     = Dif*Csupsat/Xs/(C0(1)-490000.);
end
if (V <-3.e-3)
    V = -3.e-3;
end
if(Temp<Tsolidus)
Dif = 1.e-10;
end
%boundary conditions for i=1 and n
if(Xs > 1.e-5)
D(1) = 1;
B(1) = 0;
F(1) = Csat;
else
D(1) = 1;
B(1) = -1;
F(1) = 0;
V    = 0;
end
kpl = 0.; 
if (C0(n)-Csat > UCR) && (W < 0.)
    kpl = Dplag;
end 
ZrPl = kpl*C0(n);
D(n) = -Dif-W*(R(n)-R(n-1))*(Rd-Xs)*(1-kpl);
A(n) = Dif;
F(n) = 0; 

%Coefficients for Thomas method
s=Xs;
for i = 2:n-1
psi1 = R(i-1);
psi2 = R(i);
psi3 = R(i+1);
t1   = (Dif * dt);
t5   = (psi1 * Rd - psi1 * s + s) ^ 2;
t6   = psi2 - psi1;
t8   = (t5 / t6);
t12  = (Rd * psi2);
t14  = ((-psi2 + 1) * s + t12) ^ 2;
t15  = Rd - s;
t20  = (-W + V) * psi2 - V;
A(i) = -t14 * t15 * dt * psi2 * t20 - t1 * t8;
t25  = (-psi2 * s + s + t12) ^ 2;
t28  = t25 / (psi3 - psi2);
B(i) = -t1 * t28;
t32  = -t15;
t33  = t32 ^ 2;
t34  = -t6;
t38  = (t32 * psi2 - s) ^ 2;
D(i) = -t1 * (-t28 - t8) - t33 * t34 * t38 - t20 * psi2 * dt * t38 * t32;
t44  = t15 ^ 2;
t48  = (t15 * psi2 + s) ^ 2;
F(i) = -t34 * t44 * t48 * C0(i);
end

%Forward Thomas path
alpha1(2) = -B(1)/D(1);
beta(2)   = F(1)/D(1);
for i = 2:n-1
    alpha1(i+1) = -B(i)/(A(i)*alpha1(i)+D(i));
    beta(i+1)   = (F(i)-A(i)*beta(i))/(A(i)*alpha1(i)+D(i));
end
%Backward Thomas path 
C(n) = (F(n)-A(n)*beta(n))/(D(n)+A(n)*alpha1(n));
for i = n-1:-1:1
    C(i) = C(i+1)*alpha1(i+1)+beta(i+1);
end
refresh;
end
function C=HF_progonka(C0,dt)
global  n R  Rd 
global A B D F alpha1 beta Xs Temp XH2O V W  Dscale  Kd Zmlt n_trace d_Dif
for k = 1:n_trace
    Difu = exp(TraceDiff(Temp,XH2O))/Dscale; %see below Diff Coeff dependednt on water and T in cm2/s
    Dif  = Difu(k)*1e4; 
    if(Xs>1e-5)
        K_d  = Kdd(Temp);
        Kd   = K_d(k);
        D(1) = -Dif+V*(R(2)-R(1))*(Rd-Xs)*(Kd-1);
        B(1) = Dif;
        F(1) = 0;
    else
        D(1) = 1;
        B(1) = -1;
        F(1) = 0;
        V    = 0;
    end
    D(n) = -1;
    A(n) = 1;
    F(n) = 0;
%     W=12*V*Xs^2/Rd^2; %800
%     kgm=[1.3,39,0.85,0.34,2.8]; %800
%    
%     W=12*V*Xs^2/Rd^2; %900
%     kgm=[0.7,17,0.,0.,0.64]; %900
%'Hf Y U Th Sm Dy Yb P'
    kgm  = [19,633,7,19,930,935,393,0]*0.012+[0,13.5,0,0,15.5,17.6,9.6,0]*0.07; %900
    kpl  = kgm(k)*0;
    D(n) = -Dif-W*(R(n)-R(n-1))*(Rd-Xs)*(1-kpl);
    A(n) = Dif;
    
    % D(n)=1;
    % A(n)=0;
    % F(n)=C0(n);
    %Coefficients for Thomas method
    s = Xs;
    for i = 2:n-1
        psi1 = R(i-1);
        psi2 = R(i);
        psi3 = R(i+1);
        t1   = (Dif * dt);
        t5   = (psi1 * Rd - psi1 * s + s) ^ 2;
        t6   = psi2 - psi1;
        t8   = (t5 / t6);
        t12  = (Rd * psi2);
        t14  = ((-psi2 + 1) * s + t12) ^ 2;
        t15  = Rd - s;
        t20  = (-W + V) * psi2 - V;
        A(i) = -t14 * t15 * dt * psi2 * t20 - t1 * t8;
        t25  = (-psi2 * s + s + t12) ^ 2;
        t28  = t25 / (psi3 - psi2);
        B(i) = -t1 * t28;
        t32  = -t15;
        t33  = t32 ^ 2;
        t34  = -t6;
        t38  = (t32 * psi2 - s) ^ 2;
        D(i) = -t1 * (-t28 - t8) - t33 * t34 * t38 - t20 * psi2 * dt * t38 * t32;
        t44  = t15 ^ 2;
        t48  = (t15 * psi2 + s) ^ 2;
        F(i) = -t34 * t44 * t48 * C0(i,k);
    end
    
    %Forward Thomas path
    alpha1(2) = -B(1)/D(1);
    beta(2)   = F(1)/D(1);
    for i = 2:n-1
        alpha1(i+1) = -B(i)/(A(i)*alpha1(i)+D(i));
        beta(i+1)   = (F(i)-A(i)*beta(i))/(A(i)*alpha1(i)+D(i));
    end
    %Backward Thomas path
    C(n,k) = (F(n)-A(n)*beta(n))/(D(n)+A(n)*alpha1(n));
    for i = n-1:-1:1
        C(i,k) = C(i+1,k)*alpha1(i+1)+beta(i+1);
    end
end
end


 
function Dif=DiffusionCoefficient(T,x) %X-wt%water
theta = 1000/T;
lnD   = -(11.4*x+3.13)/(0.84*x+1)-(21.4*x+47)/(1.06*x+1)*theta;%best fit of Zr Diff coefficients (several workers) and WH dependence on XH2O
Dif   = exp(lnD)*1e4;% in cm2/s
%Dif=0.1*exp(-235980/8.31/T); %Watson 96, Eq 2 in cm2/s; switch off for our
end

function Dif=TraceDiff(T,x)
% lnD_Hf=(-8.620340372*T-42705.17449-.318918919*x*T+4049.500765*x)/T;
lnD_Hf = log(DiffusionCoefficient(T,x)/3e4);
lnD_Y  = -9.925340370-(35587.24428-2615.214492*x)/T;
P      = 0.01; %Pressure in GPa
lnD_U  = -6.37-2.65*sqrt(x) - (44729-1093*P-8944*sqrt(x))/T;
lnD_Th = -7.02-1.30*sqrt(x) - (44682-1370*P-7281*sqrt(x))/T;
lnD_P  =((0.9469e1 * x - 0.1815e3) / (x + 0.2374e1) * ((1000 / T) - 0.54e0) - 0.267e2);
Dif    = [lnD_Hf,lnD_Y,lnD_U,lnD_Th,lnD_Y,lnD_Y,lnD_Y,lnD_P];
end

function Csat=ZrSaturation(T)
%Csat=4.414e7/exp(13352/T)/2; %Watson 96, Eq 1, in ppm Zr (divide by 1),or mol Zr (divide by 2) 
%Mfactor = 0.0000048*(T)^2-0.0083626*(T)+4.8484463; % empirical relations from magma
%differentiation calc (file  M_factorsforOleg.xlsx
Mfactor = 1.62;
Csat    = 490000/exp(10108/T+1.16*(Mfactor-1)-1.48); % Boehnkeetal2013ChemGeol351,324 assuming 490000ppm in Zircon
if Csat<10
    Csat = 10
end
%Csat=490000/(exp(12900/T-0.85*(Mfactor-1)-3.80));
end

 %{
function Th=TemperatureHistory(t)
global XH2O
%Th=900+273; %remove for proper T history
global  tscale tyear tfin
%return
tt = 0;
if(t>0) 
     tt = t*tscale;
end
Td   = 840+273; %T of the dike
T0   = 720+273; %T of country rocks
% Th=min(Td,T0+0.1*tt/tyear);
%Th=max(T0,Td-0.5*tt/tyear);
dTdt = (Td-T0)/tfin;
% f=@(x) (1.5+pi*(-sawtooth(x,0.5)-1)/2)/1.5;
%Th   = max(T0,Td+abs(dTdt*tt-5*cos(30*2*pi*tt/tfin)-(Td - T0)/2)-(Td - T0)/2);
%Th   = max(T0,Td-(Td-T0)/2+(Td-T0)/4*cos(-3.14*dTdt*tt/25)+0*cos(-3.14*dTdt*tt));
Th    = max(T0,Td-abs(dTdt*tt+30*cos(5*2*pi*tt/tfin)));
%Th    = max(T0,T0+(Td-T0)/2-(Td-T0)/4*sin(dTdt*tt)-0*cos(30*2*pi*tt/tfin));
%XH2O=3.*(1+0.25*sin(30.*2*pi*tt/tfin));

% 
% %Dyke
% h=20; %half thickness of the dike
% x=0; % coordinate for T evolution
% k=0.6e-6; %temeparature conduction coeff in m2/s
% Th=(Td - T0) * (erf((h - x) * (k * tt)^(-0.5)/2)/2 + erf((h + x) * (k * tt)^(-0.5)/2)/2) + T0;
end
%}

function Xpl=Xplag(T)
x   = T-273;
Xpl = 1.94319815919224E-07*x^4+-6.44583169395759E-04*x^3+0.799105630147523*x^2+-439.159560503417*x+90385.6027471822;
Xpl = Xpl/100;
Xpl = 0;
return
Tsol = 750+273;
Tliq = 900+273;
beta = -(T-Tliq)/(Tliq-Tsol);
xxx  = beta;%tanh(3.*(1.-beta));
Xpl  = min(0.99,max(0.,xxx));
%Xpl=0.0; %- no crystal growth, no motion
end

function kdd=Kdd(T)
X     = 1000./T;
KD_Hf = 11.29 .* X - 2.275;
KD_Y  = 19.47 .* X - 13.04;
KD_U  = 15.32 .* X - 9.17;
KD_Th = 13.02 .* X - 8.354;
Csat  = ZrSaturation(T);
KD_Sm = log(13.338*Csat^(-0.622));
KD_Dy = log(2460.0*Csat^(-0.867));
KD_Yb = log(33460.*Csat^(-1.040));
KD_P  = 7.646 .* X - 5.047;
kdd   = exp([KD_Hf,KD_Y,KD_U,KD_Th,KD_Sm,KD_Dy,KD_Yb,KD_P]);
end

