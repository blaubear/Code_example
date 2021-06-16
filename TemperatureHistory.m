function [time1, Temp1,CrFrac1,x1,y1]=TemperatureHistory(fname)
global isimi1
% clear all
% close all
% fname='d:\_Zircon\sim\vert_chamberQ1\post_process\markers\histories\marker.121414.magma.csv';
irock=true;
Tinit=950;
tb=readtable(fname);
A=table2array(tb);
imi=1;
ima=size(A,1);
ifi=ima;
x1=A(ifi,2);
y1=A(ifi,3);
imax=ima;
for i=ima:-1:imi
    if A(i,5) ==0
        imax=i;
    else
        break
    end
    ima=imax;
end
    imi=1;
    for i=1:ima-1
        if A(i,5) == 0 && A(i+1,5) >0 , imi=i+1; end
    end
    imi=max(1,imi-1);
    A(imi,5)=min(0.01,A(imi+1,5));
    A(imi,4)=700.1031927;
if irock
%     x=A(imi:ima,2);
%     y=A(imi:ima,3);
    time=A(imi:ima,1);
    Temp=A(imi:ima,4)+273.15;
    CrFrac=1-A(imi:ima,5);
else
    isimi1=0;
    if imi == 1
%         x=[A(imi,2), A(imi:ima,2)'];
%         y=[A(imi,2), A(imi:ima,3)'];
        time=[A(imi,1)-10, A(imi:ima,1)'];
        Temp= [Tinit,A(imi:ima,4)']+273.15; %Tinit - (Tinit-700)*(time(imi:ima)-time(imi))/(time(ima)-time(imi))+273.15; %
        CrFrac=1-[1,A(imi:ima,5)'];
        isimi1=1;
    else
        imi=imi-1;
%         x=A(imi:ima,2);
%         y=A(imi:ima,3);
        time=A(imi:ima,1);
        Temp=A(imi:ima,4)+273.15;
        CrFrac=1-A(imi:ima,5);
    end
end
is=1;
ima=max(size(Temp));
if ima<3
    time1=-10000;
    Temp1=-10000;
%     x1=x(1);
%     y1=y(1);
    CrFrac1=1;
else
     for i=2:ima
        np=round(abs(Temp(i)-Temp(i-1))/2)+1;
        for j=1:np
            if np>1
%                 xx=x(i-1)+(j-1)*(x(i)-x(i-1))/(np);
%                 yy=y(i-1)+(j-1)*(y(i)-y(i-1))/(np);
                
                TT=Temp(i-1)+(j-1)*(Temp(i)-Temp(i-1))/(np);
                tt=time(i-1)+(j-1)*(time(i)-time(i-1))/(np);
            else
%                 xx=x(i-1);
%                 yy=y(i-1);
                TT=Temp(i-1);
                tt=time(i-1);
            end
            time1(is+j-1)=tt;
            Temp1(is+j-1)=TT;
%             x1(is+j-1)=xx;
%             y1(is+j-1)=yy;
            CrFrac1(is+j-1)=mf_magma(TT);
        end
        is=is+np;
    end
end
end