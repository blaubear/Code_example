tic
clc, clear
cef_end=10;

N_cryst    = 1000;              %верхняя граница колличества кристаллов;
nn         = 1;
Crystal_size_global = zeros(20,N_cryst);
str_Speciment = ('HRT-A');
%str_Speciment = ('Yellowstone YL-4');
for cef=20:20
   
    Spec           = 1;
    FLAG_log_graph = 1;
    
    if(Spec==1)
T_start   = 897+273;            %начальная температура в кельвинах;
T_end     = 888.8804+273;            %конечная температура в кельвинах;
M         = 1.2607;
t_end_yrs = 10000;               %колличество лет;
str_Speciment = ('Yellowstone HRT-C');

X_H20          = 3;             %массовое содержание воды в расплаве;

L              = 2;

J_0            = L*t_end_yrs*1e+6;
A              = 9.e+7;
B              = 8.5e+3;

Coef_delta_C   = 3;
Coef_N         = 3.95;
Coef_Nucl_mult = L*t_end_yrs*9e-3;        %коэффициент в формуле для скорости нуклеации перед экспоентами; при увелчении значения график скорости нуклеации едет влево и график распределения болле пологий
Coef_nucl_c    = 0.06e-2;       %коэффициент в формуле для скорости нуклеации в экспненте с разницей концентраций; при увелчении значения график скорости нуклеации едет влево, пик растёт и больше больших кристаллов
Coef_nucl_T    = 0.5e+4;     %коэффициент в формуле для скорости нуклеации в экспоненте с разницей температур; при увелчении значения угол наклона графика распределения вначале уменьшается, значит меньше маленьких кристаллов
    elseif (Spec == 2)
T_start   = 801+273;            %начальная температура в кельвинах;
T_end     = 775.3080+273;            %конечная температура в кельвинах;
X_H20     = 6; %массовое содержание воды в расплаве;
M         = 1.4792;
t_end_yrs = 1000;               %колличество лет;
str_Speciment = ('Bishop Tuff LV-13');


L              = 3;            %характерная длина (в метрах) ;   

J_0            = L*t_end_yrs*1e+1;
A              = 2.e+7;
B              = 5.e+2;

Coef_delta_C   = 17.5;
Coef_N         = 8.95;
Coef_Nucl_mult = L*t_end_yrs*9e-3;        %коэффициент в формуле для скорости нуклеации перед экспоентами; при увелчении значения график скорости нуклеации едет влево и график распределения болле пологий
Coef_nucl_c    = 0.06e-2;       %коэффициент в формуле для скорости нуклеации в экспненте с разницей концентраций; при увелчении значения график скорости нуклеации едет влево, пик растёт и больше больших кристаллов
Coef_nucl_T    = 0.2e+4;     %коэффициент в формуле для скорости нуклеации в экспоненте с разницей температур; при увелчении значения угол наклона графика распределения вначале уменьшается, значит меньше маленьких кристалло
    elseif (Spec ==3)
T_start   = 800+273;            %начальная температура в кельвинах;
T_end     = 770.0811+273;            %конечная температура в кельвинах; 
t_end_yrs = 10000;               %колличество лет;
str_Speciment = ('Yellowstone MFT-1');

L              = 2;

J_0            = L*t_end_yrs*1e+1;
A              = 2.e+7;
B              = 5.e+2;

X_H20          = 3; %массовое содержание воды в расплаве;
M              = 1.2873;
Coef_delta_C   = 2.2;
Coef_N         = 2.9;
Coef_Nucl_mult = L*t_end_yrs*9e-3;        %коэффициент в формуле для скорости нуклеации перед экспоентами; при увелчении значения график скорости нуклеации едет влево и график распределения болле пологий
Coef_nucl_c    = 0.06e-2;       %коэффициент в формуле для скорости нуклеации в экспненте с разницей концентраций; при увелчении значения график скорости нуклеации едет влево, пик растёт и больше больших кристаллов
Coef_nucl_T    = 0.4e+4;     %коэффициент в формуле для скорости нуклеации в экспоненте с разницей температур; при увелчении значения угол наклона графика распределения вначале уменьшается, значит меньше маленьких кристаллов
    elseif (Spec ==4)
T_start   = 863+273;            %начальная температура в кельвинах;
T_end     = 855.1239+273;            %конечная температура в кельвинах; 
t_end_yrs = 10000;               %колличество лет;
str_Speciment = ('Yellowstone HRT-A');  

L              = 4;            %характерная длина (в метрах) ;  

J_0            = L*t_end_yrs*1e+3;
A              = 6.2e+6;
B              = 1.e+4;

X_H20           = 3;                 %массовое содержание воды в расплаве;
M              = 1.29;
Coef_delta_C   = 1.32;
Coef_N         = 3.1;
Coef_Nucl_mult = L*t_end_yrs*9e-3;        %коэффициент в формуле для скорости нуклеации перед экспоентами; при увелчении значения график скорости нуклеации едет влево и график распределения болле пологий
Coef_nucl_c    = 0.06e-2;       %коэффициент в формуле для скорости нуклеации в экспненте с разницей концентраций; при увелчении значения график скорости нуклеации едет влево, пик растёт и больше больших кристаллов
Coef_nucl_T    = 0.4e+4;     %коэффициент в формуле для скорости нуклеации в экспоненте с разницей температур; при увелчении значения угол наклона графика распределения вначале уменьшается, значит меньше маленьких кристаллов
    elseif (Spec ==5)
T_start   = 767+273;            %начальная температура в кельвинах;
T_end     = 749.0599+273;            %конечная температура в кельвинах; 
X_H20     = 6; %массовое содержание воды в расплаве;
M         = 1.5046;
t_end_yrs = 10000;               %колличество лет;
str_Speciment = ('Bishop Tuff LV-3');  

L              = 2;            %характерная длина (в метрах) ;  

J_0            = L*t_end_yrs*1e+6;
A              = 7.4e+6;
B              = 1.e+3;

Coef_delta_C   = 0.81;
Coef_N         = 2.8;
Coef_Nucl_mult = L*t_end_yrs*9e-3;        %коэффициент в формуле для скорости нуклеации перед экспоентами; при увелчении значения график скорости нуклеации едет влево и график распределения болле пологий
Coef_nucl_c    = 0.06e-2;       %коэффициент в формуле для скорости нуклеации в экспненте с разницей концентраций; при увелчении значения график скорости нуклеации едет влево, пик растёт и больше больших кристаллов
Coef_nucl_T    = 0.7e+4;     %коэффициент в формуле для скорости нуклеации в экспоненте с разницей температур; при увелчении значения угол наклона графика распределения вначале уменьшается, значит меньше маленьких кристаллов
 elseif (Spec ==6)
T_start   = 809+273;            %начальная температура в кельвинах;
T_end     = 799.2912+273;            %конечная температура в кельвинах; 
X_H20     = 3; %массовое содержание воды в расплаве;
M         = 1.4598;
t_end_yrs = 10000;               %колличество лет;
str_Speciment = ('Timber mountain TM-15');  


L              = 3;            %характерная длина (в метрах) ; 
J_0            = L*t_end_yrs*1e+6;
A              = 7.4e+6;
B              = 1.e+3;

Coef_delta_C   = 7.6;
Coef_N         = 9;
Coef_Nucl_mult = L*t_end_yrs*9e-3;        %коэффициент в формуле для скорости нуклеации перед экспоентами; при увелчении значения график скорости нуклеации едет влево и график распределения болле пологий
Coef_nucl_c    = 0.06e-2;       %коэффициент в формуле для скорости нуклеации в экспненте с разницей концентраций; при увелчении значения график скорости нуклеации едет влево, пик растёт и больше больших кристаллов
Coef_nucl_T    = 0.32e+4;     %коэффициент в формуле для скорости нуклеации в экспоненте с разницей температур; при увелчении значения угол наклона графика распределения вначале уменьшается, значит меньше маленьких кристалл
    end
    
    
    
   %L     = 0.1;
time_coef = (T_start-T_end);    %вспомогательный коэфициент
N_x       = 100001;              %количество узлов в сетке по координате;
N_t       = 1001;               %количество узлов в сетке по времени;

otstup    = 0;                  %отступ от границы кристалла;

%t_end_yrs = 1000;               %колличество лет;
t_end     = t_end_yrs*365*24*60*60; %характерное время в секундах;
C_cryst   = 490000;             %концентрация в кристалле;
C_limit   = 0;                 %предельное значение концентрации, когда кристаллы не появляются;

Drob                = 21;       %константа, отечающая за то, сколько точек будет на графике с распределениями;
Drob_figuer_counter = 0.1;      %константа, отвечающая за то, как часто обновлять графики;

raspr      = zeros(cef_end+1,Drob);
raspr_time = [0.02,0.2,0.5,1];

eps                 = 5;
f                   = 1; %вспомогательная величина для обновления графиков;
Numb_cryst_begin    = 0; %колличество кристаллов в начальный момент времени, не считая граничных;
%C_sat = C_cryst/(exp(10108/T+1.16*(M-1)-1.48));
C_sat               = 490000/(exp(12900/T_start-0.85*(M-1)-3.80));
C_bound             = C_sat;
g                   = 0; %мощность источника

t_end_nd            = 1; %окончание по безразмерному времени;
L_nd                = 1; %окончание по безразмерной координате;
S0                  = 5*1.e-6/L; %начальный радиус кристалла в безразмерных координатах;
C_sat_start         = C_sat;
h_nd   = L_nd/((N_x-1));     % расчетный шаг сетки по безразмерной пространственной координате;
tau_nd = t_end_nd/((N_t)-1); % расчетный шаг сетки по безразмерному времени;


N_changed_x         = ceil(S0/h_nd)+9;
%N_changed_x        = 10;         %колличетво первоначальных ячеек, которые отступаются от средины кристалла для измельчения (в одну сторону); 
N_lost_steps        = 0;
N_razb_x            = 20;        %колличетво ячеек, на которые раздрабливаются исходные; 
N_changed_t         = 10;         %колличество первоначальных ячеек по времени до того, как перед ними  рождается кристалл;
N_razb_t            = 20;         %колличество первоначальных ячеек по времени после того, как перед ними  рождается кристалл;
N_t_new             = N_t+N_razb_t-N_changed_t; %колличество получившехся узлов по времени;
time                = 0;          %безразмерное время;
V_Nucleation        = 0;          %скорость нуклеации;



FLAG_right_progonka         = 1;
FLAG_SIDECRYSTALLS_INCLUDED = 0;
MAX_CRYST_ONETIME           = 100; %максимальное количество кристаллов за шаг по времени;
FLAG_v_limits               = 1;   %ограничитель скорости роста кристаллов;
time_bound_coef             = 4;   %коэфициент для ограничения скорости роста кристалла в момент рождения;
FLAG_v_nucleation_limit     = 1;   %ограничитель скорости нуклеации кристаллов;
FLAG_Time_Limit             = 1;   %флаг, отвечающий за дробление по фремени во время появления кристаллов;
FLAG_drob_x                 = 1;   %флаг, отвечающий за дробление по координате во время появления кристаллов;
FLAG_1                      = 0;   %вспомогательный флаг;
FLAG_2                      = 0;   %вспомогательный флаг;
Volume_Percent              = 0;   %объёмная доля кристаллов во всём расплаве;
stability_coef_bound        = 0;   %коэффициент отвечает за то, как быстро падает концентрация источников вне границы кристаллов;
x_C                         = linspace(0,200,1000); %вспомогательный массив для потсроения графика;

Graph = figure('Units','normalized','Position',[.10,.10,.8,.8]);
ha = axes('Units','normalized','Position',[.05,.48,.40,0.47]);
      
set(Graph,'Position',[.10,.10,.8,.8]);
box on
axes_raspr = axes('Units','normalized','Position',[0.05,0.15,0.40,0.23]); box on
axes_vel   = axes('Units','normalized','Position',[0.5,0.65,0.40,0.30]); box on
axes_tempr = axes('Units','normalized','Position',[0.5,0.4,0.40,0.15]); box on
axes_V     = axes('Units','normalized','Position',[0.5,0.15,0.40,0.15]); box on
      
xlabel(axes_vel,'Time (years)')
ylabel(axes_vel,'Zr radius')

xlabel(axes_tempr,'Time (years)')
ylabel(axes_tempr,'Growth Rate, cm.s^{-1}')

xlabel(axes_V,'Time (years)')
ylabel(axes_V,'Temperature T, ^oC')

set(Graph,'Name','DIFFUSOR')
movegui(Graph,'center')
set(Graph,'NumberTitle','off')




T_Array              = zeros(1,1);
V_Nucleation_Array   = zeros(1,1);
Volume_Percent_Array = zeros(1,N_t);
ti                   = zeros(1,N_t); %массив координат по времени;
ti_nd                = 0*(0:0)*t_end_yrs; %массив координат по времени обезразмеренной;
x_nd                 = h_nd*(0:N_x-1); %массив безразмерных пространственных координат;
S_t                  = zeros(1,N_t); %массив для отслеживания численного значения S(t);
S_t_teor             = zeros(1,N_t); %массив для отслеживания теоретического значения S(t);
C                    = zeros(1,N_x); %массив концентрации;

alfa = zeros(1,N_x); %массив для прогоночных коффициентов;
beta = zeros(1,N_x); %массив для прогоночных коффициентов;

t_buffer    = zeros(1,N_razb_t);
p           = ones(1,N_razb_t)*tau_nd;
p(N_razb_t) = -(N_changed_t-1)*tau_nd;
t_roots     = roots(p);
b           = t_roots(~logical(imag(t_roots)));
b           = b(b>0);

%высчитывание длин тех отрезков, которыми заменяем
for i=1:N_razb_t
    t_buffer(i) = tau_nd*b^(N_razb_t-i);
end

S_bound_count_right = zeros(1,N_cryst);
S_bound_count_left  = zeros(1,N_cryst);
x_buffer            = zeros(1,N_razb_x);
S_cent              = zeros(1,N_cryst); %массив для левых границ кристаллов;
S_right             = zeros(1,N_cryst); %массив для левых границ кристаллов;
S_left              = zeros(1,N_cryst); %массив для правых границ кристаллов;
dCpodX_right        = zeros(1,N_cryst); %массив градиентов для левых границ кристаллов;
dCpodX_left         = zeros(1,N_cryst); %массив градиентов для правых границ кристаллов;
g_s                 = zeros(1,N_cryst); %массив для плотностей источников;
C_after_right       = zeros(1,N_cryst);
C_befor_right       = zeros(1,N_cryst);
C_after_left        = zeros(1,N_cryst);
C_befor_left        = zeros(1,N_cryst);
V_last_right        = zeros(1,N_cryst);
V_last_left         = zeros(1,N_cryst);
Crystal_next        = (1:N_cryst);
max_num_array       = zeros(1,N_cryst);
max_num_array(1)    = 1;
max_num_array(2)    = N_x;
S_check             = zeros(1,4);
S_last_right        = zeros(1,N_cryst);
S_last_left         = zeros(1,N_cryst);
N_steps             = zeros(1,N_cryst);
Coef_1              = zeros(1,N_cryst);
Coef_2              = zeros(1,N_cryst);
Crystal_growth      = zeros(1,N_cryst); %массив для размеров кристаллов;
Crystal_size        = zeros(1,N_cryst);
Crystal_size_x      = linspace(0,1,Drob);
Crystal_size_y      = linspace(0,1,Drob);
x_teor_drob         = linspace(0,1,Drob);
y_teor_drob         = linspace(0,1,Drob);
raspr_prom          = zeros(4,Drob);

%{
y_teor_drob(1)  = -5;
y_teor_drob(2)  = -1.25;
y_teor_drob(3)  = -0.5;
y_teor_drob(4)  = -0.25;
y_teor_drob(5)  = 0;
y_teor_drob(6)  = -0.25;
y_teor_drob(7)  = -0.5;
y_teor_drob(8)  = -0.75;
y_teor_drob(9)  = -1;
y_teor_drob(10) = -1.5;
y_teor_drob(11) = -2.5;
%}


y_teor_drob(1)  = -5;
y_teor_drob(2)  = -3.25;
y_teor_drob(3)  = -1.5;
y_teor_drob(4)  = -0.5;
y_teor_drob(5)  = 0;
y_teor_drob(6)  = -0.5;
y_teor_drob(7)  = -1.8;
y_teor_drob(8)  = -1.75;
y_teor_drob(9)  = -3.25;
y_teor_drob(10) = -3.25;
y_teor_drob(11) = -3.75;

HRT_C_raspr =[0.11088,	0.158996,	0.239538,	0.285264,	0.369737,	0.408577,	0.451856,	0.49503,	0.54475,	0.581224,	0.617807,	0.704456,	0.754303,	0.788375,	0.833986,	0.876848,	0.920427,	1.00239;...
             -3.25333,	-1.992,	-0.833333,	-0.144,	0.00266667,	-0.305333,	-0.334667,	-0.628,	-0.789333,	-1.61067,	-2.15333,	-1.97733,	-1.816,	-3.23867,	-2.84267,	-3.928,	-3.19467,	-3.928];
MFT_1_raspr       = [0.0525426,	0.0943261,	0.137887,	0.181224,	0.222326,	0.270061,	0.315608,	0.358926,	0.404387,	0.438984,	0.48228,	0.529824,	0.568791,	0.614136,	0.659544,	0.696294,	0.741755,	0.7871,	0.826084,	0.871377,	0.91032,	0.957835,	1.00084;...
              -2.98933,	-1.31733,	-0.628,	-0.510667,	-0.569333,	-0.276,	-0.0413333,	0.032,	0.0466667,	-0.0413333,	-0.0266667,	-0.217333,	-0.202667,	-0.481333,	-0.598667,	-0.716,	-0.701333,	-0.98,	-0.921333,	-1.332,	-1.376,	-1.64,	-2.35867];

HRT_A_raspr = [21.7068645640070,32.8385899814470,43.4137291280150,54.5454545454540,74.5825602968460,89.6103896103900,95.1762523191090,106.307977736550,120.222634508350,128.571428571430,136.363636363640,150.834879406310,161.966604823750,170.871985157700,184.230055658630,204.823747680890,217.068645640070;8.91585552790700,9.76402383443680,9.80327503685660,9.43817320120640,8.80879400446090,8.52001395129360,8.30586745684920,8.24408315674390,8.01835357150180,7.81719702098850,8.00775527789750,7.82001072725510,7.55601473678660,7.19063153050970,6.80053460419070,6.56301090018080,6.19804974984390];

HRT_C_raspr_try = [24.4821,	34.0386,	42.7092,	55.0554,	64.5076,	76.7087,	83.9091,	93.1785,	106.097,	114.661,	123.962,	136.018,	145.262,	155.333,	164.663,	175.449,	183.821,	198.242,	205.225,	215.991,	234.651;
                   6.49557,	7.7572,	8.08137,	8.89078,	9.57059,	9.5719,	9.70197,	9.36358,	9.36497,	9.09114,	8.93053,	8.12375,	7.63989,	7.77027,	7.77128,	7.88557,	6.54506,	6.91833,	5.83625,	5.83741,	5.83942];
LV_13 = [6.50759,	14.9675,	26.6811,	45.5531,	55.9653,	65.7267,	78.0911,	86.551,	    96.3124,	101.518,	107.375,	115.835,	125.597,	136.659,	156.182,	165.944,	185.466,	195.228;
                   8.47928,	10.0418,	9.93109,	8.37062,	8.03671,	7.36807,	7.48055,	6.76718,	6.43322,	5.5188,	    5.13993,	4.8505,	    4.20416,	3.84799, 	 4.1395,    3.24773,	3.24918,	2.73673];
 
LV_3  = [     15,	  23.4783,	34.5652,	43.0435,	54.1304,	62.6087,	75,	83.4783,	94.5652,	101.739,	114.783,	123.913,	133.696,	143.478,	155.217,	165,	174.783,	184.565,	194.348,	205.435;...
         3.95935,	6.49593,	6.23577,	7.34146,	7.66667,	8.04065,	7.87805,	8.12195,	8.07317,	7.79675,	7.53659,	7.34146,	7.09756,	6.85366,	6.73984,	6.56098,	5.92683,	5.42276,	4.70732,	4.69106];

TM_15 = [13.0435,	23.4783,	35.2174,	44.3478,	56.087,	63.913,	74.3478,	84.7826,	93.913,	104.348,	116.087,	122.609,	134.348,	143.478,	155.217,	164.348,	174.783,	185.217,	203.478,	284.348,	293.478;...
        6.41463,	7.13008,	7.35772,	8.82114,	9.47154,	9.47154,	8.9187,	8.85366,	7.84553,	7.81301,	7.71545,	7.09756,	6.6748,	6.5122,	6.70732,	6.41463,	5.27642,	5.69919,	4.23577,	3.19512,	3.09756];
    
YL_4 = [15.6398,	24.8815,	54.7393,	76.0664,	83.8863,	90.9953,	100.948,	113.033,	125.118,	135.782,	145.735,	157.82,	171.327,	177.014,	196.209,	228.199,	274.408;...
        4.35603,	6.11112,	7.1648,	7.58525,	7.5844,	7.43729,	6.98093,	6.65441,	6.7344,	6.48935,	5.59396,	5.52761,	5.15216,	4.68,	4.22263,	4.0403,	3.07594];

    
LCT_3a = [23.4597,	34.1232,	42.654,	65.4028,	73.2227,	81.7536,	92.4171,	105.213,	113.744,	122.986,	132.938,	144.313,	154.265,	163.507,	175.592,	184.834,	193.365,	204.028,	214.692,	245.261;...
          6.3877,	5.88248,	5.93034,	8.17177,	8.36605,	8.98301,	8.36397,	8.4764,	8.03645,	7.58016,	7.00998,    6.84614,	5.88572,	5.51073,	4.77771,	5.42712,	5.45871,	4.80715,	4.74095,	3.59943];
    
HRT_A_raspr(1,:) = HRT_A_raspr(1,:)/max(HRT_A_raspr(1,:));
HRT_A_raspr(2,:) = log(exp(HRT_A_raspr(2,:))/exp(max(HRT_A_raspr(2,:))));       
               
LV_13(1,:) = LV_13(1,:)/max(LV_13(1,:));
LV_13(2,:) = log(exp(LV_13(2,:))/exp(max(LV_13(2,:))));              
               
HRT_C_raspr_try(1,:) = HRT_C_raspr_try(1,:)/max(HRT_C_raspr_try(1,:));
HRT_C_raspr_try(2,:) = log(exp(HRT_C_raspr_try(2,:))/exp(max(HRT_C_raspr_try(2,:))));

LV_3(1,:) = LV_3(1,:)/max(LV_3(1,:));
LV_3(2,:) = log(exp(LV_3(2,:))/exp(max(LV_3(2,:)))); 

TM_15(1,:) = TM_15(1,:)/max(TM_15(1,:));
TM_15(2,:) = log(exp(TM_15(2,:))/exp(max(TM_15(2,:)))); 

YL_4(1,:) = YL_4(1,:)/max(YL_4(1,:));
YL_4(2,:) = log(exp(YL_4(2,:))/exp(max(YL_4(2,:)))); 

LCT_3a(1,:) = LCT_3a(1,:)/max(LCT_3a(1,:));
LCT_3a(2,:) = log(exp(LCT_3a(2,:))/exp(max(LCT_3a(2,:)))); 

p = polyfit(LV_13(1,:), LV_13(2,:), 10);
y1 = polyval(p,Crystal_size_x);

a      = fzero(@(p_nd) pi^(1/2)*p_nd*exp(p_nd^2)*erfc(p_nd)-(C_bound-C_sat)/(C_cryst-C_sat),0); %вычисление промежуточной величины для аналитического решения;
D      = exp(-(11.4*X_H20+3.13)/(0.84*X_H20+1)-((21.4*X_H20+47)/(1.06*X_H20+1))*(1000)/T_start);
D_nd   = D*t_end/L^2;

time_bound = time_bound_coef*tau_nd;
V_bound    = a*(D_nd/time_bound)^(1/2);
V_right    = zeros(1,N_cryst); %массив скоростей для левых границ кристаллов;
V_left     = zeros(1,N_cryst); %массив скоростей для правых границ кристаллов;

p_x           = ones(1,N_razb_x)*h_nd;
p_x(N_razb_x) = -(N_changed_x-1)*h_nd;
x_roots       = roots(p_x);
b_x           = x_roots(~logical(imag(x_roots)));
b_x           = b_x(b>0);

  
for i=1:N_razb_x
    x_buffer(i) = h_nd*b_x^(N_razb_x-i);
end
    
for i=2:N_razb_x
    x_buffer(i) = x_buffer(i-1)+x_buffer(i);
end 



%начальные условия
for j=1:(N_x)
 C(j) = C_bound;
end

S_cent(1)  = 0;
S_right(1) = S0;
S_left(1)  = 0;

S_cent(2)  = 1;
S_right(2) = 1;
S_left(2)  = 1-S0;

%начальные условия
for j=1:2
k = 1;
        while x_nd(k)<S_right(j)
        k = k+1;
        end
      k=k-1;  
        
l = 1;
        while x_nd(l)<S_left(j)
        l=l+1;
        end
        
for s=l:k
    C(s) = C_sat;
end
end

        n=1;

        if(FLAG_drob_x==1) 
        k = 0;
        while S_cent(n)+h_nd*k<=S_right(n)
        k = k+1;
        end  
        k = k+1;
  
        for k = (N_x+(N_razb_x-N_changed_x)):-1:max_num_array(n)+(N_razb_x-N_changed_x)
            x_nd(k) = x_nd(k-(N_razb_x-N_changed_x));
        end 

        for k=max_num_array(n)+1:max_num_array(n)+N_razb_x
            x_nd(k) = x_nd(max_num_array(n))+x_buffer(k-max_num_array(n));
        end
        
        % for k=max_num_array(n)-1:-1:max_num_array(n)-N_razb_x
        %    x_nd(k) = x_nd(max_num_array(n))-x_buffer(max_num_array(n)-k);
        % end
        
        c_buffer = zeros(1,N_changed_x+1);
        
        for k=1:N_changed_x+1
          c_buffer(k) = C(max_num_array(n)-(-N_changed_x)-N_changed_x+k-1);
        end
        
        for k = (N_x+(N_razb_x-N_changed_x)):-1:max_num_array(n)+N_razb_x-N_changed_x
            C(k) = C(k-(N_razb_x-N_changed_x));
        end 
        
        k=max_num_array(n);
        for p=1:N_changed_x     
            while x_nd(k)<=x_nd(max_num_array(n))+h_nd*(p-N_changed_x) && x_nd(k)>=x_nd(max_num_array(n))+h_nd*(p-1-N_changed_x)
                C(k) = c_buffer(p)+(c_buffer(p+1)-c_buffer(p))*(x_nd(k)-(x_nd(max_num_array(n))+h_nd*(p-1-N_changed_x)))/(h_nd);
                k    = k+1;
            end 
        end
      
        N_x = N_x+(N_razb_x-N_changed_x);

        max_num_array(2)       = max_num_array(2)+(N_razb_x-N_changed_x);
        S_bound_count_right(2) = N_x;
        
        end
        
      
        
        n=2;

        if(FLAG_drob_x==1) 
        k = 0;
        while S_cent(n)+h_nd*k<=S_left(n)
        k = k+1;
        end  
        k = k+1;
  
        for k = (N_x+(N_razb_x-N_changed_x)):-1:max_num_array(n)+(N_razb_x-N_changed_x)
            x_nd(k) = x_nd(k-(N_razb_x-N_changed_x));
        end 

        %for k=max_num_array(n)+1:max_num_array(n)+N_razb_x
        %    x_nd(k) = x_nd(max_num_array(n))+x_buffer(k-max_num_array(n));
        %end
        
        max_num_array(2)=max_num_array(2)+N_razb_x-N_changed_x;
        for k=max_num_array(n)-1:-1:max_num_array(n)-N_razb_x
            x_nd(k) = x_nd(max_num_array(n))-x_buffer(max_num_array(n)-k);
        end
        
        c_buffer = zeros(1,N_changed_x+1);
        
        for k=N_changed_x+1:-1:1
          c_buffer(k) = C(max_num_array(n)-N_razb_x+N_changed_x-k-1);
        end
        
        for k = (N_x+(N_razb_x-N_changed_x)):-1:max_num_array(n)-N_razb_x+N_changed_x
            C(k) = C(k-(N_razb_x-N_changed_x));
        end 
        
        k=max_num_array(n);
        for p=-N_changed_x:-1     
            while x_nd(k)<=x_nd(max_num_array(n))+h_nd*(p-N_changed_x) && x_nd(k)>=x_nd(max_num_array(n))+h_nd*(p-1-N_changed_x)
                C(k) = c_buffer(p)+(c_buffer(p+1)-c_buffer(p))*(x_nd(k)-(x_nd(max_num_array(n))+h_nd*(p-1-N_changed_x)))/(h_nd);
                k    = k+1;
            end 
        end
      
        N_x = N_x+(N_razb_x-N_changed_x);

        end
        %начальные условия
      for j=1:2
k = 1;
        while x_nd(k)<S_right(j)
        k = k+1;
        end
      k=k-1;  
        
l = 1;
        while x_nd(l)<S_left(j)
        l=l+1;
        end
        
for s=l:k
    C(s) = C_sat;
end
end  
        
        
     

 %граничное условие
        n     = 2;
        index = 1;
        sum_2 = 0;
        n_img        = n;
        Time_Counter = 1;
        S_bound_count_right(2) = N_x;
        
       
    while time < t_end_nd
        if(V_Nucleation>exp(13))
            1100110011
        break
        end
        
        Time_Counter = Time_Counter+1;
        %часть, отвечающая за изменеие температуры
        T     = T_start-(T_start-T_end)*time/(0.05*cef); %температура в кельвинах;
        %M     = 4.8*10^(-6)*T^2-8.4*10^(-3)*T+4.84;
        %C_sat = C_cryst/(exp(10108/T+1.16*(M-1)-1.48));
        C_sat = 490000/(exp(12900/T-0.85*(M-1)-3.80));
        a     = fzero(@(x_nnd) pi^(1/2)*x_nnd*exp(x_nnd^2)*erfc(x_nnd)-(C_bound-C_sat)/(C_cryst-C_sat) ,0); %вычисление промежуточной величины для аналитического решения;
        D     = (exp(-(11.4*X_H20+3.13)/(0.84*X_H20+1)-((21.4*X_H20+47)/(1.06*X_H20+1))*(1000)/T));
        
        D_nd  = D*t_end/L^2;
        
        if(T >= T_end)
        T_Array(Time_Counter) = T-273;
        else
        T_Array(Time_Counter) = T_end-273;
        end
        
        if (T<T_end)
            time = 1;
        end
        
        if (C_sat<C_limit)
            FLAG_v_nucleation_limit=1;
        end
        
        if(FLAG_Time_Limit==0)
        index = N_razb_t;
        end
        
        if (index<N_razb_t)
        index = index+1;
        end
        
        tau_nd=t_buffer(index);
        time = time + tau_nd;
        
        S_teor      =(2*a*(D_nd*time)^(1/2));
        V_teor      = a*(D_nd/time)^(1/2);
        dCpodX_teor = -V_teor/D_nd*(C_sat-C_cryst);
        
        max_delta=0;
        for i=1:N_x
        if(max_delta<abs(C(i)-C_sat) && C(i)>=C_sat)
        max_delta   = abs(C(i)-C_sat);
        end
        end
        
        if(n>N_cryst-1)
          FLAG_v_nucleation_limit = 1;
        end
        
        if (FLAG_v_nucleation_limit == 1)
            V_Nucleation = 0;
        else
            V_Nucleation = Coef_Nucl_mult*exp((max_delta/Coef_delta_C)^Coef_N*Coef_nucl_c)*exp(-1/T*(Coef_nucl_T));
            %V_Nucleation = V_nucleation(J_0,A,B,T,max_delta);
        end
        
        
        V_Nucleation_Array(Time_Counter) = V_Nucleation;
        Numb_cryst_begin                 = Numb_cryst_begin+V_Nucleation*tau_nd;
        if(V_Nucleation*tau_nd>1)
        k=0;  
        end
        ti_nd(Time_Counter)              = t_end_yrs*time;
        
        %%часть, создающая кристаллы;
        for index_2=1:MAX_CRYST_ONETIME
        if (Numb_cryst_begin>(n_img-1))
        index      = 1;
        max_num    = 1;
        max_dist   = 0;
        min_dist_i = 0;
        max_value  = C(1)-C_sat;
      
        %N_changed_t         = 10;         %колличество первоначальных ячеек по времени до того, как перед ними  рождается кристалл;
        %N_razb_t
        %ti_nd_buffer                 = zeros(1,size(ti_nd,2)+N_razb_t-N_changed_t);
        %ti_nd_buffer(1:Time_Counter) = ti_nd(1:Time_Counter);
        %ti_nd_buffer(Time_Counter:Time_Counter+N_razb_t) = ti_nd(Time_Counter:Time_Counter+N_razb_t);
        %нахождение точки с наибольшей разницей концентраци;
        %{
        for i=1:N_x
            if C(i)>C_bound
                C(i)=C_bound;
            end
        end
        %}
      
        
        max_value_test = max(C-C_sat)-0;
        if(n~=2)
        [max_value,max_num] = max(C);
        
        k        = 1;
        x_nd_max = x_nd;
        %C(i)>=max_value-eps && C(i)<=max_value+eps
        for i=1:N_x
            if (C(i)<=max_value-eps)
                x_nd_max(i)=0;
            end
        end
        
        S_count=zeros(n,N_x);
        for k=1:n
        S_count(k,:)=abs(S_cent(k)-x_nd_max);
        end
        zxcvb           = min(S_count);
        [S_count_max_min,max_num] = max(min(S_count));

        max_num_array(n+1) = max_num;
        else
        max_num            = ceil(N_x/2);    
        max_num_array(n+1) = max_num;
        end
       
        if(max_value+C_sat>C_limit)
        n          = n+1;
        n_img      = n_img+1;
        S_cent(n)  = x_nd(max_num);
        S_right(n) = S_cent(n)+S0;
        S_left(n)  = S_cent(n)-S0;
        %вычисление координаты, наиболее близкой к правой границе слева 
        k          = max_num;
        while x_nd(k)<S_right(n)
        k = k+1;
        end
        k = k-1;
        
        %вычисление координаты, наиболее близкой к левой границе справа
        l = max_num;
        while x_nd(l)>S_left(n)
        l = l-1;
        end
        
        %{
        for r=l:k
        C(r)= C_sat;
        end
        %}
        
        for r=1:2*n
        for p=1:n-1       
        if (S_cent(Crystal_next(p))>S_cent(Crystal_next(p+1)))
            m                 = Crystal_next(p);
            Crystal_next(p)   = Crystal_next(p+1);
            Crystal_next(p+1) = m;
        end
        end
        end
        else
            n_img=n_img+1;
        end
        
        if(FLAG_drob_x==1) 
        k = 0;
        while S_cent(n)+h_nd*k<=S_right(n)
        k = k+1;
        end  
        k = k+1;
  
        for k = (N_x+2*(N_razb_x-N_changed_x)):-1:max_num_array(n)+2*(N_razb_x-N_changed_x)
            x_nd(k) = x_nd(k-2*(N_razb_x-N_changed_x));
        end 
        
        max_num_array(n)       = max_num_array(n)+N_razb_x-N_changed_x;
        x_nd(max_num_array(n)) = x_nd(max_num_array(n)-(N_razb_x-N_changed_x));
 
        for k=max_num_array(n)+1:max_num_array(n)+N_razb_x
            x_nd(k) = x_nd(max_num_array(n))+x_buffer(k-max_num_array(n));
        end
        
         for k=max_num_array(n)-1:-1:max_num_array(n)-N_razb_x
            x_nd(k) = x_nd(max_num_array(n))-x_buffer(max_num_array(n)-k);
         end
        
        c_buffer = zeros(1,2*N_changed_x+1);
        
        for k=1:2*N_changed_x+1
          c_buffer(k) = C(max_num_array(n)-(N_razb_x-N_changed_x)-N_changed_x+k-1);
        end
        
        for k = (N_x+2*(N_razb_x-N_changed_x)):-1:max_num_array(n)+N_razb_x-N_changed_x
            C(k) = C(k-2*(N_razb_x-N_changed_x));
        end 
        k=max_num_array(n)-N_changed_x-(N_razb_x-N_changed_x);
         
        %{
        for p=1:2*N_changed_x+1     
            while x_nd(k)<=x_nd(max_num_array(n))+h_nd*(p-N_changed_x) && x_nd(k)>=x_nd(max_num_array(n))+h_nd*(p-1-N_changed_x)
                C(k) = c_buffer(p)+(c_buffer(p+1)-c_buffer(p))*(x_nd(k)-(x_nd(max_num_array(n))+h_nd*(p-1-N_changed_x)))/(h_nd);
                k    = k+1;
            end
        end
        %}
            
        N_x = N_x+2*(N_razb_x-N_changed_x);
        
        for i=1:N_x
           if (x_nd(i) <= S_cent(n)+S0&& x_nd(i)>= S_cent(n))
               C(i)=(c_buffer(1)-C_sat)*(x_nd(i)-S_cent(n))/(S0)+C_sat;
           end
           if (x_nd(i) >= S_cent(n)-S0&& x_nd(i)<= S_cent(n))
               C(i)=(c_buffer(1)-C_sat)*(S_cent(n)-x_nd(i))/(S0)+C_sat;
           end
        end
    
       %{
        for k=1:N_x
        if(x_nd(k)<=S_cent(n)+S0 && x_nd(k)>=S_cent(n)-S0)
               C(k)=C_sat;
        end
        end
        %}
        
        r   = n;
        for r=1:n
            if(Crystal_next(r)==n)
                numb=r;
            end
        end
        numb=numb+1;
        while Crystal_next(numb)~=2
            max_num_array(Crystal_next(numb))=max_num_array(Crystal_next(numb))+2*(N_razb_x-N_changed_x);
            numb=numb+1;
        end
        max_num_array(2)       = max_num_array(2)+2*(N_razb_x-N_changed_x);
        S_bound_count_right(2) = N_x;
        end
        end
   
        
        end
 
        %часть, задающая источник во внутренних кристаллах
        
        for j=3:n
   
        V_last_right(j) = V_right(j);
        S_last_right(j) = S_right(j);    
        S_right(j)      = S_last_right(j)+V_last_right(j)*tau_nd;
         
        V_last_left(j)  = V_left(j);
        S_last_left(j)  = S_left(j);    
        S_left(j)       = S_last_left(j)+V_last_left(j)*tau_nd;    
            
        k=max_num_array(j);
        while x_nd(k)<S_right(j)
        k=k+1;
        end
        S_bound_count_right(j) = k;
        k                      = k+otstup;
        
        C_after_right(j) = C(k);
        C_befor_right(j) = C(k-1);
        
        m=max_num_array(j);
        while x_nd(m)>S_left(j)
        m = m-1;
        end
        m = m+1;
        S_bound_count_left(j) = m-1;
        m = m-otstup;

        C_after_left(j) = C(m);
        C_befor_left(j) = C(m-1);
        dCpodX_left(j)  = (C_after_left(j)-C_befor_left(j))/(x_nd(m)-x_nd(m-1));
        dCpodX_right(j) = (C_after_right(j) - C_befor_right(j))/(x_nd(k) - x_nd(k-1));
        
        if (N_steps(j)>N_lost_steps)
        V_left(j)  = -D_nd*dCpodX_left(j)/(C_sat-C_cryst);
        V_right(j) = -D_nd*dCpodX_right(j)/(C_sat-C_cryst);
        else
        V_left(j)  = 0;
        V_right(j) = 0;
        end
        if j==2
            V_right(j) = 0;
        end
        if(FLAG_v_limits==1)
            N_steps(j)=N_steps(j)+1;
        if (V_left(j)<-a*(D_nd/time_bound)^(1/2)||V_left(j)>0)
            V_left(j)      = -a*(D_nd/time_bound)^(1/2);
            dCpodX_left(j) = -V_left(j)*(C_sat-C_cryst)/D_nd;
        end
        if (V_right(j)>a*(D_nd/time_bound)^(1/2)||V_right(j)<0)
            V_right(j)      = a*(D_nd/time_bound)^(1/2);
            dCpodX_right(j) = -V_right(j)*(C_sat-C_cryst)/D_nd;
        end
        end
       
        Coef_1(j)=-6 * D_nd * (dCpodX_right(j) + dCpodX_left(j)) / (S_right(j) ^ 2 - 2 * S_right(j) * S_left(j) + S_left(j) ^ 2);
        Coef_2(j)=2*(S_left(j)*dCpodX_left(j) + 2*S_right(j)*dCpodX_left(j) + 2*S_left(j)*dCpodX_right(j) + S_right(j)*dCpodX_right(j))*D_nd/(S_left(j) - S_right(j))^2;
        end
        
        %часть, задающая плотность источников в граничных кристаллах;  
        
        V_last_right(1) = V_right(1);
        S_last_right(1) = S_right(1);
        S_right(1)      = S_last_right(1)+V_last_right(1)*tau_nd; 
        
        k = 1;
        while x_nd(k)<=S_right(1)
        k = k+1;
        end
        S_bound_count_right(1) = k;
        k = k+otstup;
        
        C_after_right(1) = C(k);
        C_befor_right(1) = C(k-1);
               
        dCpodX_right(1) = ((C_after_right(1) - C_befor_right(1))/(x_nd(k) - x_nd(k-1)));
        V_right(1)      = -D_nd*dCpodX_right(1)/(C_sat-C_cryst);
        if(FLAG_v_limits==1)
        if (V_right(1)>a*(D_nd/time_bound)^(1/2)||V_right(1)<0)
            V_right(1)      = a*(D_nd/time_bound)^(1/2);
            dCpodX_right(1) = -V_right(1)*(C_sat-C_cryst)/D_nd;
            FLAG_1 = 1;
        else
            FLAG_1 = 0;
        end
        end
        g_s(1)    = (-C_after_right(1) + C_sat)*D_nd/((-x_nd(k) + S_right(1))*S_right(1));
        g_s(1)    = -dCpodX_right(1)*D_nd/(S_right(1));
        g         = g_s(1);
        Coef_2(1) = g_s(1);
     
        V_last_left(2) = V_left(2);
        S_last_left(2) = S_left(2);
        S_left(2)      = S_last_left(2)+V_last_left(2)*tau_nd; 
        
        k=1;
        while x_nd(k)<=S_left(2)
        k = k+1;
        end
        S_bound_count_left(2) = k-1;
        k               = k-otstup;
        C_after_left(2) = C(k);
        C_befor_left(2) = C(k-1);
        
        dCpodX_left(2) = (C_after_left(2)-C_befor_left(2))/(x_nd(k)-x_nd(k-1));
        V_left(2)      = -D_nd*dCpodX_left(2)/(C_sat-C_cryst);
        if(FLAG_v_limits==1)
        if (V_left(2)<-a*(D_nd/time_bound)^(1/2)||V_left(2)>0)
            V_left(2)      = -a*(D_nd/time_bound)^(1/2);
            dCpodX_left(2) = -V_left(2)*(C_sat-C_cryst)/D_nd;
            FLAG_2         = 1;
        else
            FLAG_2         = 0;
        end
        end
        g_s(2)    = -dCpodX_left(2)*D_nd/(S_left(2) - 1);
        Coef_2(2) = g_s(2);
       
      
        %граничное условие 
        Crystal_growth(1) = (S_right(1)*2-2*S0)*L;
        Crystal_growth(2) = (2*(1-S_left(2))-2*S0)*L;
        
    for i=3:N_cryst
        Crystal_growth(i) = ((S_right(i)-S_left(i))-2*S0)*L;
    end
    
        Crystal_size(1) = (S_right(1)*2-2*S0);
        Crystal_size(2) = (2*(1-S_left(2))-2*S0);

    for i=3:n
       Crystal_size(i)  = abs(((S_right(i)-S_left(i)))-2*S0);
    end
    
       Crystal_size_y   = zeros(1,Drob);
       Crystal_size_max = Crystal_size(1);
       first_cryst      = 1;
    
    if(FLAG_SIDECRYSTALLS_INCLUDED==0)
    if(n>=3)
       Crystal_size_max = Crystal_size(3);
    end
    first_cryst      = 3;
    end
   
    for m=first_cryst:n
        if(Crystal_size_max<Crystal_size(m))
            Crystal_size_max = Crystal_size(m);
        end
    end
    
    qsc=Crystal_size/Crystal_size_max;
    
    %Crystal_size_max=9.208604625814420e-06;
    
    for i=first_cryst:n
        r = 1;
     while Crystal_size_x(r)<(Crystal_size(i)/Crystal_size_max)
        r = r+1;
     end
        Crystal_size_y(r)  = Crystal_size_y(r)+1;
    end
    
        Number_crystals_max = 1;
    for j=1:Drob
        if (Crystal_size_y(j)>=Number_crystals_max)
            Number_crystals_max = Crystal_size_y(j);
        end
    end
    %Number_crystals_max = 19;
    
         %Crystal_size_y      = log(abs(Crystal_size_y/Number_crystals_max));
         Crystal_size_y_test = Crystal_size_y;
   
     
       
            Crystal_size_sum = 0;
            Proc_Cryst = 1.13656*(1.e-9)*(T-273)^4 - 3.8926*(1.e-6)*(T-273)^3+4.9848*(1.e-3)*(T-273)^2-2.832*(T-273)+602.99;
            
        for i=1:n
            Crystal_size_sum = Crystal_size_sum+Crystal_size(i)+S0*2*L;
        end
            Volume_Percent                     = Crystal_size_sum*100/L;
            Volume_Percent_Array(Time_Counter) = Volume_Percent;
            j                                  = 1;
            FLAG_coef_prog                     = 0;
           
     for q=1:n-1
         
            alfa(max_num_array(Crystal_next(q))) = 0;
            beta(max_num_array(Crystal_next(q))) = -Coef_1(Crystal_next(q))*S_cent(Crystal_next(q))^3/(6*D_nd) - Coef_2(Crystal_next(q))*S_cent(Crystal_next(q))^2/(2*D_nd) + ((Coef_1(Crystal_next(q))*S_left(Crystal_next(q))^2 + 2*dCpodX_left(Crystal_next(q))*D_nd + 2*Coef_2(Crystal_next(q))*S_left(Crystal_next(q)))*S_cent(Crystal_next(q)))/(2*D_nd) + (-2*Coef_1(Crystal_next(q))*S_left(Crystal_next(q))^3 - 6*D_nd*S_left(Crystal_next(q))*dCpodX_left(Crystal_next(q)) - 3*Coef_2(Crystal_next(q))*S_left(Crystal_next(q))^2 + 6*C_sat*D_nd)/(6*D_nd);
            %alfa(max_num_array(Crystal_next(q+1))) = 0;
            %beta(max_num_array(Crystal_next(q+1))) = -Coef_1(Crystal_next(q + 1))*S_cent(Crystal_next(q + 1))^3/(6*D_nd) - Coef_2(Crystal_next(q + 1))*S_cent(Crystal_next(q + 1))^2/(2*D_nd) + ((Coef_1(Crystal_next(q + 1))*S_left(Crystal_next(q + 1))^2 + 2*dCpodX_left(Crystal_next(q + 1))*D_nd + 2*Coef_2(Crystal_next(q + 1))*S_left(Crystal_next(q + 1)))*S_cent(Crystal_next(q + 1)))/(2*D_nd) + (-2*Coef_1(Crystal_next(q + 1))*S_left(Crystal_next(q + 1))^3 - 6*D_nd*S_left(Crystal_next(q + 1))*dCpodX_left(Crystal_next(q + 1)) - 3*Coef_2(Crystal_next(q + 1))*S_left(Crystal_next(q + 1))^2 + 6*C_sat*D_nd)/(6*D_nd);
            C(max_num_array(Crystal_next(q)))    = -Coef_1(Crystal_next(q))*S_cent(Crystal_next(q))^3/(6*D_nd) - Coef_2(Crystal_next(q))*S_cent(Crystal_next(q))^2/(2*D_nd) + ((Coef_1(Crystal_next(q))*S_left(Crystal_next(q))^2 + 2*dCpodX_left(Crystal_next(q))*D_nd + 2*Coef_2(Crystal_next(q))*S_left(Crystal_next(q)))*S_cent(Crystal_next(q)))/(2*D_nd) + (-2*Coef_1(Crystal_next(q))*S_left(Crystal_next(q))^3 - 6*D_nd*S_left(Crystal_next(q))*dCpodX_left(Crystal_next(q)) - 3*Coef_2(Crystal_next(q))*S_left(Crystal_next(q))^2 + 6*C_sat*D_nd)/(6*D_nd);
            C(max_num_array(Crystal_next(q+1)))  = -Coef_1(Crystal_next(q + 1))*S_cent(Crystal_next(q + 1))^3/(6*D_nd) - Coef_2(Crystal_next(q + 1))*S_cent(Crystal_next(q + 1))^2/(2*D_nd) + ((Coef_1(Crystal_next(q + 1))*S_left(Crystal_next(q + 1))^2 + 2*dCpodX_left(Crystal_next(q + 1))*D_nd + 2*Coef_2(Crystal_next(q + 1))*S_left(Crystal_next(q + 1)))*S_cent(Crystal_next(q + 1)))/(2*D_nd) + (-2*Coef_1(Crystal_next(q + 1))*S_left(Crystal_next(q + 1))^3 - 6*D_nd*S_left(Crystal_next(q + 1))*dCpodX_left(Crystal_next(q + 1)) - 3*Coef_2(Crystal_next(q + 1))*S_left(Crystal_next(q + 1))^2 + 6*C_sat*D_nd)/(6*D_nd);
       for i=max_num_array(Crystal_next(q))+1:max_num_array(Crystal_next(q+1))-1
           FLAG = 0;
           if (x_nd(i)>=S_left(Crystal_next(j))-stability_coef_bound*h_nd) && (x_nd(i)<=S_right(Crystal_next(j))+stability_coef_bound*h_nd)   
                for m=1:n
                    g_s(m) = x_nd(i)*Coef_1(m)+Coef_2(m);
                end
                    %g_s(Crystal_next(j)) = x_nd(i)*Coef_1(m)+Coef_2(m);
                    g      = x_nd(i)*Coef_1(Crystal_next(j))+Coef_2(Crystal_next(j));
                    FLAG   = 1;          
                     
                if (x_nd(i)<S_left(Crystal_next(j)) && x_nd(i)>=S_left(Crystal_next(j))-stability_coef_bound*h_nd)
                    g      = g_s(Crystal_next(j))*(x_nd(i)-S_left(Crystal_next(j))-stability_coef_bound*h_nd)/stability_coef_bound*h_nd;    
                end
                if(x_nd(i)>S_right(Crystal_next(j)) && x_nd(i)<=S_right(Crystal_next(j))+stability_coef_bound*h_nd)
                    g      = g_s(Crystal_next(j))*(S_right(Crystal_next(j))+stability_coef_bound*h_nd-x_nd(i))/stability_coef_bound*h_nd;
                end 
                    FLAG_1 = 1;
           else
                    
                    Proc_Cryst = 0.2;
                    g          = C_bound/(1-Proc_Cryst)-C_bound;
                    %g          = 0;
                    %g          =;
                    %g         = 0;
                if(j~=N_cryst&&FLAG_1==1)
                    j      = j+1;
                    FLAG_1 = 0;
                end
           end
            h_befor = x_nd(i)-x_nd(i-1);
            h_after = x_nd(i+1)-x_nd(i);
            ci      = (D_nd*(2/(h_after+h_befor)))/h_after;
            bi      = -((h_befor+h_after)*(2/(h_after+h_befor))*D_nd/(h_befor*h_after)+1/tau_nd);
            ai      = (2/(h_after+h_befor))*D_nd/(h_befor);
            fi      = -g-C(i)/tau_nd;
            
            alfa(i) = -ci/(bi+ai*alfa(i-1)); % прогоночные коэффициенты
            beta(i) = -(ai*beta(i-1)-fi)/(bi+ai*alfa(i-1));% прогоночные коэффициенты  
       end
        
        for k       = max_num_array(Crystal_next(q+1))-1:-1:max_num_array(Crystal_next(q))+1
            C(k)    = alfa(k)*C(k+1)+ beta(k);
        end 
        
     end
      
     if (time>Drob_figuer_counter*f||(time>1-1.e-8&&time<1+1.e-8))
          
      str_Time           = {'Time(years)=',num2str(time*t_end_yrs)};
      str_Crystal_Numb   = {'Crystal Number=',num2str(n)};
  
      str_C_limit        = {'C_{limit}=',num2str(C_limit)};
      str_L              = {'L=',num2str(L)};
      str_Coef_nucl_c    = {'J_0=',num2str(J_0)};
      str_Coef_nucl_T    = {'A=',num2str(A)};
      str_Coef_nucl_mult = {'B=',num2str(B)};
      str_Coef_N         = {'N',num2str(Coef_N)};
      str_Coef_delta_C   = {'\Delta C_{*}',num2str(Coef_delta_C)};
      set(Graph,'CurrentAxes',ha);
      
      x_n=x_nd*C_bound;
      
      plot(x_nd,C);
 
      
      set(Graph,'CurrentAxes',axes_raspr);
      %set(axes_raspr, 'YScale', 'log')
      switch (Spec)
          case (1)
      plot(HRT_C_raspr(1,:),HRT_C_raspr(2,:),Crystal_size_x,log(Crystal_size_y/Number_crystals_max),'o');
      legend('HRT-C','Calculated');
          case(2)
      plot(LV_13(1,:),LV_13(2,:),Crystal_size_x,log(Crystal_size_y/Number_crystals_max),'o');
      legend('LV-13','Calculated');
          case(3)
      plot(MFT_1_raspr(1,:),MFT_1_raspr(2,:),Crystal_size_x,log(Crystal_size_y/Number_crystals_max),'o');
      legend('MFT-1','Calculated');
          case(4)
      plot(HRT_A_raspr(1,:),HRT_A_raspr(2,:),Crystal_size_x,log(Crystal_size_y/Number_crystals_max),'o');
      legend('HRT-A','Calculated');
          case(5)
      plot(LV_3(1,:),LV_3(2,:),Crystal_size_x,log(Crystal_size_y/Number_crystals_max),'o');
      legend('LV-3','Calculated');
          case(6)
      plot(TM_15(1,:),TM_15(2,:),Crystal_size_x,log(Crystal_size_y/Number_crystals_max),'o');
      legend('MT-15','Calculated');
      end
      
      set(axes_raspr,'xlim',[0 1]);
      set(axes_raspr,'ylim',[-8 1]);
      
      %set(axes_raspr,'xlim',[0 1]);
      %set(axes_raspr,'ylim',[-6 1]);
      
      text(1.05,-11,str_Time);
      text(1.35,-11,str_Crystal_Numb);
      text(1.75,-11,str_C_limit);
      text(2.15,-11,str_L);
      
      text(0,-11.0,str_Coef_nucl_c);
      text(0.2,-11,str_Coef_nucl_T);
      text(0.4,-11,str_Coef_nucl_mult);
      %text(0.6,-11,str_Coef_N);
      %text(0.8,-11,str_Coef_delta_C);
      
      text(1.05,24.5,str_Speciment);
      
      set(axes_raspr,'xlim',[0 1]);
      set(axes_raspr,'ylim',[-8 1]);
      
      set(Graph,'CurrentAxes',axes_vel)
      
      if(FLAG_log_graph)
      plot(ti_nd,log(V_Nucleation_Array));      
      set(axes_vel,'xlim',[0 t_end_yrs]);
      set(axes_vel,'ylim',[-3 15]);
      else
      plot(ti_nd,V_Nucleation_Array);      
      set(axes_vel,'xlim',[0 t_end_yrs]);
      %set(axes_vel,'ylim',[-3 15]);
      end
      
      set(Graph,'CurrentAxes',axes_tempr)
      plot(ti_nd,T_Array);
      set(axes_tempr,'xlim',[0 t_end_yrs]);
      set(axes_tempr,'ylim',[T_end-273 T_start-273]);
      xlabel(axes_tempr,'Time (years)')
      ylabel(axes_tempr,'Temperature T, ^oC')
      
      set(Graph,'CurrentAxes',axes_V)
      %Coef_Nucl_mult*exp((abs(C(i)-C_sat))*Coef_nucl_c)*exp(-1/T*(Coef_nucl_T));
      plot(x_nd,V_nucleation(J_0,A,B,T,C-C_sat));
      
     
   
      xlabel(ha,'x')
      ylabel(ha,'C')     
      
      %xlabel(axes_raspr,'D/D_{max}')
     % ylabel(axes_raspr,'ln(N/N_{max})') 
      
      xlabel(axes_vel,'Time(years)')
      ylabel(axes_vel,'J') 
      
      xlabel(axes_V,'\Delta C')
      ylabel(axes_V,'J')
      
      drawnow;
      f = f+1;

     end 
     
    if(time>=raspr_time(nn))
    raspr_prom(nn,:) = Crystal_size_y(1:end);
    nn               = nn+1;
    end
     
    end

   raspr(cef,1:end)           = Crystal_size_y(1:end);
   Crystal_size_global(cef,:) = Crystal_size(:);
   
end
raspr(21,1:end)  = y_teor_drob;
toc
Crystalwe        = zeros(2,Drob);
Crystalwe(1,:)   = Crystal_size_y(1:Drob);
for i=0:Drob-1
Crystalwe(2,i+1) = Crystal_size_max*i/10;
end
C_sat_end         = C_sat;
%C_kek=C(C>C_bound)
delta_C = C_sat_start-C_sat_end
%close all