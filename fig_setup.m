f = figure('Units','normalized','Visible','off','Position',[0,0,1,1],'Color',[1 1 1]);

ha = axes('Units','normalized','Position',[0.05,.38,.40,0.57]);
xlim([0 1]); box on
ha1 = axes('Units','normalized','Position',[0.05,0.1,0.40,0.23]);
%xlim([0 1]); 
box on

ha2 = axes('Units','normalized','Position',[0.5,0.8,0.40,0.15]); box on
ha3 = axes('Units','normalized','Position',[0.5,0.56,0.40,0.15]); box on
ha4 = axes('Units','normalized','Position',[0.5,0.32,0.40,0.15]); box on
ha5 = axes('Units','normalized','Position',[0.5,0.1,0.40,0.15]); box on

xlabel(ha1,'Distance, um')
ylabel(ha1,'Zr concentration')
xlabel(ha2,'Time (years)')
ylabel(ha2,'Zr radius')
xlabel(ha3,'Time (years)')
ylabel(ha3,'Growth Rate, cm.s^{-1}')
xlabel(ha4,'Time (years)')
ylabel(ha4,'Temperature T, ^oC')
xlabel(ha5,'Zr radius')
ylabel(ha5,'Trace(norm)')


set(f,'CurrentAxes',ha)
set(f,'Name','DIFFUSOR')
movegui(f,'center')
set(f,'Visible','on');

hold on
tyear=3600*365.25*24;
