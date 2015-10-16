function []=WHostParamFit8
clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'FigCurve_fitting_dormancy_host_demography';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[216 41 606 756]);
tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

%%
%loadstruct=load('Data0011.mat');
%Sdata=loadstruct.Data0011;
Sdata=[12,8.9099*10^8;24,8.9401*10^8;48,8.9448*10^8;72,8.9606*10^8]; % [time, data;time,data;time,data]
data=[Sdata];
tdat=data(:,1)
ydat=log10(data(:,2))
%tdat=round(tdat);
length(tdat)
length(ydat)


%%
function f = funsys(t,X,p) 
f=zeros(1,1);
r0=p(1);
K0=p(2);
mu=p(3);

S1 =X(1);

f(1)= r0*S1*(1-(S1/K0))-mu*S1;
end

 
function[z] = Solve_mybird(param,tdat)
   
      
    init_cond = [8.3*10^8];
    %tspan=min(tdat):max(tdat);
    [T,x] = ode15s(@(t,X) funsys(t,X,param),tdat,init_cond);
    P=x(:,1);
    %z1=P()
    z=log10(P);
    
end
paramguess=[0.239,(10^8)*(9.4155),1/24];

initi_cond =[8.3*10^8];
tdat1=0:0.1:75;
[t1,x1] = ode15s(@(t,X) funsys(t,X,paramguess),tdat1,initi_cond);
u1=x1(:,1);

%figure (2)
%plot(t1,u1,'r-');

%%
[paramfit,resnorm] = lsqcurvefit(@Solve_mybird,paramguess,tdat,ydat,[0.0001,ydat(end),0.0001],[0.5,10^10,0.5]);
paramfit(1)
paramfit(2)
paramfit(3)
resnorm

initi_cond = [8.3*10^8]
tdat1=0:0.1:75;
[t1,x1] = ode15s(@(t,X) funsys(t,X,paramfit),tdat1,initi_cond);
u1=x1(:,1);



%figure (3)
%plot(t1,u1,'r-',Sdata(:,1),Sdata(:,2),'o');


%figure (4)
%plot(t1,log10(u1),'r-',tdat,ydat,'o');

tmph=plot(t1,u1,'r-',Sdata(:,1),Sdata(:,2),'o')
set(tmph,'linewidth',2);
set(gca,'xtick');
set(gca,'ytick');
set(gca,'xticklabel');
set(gca,'yticklabel');
xlabel('Hours postinfection, $t$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Total cell density, $N(t)$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
title('Total cell density $N(t)$ at $p=0,0.5,1$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
set(gca,'fontsize',20);
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
psprintc(tmpfilenoname);
psprint(tmpfilebwname);
tmpt = pwd;
tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
%datenamer(1.1,.05,90);
psprintc(tmpfilename);
clear tmp*



end