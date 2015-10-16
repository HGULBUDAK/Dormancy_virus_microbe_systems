function []=DYNvsEpsilon1
clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'Dormancy_Dyn_host_demography_N';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[216 41 606 756]);
tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

  
phi=10^(-9); %cells/(ml*hr)
S0=8.3*10^8;
V0=S0*0.04;
delta=45;
mu=1/24
d=0.0866

r=0.3390
K=8.9470e+08


pvec=linspace(0,1/24,4) %1/24 approximately 0.042
Svec=zeros(length(pvec),length(0:0.1:72));
Dvec=zeros(length(pvec),length(0:0.1:72));
Ivec=zeros(length(pvec),length(0:0.1:72));
Nvec=zeros(length(pvec),length(0:0.1:72));
    function f = funsys(t,X,p)
f=zeros(4,1); 
S1 =X(1);
D1=X(2);
I1=X(3);
V1=X(4);

f(1)= -phi*S1*V1*(1+delta)+p*(D1+I1)+r*S1*(1-(S1+D1+I1)/K);
f(2)=phi*S1*V1*delta-p*D1-mu*D1;
f(3)=phi*S1*V1-p*I1-mu*I1;
f(4)=-phi*S1*V1-d*V1;
end

for i=1:length(pvec)
        p=pvec(i);

[T,Xv]=ode45(@(t,X)funsys(t,X,p),0:0.1:72,[S0;0;0;V0]); 
S=Xv(:,1);
D=Xv(:,2);
I=Xv(:,3);
V=Xv(:,4);
T;
TL=length(T);

Svec(i,:)=S;
Dvec(i,:)=D;
Ivec(i,:)=I;
Nvec(i,:)=S+D+I;

end



tmph=plot(T,Nvec(1,:),'r-',T,Nvec(2,:),'r--',T,Nvec(3,:),'r:',T,Nvec(4,:),'r-.')
legend({'\gamma=0','\gamma\approx 0.014','\gamma \approx 0.028  ','\gamma \approx 0.042'},'Position',[0.65,0.7,0.25,0.1],'FontSize',15,'FontWeight','bold')
set(tmph,'linewidth',2);
set(gca,'xtick');
set(gca,'ytick');
set(gca,'xticklabel');
set(gca,'yticklabel');
xlabel('Hours postinfection (hpi), $t$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Total cell density, $N(t)$ (cells/(ml*hr))','fontsize',20,'verticalalignment','bottom','interpreter','latex');
%title('Total cell density $N(t)$ at varying $\gamma$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
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


