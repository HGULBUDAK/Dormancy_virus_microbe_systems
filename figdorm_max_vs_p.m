function []=DYNvsEpsilon1
clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'DormFrac_max_vs_p';
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

  
phi=10^(-7); %cells/(ml*hr)
S0=10^7;
V0=S0*0.022;
delta=45;

r=1
K=10^(10)

pvec=linspace(0,1,20); %1/24 approximately 0.042
Dvec=zeros(length(pvec),length(0:0.1:24))


    function f = funsys(t,X,p)
f=zeros(4,1); 
S1 =X(1);
D1=X(2);
I1=X(3);
V1=X(4);

f(1)= -phi*S1*V1*(1+delta)+p*(D1+I1)+r*S1*(1-(S1+D1+I1)/K);
f(2)=phi*S1*V1*delta-p*D1;
f(3)=phi*S1*V1-p*I1;
f(4)=-phi*S1*V1;
end

for i=1:length(pvec)
        p=pvec(i);

[T,Xv]=ode45(@(t,X)funsys(t,X,p),0:0.1:24,[S0;0;0;V0]); 
S=Xv(:,1);
D=Xv(:,2);
I=Xv(:,3);
V=Xv(:,4);
T;
TL=length(T);

Dvec(i,:)=D

[maxfrac(i),inx(i)]=max(Dvec(i,:)./S0)
T1(i)=T(inx(i))

end
T(2)


tmph=plot(pvec,T1,'b:*')
%,T,Dvec(2,:)./S0,'r--',T,Dvec(3,:)./S0,'r:',T,Dvec(4,:)./S0,'r-.',T,Dvec(5,:)./S0,'r--.')
set(tmph,'linewidth',2);

%tmph=plot(pvec,maxfrac,'ko');
%set(tmph,'markerfacecolor','k','markersize',14);
%tmpt=text(pvec,maxfrac,'(pvec,maxfrac)');
%set(tmpt,'fontsize',14,'interpreter','latex');

set(gca,'xtick');
set(gca,'ytick');
set(gca,'xticklabel');
set(gca,'yticklabel');
xlabel('Clearance rate, $p$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('$D_{max}(t)/S_0$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
%title({'Fraction of max. dormant cell density to',' initial susceptible cell density,  $D_{max}/S_0$},'fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
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


