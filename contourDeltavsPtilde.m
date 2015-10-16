function []=CtourDeltavsEpsilon
clf;
clear all
set(0,'DefaultTextInterpreter', 'latex')
tmpfilename = 'figdormfrac_v2';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);
tmpxfigfilename = sprintf('x%s',tmpfilename);
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[12 138 763 655]);
tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

set(gcf,'DefaultLineMarkerSize',10);
phi=10^(-7); %cells/(ml*hr)
S0=10^7;
V0=S0*0.022;
epsilonvec=linspace(0,1,100);
intVvec=zeros(1,length(epsilonvec));
Sinfvec=zeros(1,length(epsilonvec));
Sinfsolvec=zeros(1,length(epsilonvec));

deltavec=linspace(0.01,50);
DinfS0vec=zeros(length(deltavec),length(epsilonvec));

function f = funsys(t,X,delta,epsilon)
f=zeros(4,1); 
S1 =X(1);
D1=X(2);
I1=X(3);
V1=X(4);

f(1)= -phi*S1*V1*(1+delta);
f(2)=phi*S1*V1*delta;
f(3)=phi*S1*V1;
f(4)=-phi*S1*V1-epsilon*phi*D1*V1;
end
for i=1:length(deltavec)
delta=deltavec(i);
for k=1:length(epsilonvec)
        epsilon=epsilonvec(k);

[T,Xv]=ode45(@(t,X)funsys(t,X,delta,epsilon),0:0.01:20,[S0;0;0;V0]); 
S=Xv(:,1);
D=Xv(:,2);
I=Xv(:,3);
V=Xv(:,4);
T
TL=length(T);
S(length(T))

intVvec(k,i)=trapz(T,V);
Sinfvec(k,i)=S0-((1+delta)/(1-(epsilon.*delta)/(1+delta))).*V0+((phi.*epsilon*delta*S0)./((1-(epsilon.*delta)/(1+delta)))).*intVvec(k,i);
Sinfsolvec(k,i)=S(length(T));

DinfvsS0vec(k,i)=(delta./(1+delta)).*(1-(Sinfvec(k,i)./S0));
end
end
DinfvsS0vec(length(epsilonvec),length(deltavec))
imagesc(deltavec,epsilonvec,flipud(DinfvsS0vec));
tmph=colorbar
%set(tmph,'yticklabel');
%hold on
[X,Y]=meshgrid(deltavec,epsilonvec);
contourf(X,Y,DinfvsS0vec)
%colormap(flipud(gray))
colorbar
%surf(X,Y,DinfvsS0vec)
%xlabel('Individual cell fate ratio,\delta');
xlabel('Individual cell fate ratio, $\delta$','fontsize',20,'verticalalignment','top','interpreter','latex');

% Create ylabel
%ylabel('Reduction in the dormant cell virus absorption,\epsilon');
ylabel('Reduction in the dormant cell virus absorption,$\epsilon$','fontsize',20,'verticalalignment','bottom','interpreter','latex');

% Create zlabel
%zlabel('The fraction of cells become dormant, D_\infty/S_0');
title('The fraction of cells become dormant, $D_{\infty}/S_0$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
set(gca,'fontsize',20);
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
tmpt = pwd;
clear tmp*
end


