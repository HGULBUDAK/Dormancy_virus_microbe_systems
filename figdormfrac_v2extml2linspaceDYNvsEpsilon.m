function []=DYNvsEpsilon
clf;
clear all
set(0,'DefaultTextInterpreter', 'latex')
tmpfilename = 'figdormfrac_v2extml2linspaceDYNvsEpsilon';
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
delta=45;
epsilonvec=linspace(0,1,101);

intVvec=zeros(3,length(epsilonvec));
Dtvec=zeros(3,length(epsilonvec));
Stvec=zeros(3,length(epsilonvec));
Itvec=zeros(3,length(epsilonvec));
Vtvec=zeros(3,length(epsilonvec));

    function f = funsys(t,X,epsilon)
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

for i=1:length(epsilonvec)
        epsilon=epsilonvec(i);

[T,Xv]=ode45(@(t,X)funsys(t,X,epsilon),0:0.1:72,[S0;0;0;V0]); 
S=Xv(:,1);
D=Xv(:,2);
I=Xv(:,3);
V=Xv(:,4);
T;
TL=length(T);
ind=[find(T==24),find(T==48),find(T==72)];



S(length(T));
for j=1:length(ind)
%intVvec(j,i)=trapz(T(1:ind(j)),V(1:ind(j)));
%Stvec(j,i)=S0-((1+delta)/(1-(epsilon.*delta)/(1+delta))).*V0+((phi.*epsilon*delta*S0)./((1-(epsilon.*delta)/(1+delta)))).*intVvec(j,i);
%Dtvec(j,i)=(delta/(1+delta)).*(S0-Stvec(j,i));
%Itvec(j,i)=(1/(1+delta)).*(S0-Stvec(j,i));

Stvec(j,i)=S(ind(j));
Dtvec(j,i)=D(ind(j));
Itvec(j,i)=I(ind(j));
Vtvec(j,i)=V(ind(j));

end

end



figure(1)
tmph=plot(epsilonvec,Stvec(1,:),'r-*',epsilonvec,Stvec(2,:),'g-*',epsilonvec,Stvec(3,:),'b-*')
set(tmph,'linewidth',2);
xlabel('Reduction in the dormant cell virus absorption, $\epsilon$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Susceptible cell density, $S(t)$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
title('Susceptible cell density $S(t)$ at $t=24,48,72$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
%set(gca,'fontsize',20);
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
tmpt = pwd;
%clear tmp*

figure(2)
plot(epsilonvec,Dtvec(1,:),'r-*',epsilonvec,Dtvec(2,:),'g-*',epsilonvec,Dtvec(3,:),'b-*')
xlabel('Reduction in the dormant cell virus absorption, $\epsilon$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Dormant cell density, $D(t)$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
title('Dormant cell density $D(t)$ at $t=24,48,72$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
set(gca,'fontsize',20);
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
tmpt = pwd;
%clear tmp*


figure(3)
plot(epsilonvec,Itvec(1,:),'r-*',epsilonvec,Itvec(2,:),'g-*',epsilonvec,Itvec(3,:),'b-*')
xlabel('Reduction in the dormant cell virus absorption, $\epsilon$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Infected cell density, $I(t)$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
title('Infected cell density $I(t)$ at $t=24,48,72$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
set(gca,'fontsize',20);
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
tmpt = pwd;
%clear tmp*

figure(4)
plot(epsilonvec,Vtvec(1,:),'r-*',epsilonvec,Vtvec(2,:),'g-*',epsilonvec,Vtvec(3,:),'b-*')
xlabel('Reduction in the dormant cell virus absorption, $\epsilon$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Virus density, $V(t)$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
title('Virus density $V(t)$ at $t=24,48,72$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
set(gca,'fontsize',20);
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
tmpt = pwd;
%clear tmp*

figure(5)
plot(epsilonvec,Dtvec(1,:)/S0,'r-*',epsilonvec,Dtvec(2,:)/S0,'g-*',epsilonvec,Dtvec(3,:)/S0,'b-*')
xlabel('Reduction in the dormant cell virus absorption, $\epsilon$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Dormant cell density fraction, $D(t)/S_0$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
title('Dormant cell density fraction $D(t)/S_0$ at $t=24,48,72$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
set(gca,'fontsize',20);
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
tmpt = pwd;
%clear tmp*

figure(6)
plot(epsilonvec,Dtvec(1,:)/V0,'r-*',epsilonvec,Dtvec(2,:)/V0,'g-*',epsilonvec,Dtvec(3,:)/V0,'b-*')
xlabel('Reduction in the dormant cell virus absorption, $\epsilon$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel({'Population dormancy enhancement,', '$D(t)/V_0$'},'fontsize',20,'verticalalignment','bottom','interpreter','latex');
title('Population dormancy enhancement, $D(t)/V_0$ at $t=24,48,72$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
set(gca,'fontsize',20);
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
tmpt = pwd;
clear tmp*

end


