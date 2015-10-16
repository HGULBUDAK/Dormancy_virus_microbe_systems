function []=figdormfrac_v2extml2linspace1
clf;
clear all
set(0,'DefaultTextInterpreter', 'latex')
tmpfilename = 'figdormfrac_v2extml2linspace1';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[12 138 763 669]);
tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

% main data goes here
phi=10^(-7); %cells/(ml*hr)
S0=10^7;
V0=S0*0.022;
epsilonvec=linspace(0.001,1,100);
intVvec=zeros(1,length(epsilonvec));
Sinfvec=zeros(1,length(epsilonvec));
Sinfsolvec=zeros(1,length(epsilonvec));

deltavec=linspace(0.01,50,100);
DinfvsS0vec=zeros(length(deltavec),length(epsilonvec));

[vDelta,epsilon]=meshgrid(deltavec,epsilonvec);
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

[T,Xv]=ode45(@(t,X)funsys(t,X,delta,epsilon),0:0.01:24,[S0;0;0;V0]); 
S=Xv(:,1);
D=Xv(:,2);
I=Xv(:,3);
V=Xv(:,4);
T;
TL=length(T);
S(length(T));

intVvec(k,i)=trapz(T,V);
Sinfvec(k,i)=S0-((1+delta)/(1-(epsilon.*delta)/(1+delta))).*V0+((phi.*epsilon*delta*S0)./((1-(epsilon.*delta)/(1+delta)))).*intVvec(k,i);
Sinfsolvec(k,i)=S(length(T));

DinfvsS0vec(k,i)=(delta./(1+delta)).*(1-(Sinfsolvec(k,i)./S0));
end
end
DinfvsS0vec(length(epsilonvec),length(deltavec))
imagesc(deltavec,epsilonvec,flipud(DinfvsS0vec));
tmph=colorbar
set(tmph,'zticklabel',{'0';'0.1';'0.2';'0.3';'0.4';'0.5';'0.6';'0.7';'0.8';'0.9';'0.95'});
hold on
% Add a single contour at the critical value
[c,h]=contour(deltavec,epsilonvec,flipud(DinfvsS0vec),[0 0])
%set(h,'linecolor','k','linewidth',3);

set(gca,'xtick',[0.01:5:50.01]);
set(gca,'ytick',[0.001:0.1:1.001]);
set(gca,'xticklabel',{'0.01';'5';'10';'15';'20';'25';'30';'35';'40';'45';'50'});
set(gca,'yticklabel',{'1';'0.9';'0.8';'0.7';'0.6';'0.5';'0.4';'0.3';'0.2';'0.01';'0.001'});



% Label
xlabel('Individual cell fate ratio, $\delta$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel({'Reduction in dormant cell', 'virus absorption,$\epsilon$'},'fontsize',20,'verticalalignment','bottom','interpreter','latex');
title({'Fraction of dormant cells $24$ hrs after',' virus initiation $D(t)/S_0$'},'fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
set(gca,'fontsize',20);

% creation of postscript for papers
%psprint(tmpxfigfilename);



tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');

% automatic creation of postscript
% without name/date
psprint(tmpfilebwname);

tmpt = pwd;
tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
datenamer(1.1,.05,90);
% datename(.5,.05);
% datename2(.5,.05); % 2 rows

% automatic creation of postscript
psprintc(tmpfilename);

% set following on if zooming of 
% plots is required
% may need to get legend up as well
%axes(tmpa1)
%axes(tmplh)
clear tmp*

end


