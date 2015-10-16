function []=DYNvsEpsilon1
clf;
% automatically create postscript whenever
% figure is drawn
set(0,'DefaultTextInterpreter', 'latex')
tmpfilename = 'figdormfinal';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
%tmpfilenoname = sprintf('%s_noname',tmpfilename);

tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[680 170 690 636]);
tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

% main data goes here
mu=1/36
p=1/48
d=0.0866

r=0.3390%0.001
K=8.9470e+08%10^(10)
%info.p=0.5;  % cond prob of dormancy given release
%info.q=0.1;  % prob of infection given contact
%info.phi=10^-7; %cells/(ml*hr)
%info.kminus=10^6*60*60; % per hr
%info.kforw = info.kminus*info.q/(1-info.q);  % per hr
%info.kplus = info.phi*info.q;
phi=10^(-9); %cells/(ml*hr)
%info.S0=10^7;  % cells per ml
S0=8.3*10^8;

info.MOI=logspace(4,10,100);  % Viral density range
info.delta_range=logspace(-2,2,5);
info.Dt_range=linspace(-2,2,5);
info.Nt_range=linspace(-2,2,5);

function f = funsys(t,X,delta)
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

for j=1:length(info.delta_range),
delta=info.delta_range(j)
  % Assign the state
  %info.delta=info.delta_range(i);
  %info.Vcrit=info.S0/(1+info.delta);

  % Plot the Positive Omega case
  %tmpi=find(info.V0<info.Vcrit);
  %tmph = loglog(info.V0(tmpi)/info.S0,info.delta*info.V0(tmpi)/info.S0,'k-');
  %set(tmph,'linewidth',3);
  %hold on
  %tmpi=find(info.V0>=info.Vcrit);
  %tmph = loglog(info.V0(tmpi)/info.S0,ones(length(tmpi),1)*info.delta*info.S0/(1+info.delta)/info.S0,'k--');
  %set(tmph,'linewidth',3);
  %hold on
%end
for i=1:length(info.MOI),
    MOI=info.MOI;
    [T,Xv]=ode45(@(t,X)funsys(t,X,delta),0:0.1:48,[S0;0;0;S0.*MOI]); 
S=Xv(:,1);
D=Xv(:,2);
I=Xv(:,3);
V=Xv(:,4);
T;
TL=length(T);

Dvec(i,:)=D;
Nvec(i,:)=S+D+I;

[maxfrac(i),inx(i)]=max(Dvec(i,:)./Nvec(i,:));
T1(i)=T(inx(i));
maxfrac1(j,i)=maxfrac(i);
end
  % Assign the state
  %info.delta=info.delta_range(i);
  %info.Vcrit=info.S0/(1+info.delta);
  % Overlay the critical points
  %tmph=loglog(MOI,maxfrac,'rd');
  %set(tmph,'markerfacecolor','r','markersize',10);
  % Write the text
  %if (i>=3)
  %tmph=text(info.V0(1)/info.S0*0.2,2.2*info.delta*info.V0(1)/info.S0,sprintf('$\\delta = %d$',info.delta));
  %elseif (i==2) 
    %tmph=text(info.V0(1)/info.S0*0.2,2.2*info.delta*info.V0(1)/info.S0,sprintf('$\\delta = %4.1f$',info.delta));
  %else (i==1) 
   % tmph=text(info.V0(1)/info.S0*0.2,2.2*info.delta*info.V0(1)/info.S0,sprintf('$\\delta = %4.2f$',info.delta));
  end
  %set(tmph,'interpreter','latex','fontsize',14);
%end

% Label
tmph=plot(MOIvec,maxfrac1(1,:),'r:*',MOIvec,maxfrac1(2,:),'b:*',MOIvec,maxfrac1(3,:),'g:*',MOIvec,maxfrac1(4,:),'c:*',MOIvec,maxfrac1(5,:),'y:*');
set(tmph,'linewidth',2);
set(gca,'xtick');
set(gca,'ytick');
set(gca,'xticklabel');
set(gca,'yticklabel');
xlabel('Initial virus-host ratio, $V_0/S_0$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel({'Final dormant cell';'ratio, $D_{\infty}/S_0$'},'fontsize',20,'verticalalignment','bottom','interpreter','latex');
set(gca,'fontsize',20);
ylim([10^-6 5]);
xlim([10^-4 10^3.5]);
set(gca,'xtick',10.^[-3:1:3]);
set(gca,'ytick',10.^[-6:1:0]);

% Legend
%tmph=loglog([10^2 10^2.3]*10^-0.4,[10^-5 10^-5],'k--');
%set(tmph,'linewidth',2);
%tmph=loglog([10^2.15]*10^-0.4,[10^-4.7],'rd');
%set(tmph,'markerfacecolor','r','markersize',10);
%tmph=loglog([10^2 10^2.3]*10^-0.4,[10^-4.4 10^-4.4],'k-');
%set(tmph,'linewidth',2);
%tmpt=text(10^2,10^-5,'$\Omega<0$');
%set(tmpt,'interpreter','latex','fontsize',20);
%tmpt=text(10^2,10^-4.7,'$\Omega=0$');
%set(tmpt,'interpreter','latex','fontsize',20);
%tmpt=text(10^2,10^-4.4,'$\Omega>0$');
%set(tmpt,'interpreter','latex','fontsize',20);

% loglog(,, '');
%
%
% Some helpful plot commands
% tmph=plot(x,y,'ko');
% set(tmph,'markersize',10,'markerfacecolor,'k');
% tmph=plot(x,y,'k-');
% set(tmph,'linewidth',2);


% for use with layered plots
% set(gca,'box','off')

% adjust limits
% tmpv = axis;
% axis([]);
% ylim([]);
% xlim([]);

% change axis line width (default is 0.5)
% set(tmpa1,'linewidth',2)

% fix up tickmarks
% set(gca,'xtick',[1 100 10^4])
% set(gca,'ytick',[1 100 10^4])

% creation of postscript for papers
% psprint(tmpxfigfilename);

% the following will usually not be printed 
% in good copy for papers
% (except for legend without labels)

% legend
% tmplh = legend('stuff',...);
% tmplh = legend('','','');
% remove box
% set(tmplh,'visible','off')
% legend('boxoff');

% title('','fontsize',24)
% 'horizontalalignment','left');

% for writing over the top
% coordinates are normalized again to (0,1.0)
tmpa2 = axes('Position', tmppos);
set(tmpa2,'visible','off');
% first two points are normalized x, y positions
% text(,,'','Fontsize',14);

% automatic creation of postscript
% without name/date
%psprintc(tmpfilenoname);
psprint(tmpfilebwname);

tmpt = pwd;
tmpnamememo = sprintf('[source=%s/%s.ps]',tmpt,tmpprintname);
text(1.05,.05,tmpnamememo,'Fontsize',6,'rotation',90);
%datenamer(1.1,.05,90);
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