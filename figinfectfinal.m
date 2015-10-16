clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figinfectfinal';
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
set(gcf,'Position',[680 170 690 636]);
tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

% main data goes here
info.p=0.5;  % cond prob of dormancy given release
info.q=0.1;  % prob of infection given contact
info.phi=10^-7; %cells/(ml*hr)
info.kminus=10^6*60*60; % per hr
info.kforw = info.kminus*info.q/(1-info.q);  % per hr
info.kplus = info.phi*info.q;
info.S0=10^7;  % cells per ml
info.V0=logspace(4,10,100);  % Viral density range
info.delta_range=logspace(-2,2,5);

% 
for i=1:length(info.delta_range),

  % Assign the state
  info.delta=info.delta_range(i);
  info.Vcrit=info.S0/(1+info.delta);

  % Plot the Positive Omega case
  tmpi=find(info.V0<info.Vcrit);
  tmph = loglog(info.V0(tmpi)/info.S0,info.V0(tmpi)/info.S0,'k-');
  set(tmph,'linewidth',3);
  hold on
  tmpi=find(info.V0>=info.Vcrit);
  tmph = loglog(info.V0(tmpi)/info.S0,ones(length(tmpi),1)*info.S0/(1+info.delta)/info.S0,'k--');
  set(tmph,'linewidth',3);
  hold on
end
for i=1:length(info.delta_range),
  % Assign the state
  info.delta=info.delta_range(i);
  info.Vcrit=info.S0/(1+info.delta);
  % Overlay the critical points
  tmph=loglog(info.Vcrit/info.S0,1/(1+info.delta),'rd');
  set(tmph,'markerfacecolor','r','markersize',10);
  % Write the text
  if (i>=3)
  tmph=text(info.V0(end)/info.S0*1.3,1/(1+info.delta),sprintf('$\\delta = %d$',info.delta));
  elseif (i==2) 
    tmph=text(info.V0(end)/info.S0*1.3,0.9/(1+info.delta),sprintf('$\\delta = %4.1f$',info.delta));
  else (i==1) 
    tmph=text(info.V0(end)/info.S0*1.3,1.2/(1+info.delta),sprintf('$\\delta = %4.2f$',info.delta));
  end
  set(tmph,'interpreter','latex','fontsize',14);
end

% Label
xlabel('Initial virus-host ratio, $V_0/S_0$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel({'Final infected cell';'ratio, $I_{\infty}/S_0$'},'fontsize',20,'verticalalignment','bottom','interpreter','latex');
set(gca,'fontsize',20);
ylim([10^-4 5]);
xlim([10^-3.5 10^4.7]);
set(gca,'xtick',10.^[-3:1:3]);
set(gca,'ytick',10.^[-6:1:0]);


% Legend
tmph=loglog([10^2 10^2.3]*10^-0.4,[10^-3.5 10^-3.5],'k--');
set(tmph,'linewidth',2);
tmph=loglog([10^2.15]*10^-0.4,[10^-3.2],'rd');
set(tmph,'markerfacecolor','r','markersize',10);
tmph=loglog([10^2 10^2.3]*10^-0.4,[10^-2.9 10^-2.9],'k-');
set(tmph,'linewidth',2);
tmpt=text(10^2,10^-3.5,'$\Omega<0$');
set(tmpt,'interpreter','latex','fontsize',20);
tmpt=text(10^2,10^-3.2,'$\Omega=0$');
set(tmpt,'interpreter','latex','fontsize',20);
tmpt=text(10^2,10^-2.9,'$\Omega>0$');
set(tmpt,'interpreter','latex','fontsize',20);

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
psprintc(tmpfilenoname);
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
