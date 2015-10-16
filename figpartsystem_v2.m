clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'figpartsystem_v2';
tmpfilebwname = sprintf('%s_noname_bw',tmpfilename);
tmpfilenoname = sprintf('%s_noname',tmpfilename);

%tmpprintname = fixunderbar(tmpfilename);
% for use with xfig and pstex
tmpxfigfilename = sprintf('x%s',tmpfilename);

set(gcf,'DefaultLineMarkerSize',10);
% set(gcf,'DefaultLineMarkerEdgeColor','k');
% set(gcf,'DefaultLineMarkerFaceColor','w');
set(gcf,'DefaultAxesLineWidth',2);

set(gcf,'PaperPositionMode','auto');
set(gcf,'Position',[216 41 606 756]);

% main data goes here
info.p=0.5;  % cond prob of dormancy given release
info.q=0.1;  % prob of infection given contact
info.phi=10^-7; %cells/(ml*hr)
info.kminus=10^6*60*60; % per hr
info.kforw = info.kminus*info.q/(1-info.q);  % per hr
info.kplus = info.phi*info.q;
info.S0=10^7;  % cells per ml
info.delta=info.p*info.kminus/info.kforw;  % Dormant cells per infected host
info.Vcrit = info.S0/(1+info.delta);  % Critical value for the transition
info.V0_range = [0.1 1 1.9]*info.Vcrit;
info.trange=0:0.05:24;
for i=1:length(info.V0_range),

  % Set the plot area
  tmppos= [0.2 0.65-(i-1)*0.25 0.7 0.25];
  tmpa1 = axes('position',tmppos);
  
  % Simulate
  info.V0=info.V0_range(i);
  info.Omega = info.S0-(1+info.delta)*info.V0;
  y0 = [info.S0 info.V0];
  options=odeset('RelTol',1e-7);
  [t,y]=ode45(@vdormant_part,info.trange,y0,options,info);
  stats(i).info=info;
  stats(i).t=t;
  stats(i).y=y;
  if (info.Omega~=0)
    stats(i).Sanal = info.Omega./(1-(1+info.delta)*info.V0/info.S0*exp(-info.phi*info.Omega*t));
    stats(i).Vanal = info.V0/info.S0*info.Omega*exp(-info.phi*info.Omega*t)./(1-(1+info.delta)*info.V0/info.S0*exp(-info.phi*info.Omega*t));
  else
    stats(i).Vanal = info.V0./(info.V0*info.phi*(1+info.delta)*t+1);
    stats(i).Sanal = (1+info.delta)*stats(i).Vanal;
  end

  % Now plot
  tmph=semilogy(t,stats(i).Sanal,'b-');
  set(tmph,'linewidth',3);
  hold on
  tmph=semilogy(t,stats(i).Vanal,'r-');
  set(tmph,'linewidth',3);
  xlim([0 24]);
 
  % Display
  ymin=max(min(stats(i).Vanal(end),stats(i).Sanal(end)),10^-2);
  ymax=max(stats(i).Vanal(1),stats(i).Sanal(1));
  ylim([0.05*ymin 20*ymax]);
  if (i<3)  
    set(gca,'xtick',[]);
  else
    set(gca,'xtick',[0:5:25]);
    xlabel('Time, hrs','fontsize',20,'verticalalignment','top','interpreter','latex');
  end
  if (i~=2)
    set(gca,'ytick',10.^[-9:2:9]);
  else
    set(gca,'ytick',10.^[-9:1:9]);
  end
  set(gca,'fontsize',20);
  ylabel('Density, ml$^{-1}$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
  tmpyplot = [10^-1 10^4.5 10^-1];
  if (i~=2)
    tmph=text(2,tmpyplot(i),sprintf('$\\frac{\\Omega}{S_0} = %4.1f$',info.Omega/info.S0));
    set(tmph,'interpreter','latex','fontsize',20);
  else
    tmph=text(2,tmpyplot(i),sprintf('$\\frac{\\Omega}{S_0} = 0$'));
    set(tmph,'interpreter','latex','fontsize',20);
  end
 
  % Overlay the types
  tmpi = find(t==20);
  soffset=[0.25 3 8];
  voffset=[8 0.4 0.1];
  tmph=text(t(tmpi),y(tmpi,1)*soffset(i),'S');
  set(tmph,'interpreter','latex','fontsize',20);
  tmph=text(t(tmpi),y(tmpi,2)*voffset(i),'V');
  set(tmph,'interpreter','latex','fontsize',20);
end
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
