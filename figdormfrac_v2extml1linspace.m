clf;
clear all
% automatically create postscript whenever
% figure is drawn
set(0,'DefaultTextInterpreter', 'latex')
tmpfilename = 'figdormfrac_v2extml1linspace';
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
set(gcf,'Position',[12 138 763 655]);
tmppos= [0.2 0.2 0.7 0.7];
tmpa1 = axes('position',tmppos);

% main data goes here
info.q=0.022;  % prob of infection given contact
info.phi=10^-7; %cells/(ml*hr)
info.kminus=10^6*60*60; % per hr
info.kforw = info.kminus*info.q/(1-info.q);  % per hr
info.kplus = info.phi*info.q;
info.S0=10^7;  % cells per ml
info.M0=(0.1+0.01)/2; %info.V0_range/info.S0; 
info.delta_range=linspace(0,50,1000);
info.ptilde_range=linspace(0,1,1000);

[vDelta,vPtilde]=meshgrid(info.delta_range,info.ptilde_range);

% 
for i=1:length(info.ptilde_range),
  info.ptilde=info.ptilde_range(i);
  for j=1:length(info.delta_range),
      info.delta=info.delta_range(j);
      info.OmegaCRIT=1-(1+info.delta)*info.M0;
     if (info.OmegaCRIT<0)
      Dout(i,j)=(info.delta+info.ptilde)/(1+info.delta); 
    else
      Dout(i,j)=(info.delta+info.ptilde)*info.M0;
    end
    DtoS(i,j) = Dout(i,j); % Dout(i,j)=D_\infty/S_0
  end
end
imagesc(info.delta_range,info.ptilde_range,flipud(DtoS));
tmph=colorbar
set(tmph,'yticklabel');
hold on
% Add a single contour at the critical value
[c,h]=contour(info.delta_range,info.ptilde_range,flipud(DtoS),[0 0]);
set(h,'linecolor','k','linewidth',3);

set(gca,'xtick',[0:5:50]);
set(gca,'ytick',[0:0.1:1]);
set(gca,'xticklabel',{'0';'5';'10';'15';'20';'25';'30';'35';'40';'45';'50'});
set(gca,'yticklabel',{'1';'0.9';'0.8';'0.7';'0.6';'0.5';'0.4';'0.3';'0.2';'0.1';'0'});



% Label
xlabel('Individual cell fate ratio, $\delta$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel({'Probability of host being dormant';'after viral entry, $\tilde{p}$'},'fontsize',20,'verticalalignment','bottom','interpreter','latex');
title('Fraction of final size of dormant cells, $D_{\infty}/S_0$','fontsize',20,'verticalalignment','bottom','interpreter','latex','horizontalalignment','center');
set(gca,'fontsize',20);



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
 %set(tmpa1,'linewidth',2)

% fix up tickmarks
% set(gca,'xtick',[1 100 10^4])
% set(gca,'ytick',[1 100 10^4])

% creation of postscript for papers
 psprint(tmpxfigfilename);

% the following will usually not be printed 
% in good copy for papers
% (except for legend without labels)

% legend
% tmplh = legend('stuff',...);
% tmplh = legend('','','');
% remove box
% set(tmplh,'visible','off')
% legend('boxoff');

 %title('','fontsize',24)
 %'horizontalalignment','left');

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

