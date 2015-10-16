function []=DYNvsEpsilon1
clf;
% automatically create postscript whenever
% figure is drawn
tmpfilename = 'DormFrac_max_host_demography_vary_delta_clearance';
%tmpfilename = 'DormFrac_max_host_demography_vary_delta';
%tmpfilename = 'DormFrac_max_host_demography_vary_delta_no_cell_death';
%tmpfilename = 'DormFrac_max_host_demography_vary_delta_no_cell_death_no_virus_decay';
%tmpfilename = 'DormFrac_max_host_demography_vary_delta_no_cell_demography_no_virus_decay';
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

  
phi=10^(-9); %cells/(ml*hr)
S0=8.3*10^8;
%V0=S0*0.04;
%delta=45;
mu=1/48
p=1/36
d=0.0866

r=0.3390%0.001
K=8.9470e+08%10^(10)

MOIvec=linspace(0.01,0.05,20); %1/24 approximately 0.042
deltavec=linspace(25,100,5)
Dvec=zeros(length(MOIvec),length(0:0.1:48));
Nvec=zeros(length(MOIvec),length(0:0.1:48));
%Dvec=zeros(1,length(0:0.1:48));
%Nvec=zeros(1,length(0:0.1:48));




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

for j=1:length(deltavec)
        delta=deltavec(j);
        
for i=1:length(MOIvec)
        MOI=MOIvec(i);
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

end


tmph=plot(MOIvec,maxfrac1(1,:),':m*',MOIvec,maxfrac1(2,:),':b*',MOIvec,maxfrac1(3,:),':g*',MOIvec,maxfrac1(4,:),':c*',MOIvec,maxfrac1(5,:),':r*');
legend({'\delta=25','\delta=43.75 ','\delta= 62.5','\delta=81.25 ','\delta= 100'},'Position',[0.62,0.35,0.25,0.1],'FontSize',20,'FontWeight','bold')

set(tmph,'linewidth',2);

%tmph=plot(pvec,maxfrac,'ko');
%set(tmph,'markerfacecolor','k','markersize',14);
%tmpt=text(pvec,maxfrac,'(pvec,maxfrac)');
%set(tmpt,'fontsize',14,'interpreter','latex');

set(gca,'xtick');
set(gca,'ytick');
set(gca,'xticklabel');
set(gca,'yticklabel');
xlabel('Initial virus-host ratio, $V_0/S_0$','fontsize',20,'verticalalignment','top','interpreter','latex');
ylabel('Maximum dormant cell ratio, $[D(t)/N(t)]_{max}$','fontsize',20,'verticalalignment','bottom','interpreter','latex');
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


