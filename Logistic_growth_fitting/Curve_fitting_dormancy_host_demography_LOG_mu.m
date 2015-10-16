function []=HostParamFit

%
%loadstruct=load('Data002ST.mat');
%Ndata=loadstruct.Data002ST;
Ndata=[12,8.8250*10^7.2;24,6.9750*10^7.2;48,5.8840*10^7.2]
    %;72,6.4949*10^7]; % [time, data;time,data;time,data]

data=[Ndata];
tdat=data(:,1)
ydat=data(:,2);
%ydat1=data(:,2)
%tdat=round(tdat);
length(tdat)
length(ydat)
phi=10^(-6); %cells/(ml*hr)
S0=6*10^7.9;
V0=S0*0.04;
delta=45;
%mu=1/24
d=0.0886
gamma=1/72

r=0.3390%0.001
K=8.9470e+08%10^(10)

%%
function f = funsys(t,X,p) 
f=zeros(4,1);
mu=p(1);


S1 =X(1);
D1=X(2);
I1=X(3);
V1=X(4);


f(1)= -phi*S1*V1*(1+delta)+gamma*(D1+I1)+r*S1*(1-(S1+D1+I1)/K);
f(2)=phi*S1*V1*delta-gamma*D1;
f(3)=phi*S1*V1-gamma*I1-mu*I1;
f(4)=-phi*S1*V1-d*V1;

end

 
function[z] = Solve_mybird(param,tdat)
   
      
    init_cond = [S0,0,0,V0];
    %tspan=min(tdat):max(tdat);
    [T,x] = ode45(@(t,X) funsys(t,X,param),tdat,init_cond);
    N1=x(:,1)+x(:,2)+x(:,3);
    %z1=P()
    z=N1;

end
paramguess=[1/72];

initi_cond =[S0,0,0,V0];
tdat1=0:0.1:48;
[t1,x1] = ode45(@(t,X) funsys(t,X,paramguess),tdat1,initi_cond);
u1=x1(:,1)+x1(:,2)+x1(:,3);

figure (2)
plot(t1,u1,'r-');

%%
[paramfit,resnorm] = lsqcurvefit(@Solve_mybird,paramguess,tdat,ydat,[0.001],[100]);
paramfit(1)
resnorm

initi_cond = [S0,0,0,V0];
tdat1=0:0.1:48;
[t1,x1] = ode45(@(t,X) funsys(t,X,paramfit),tdat1,initi_cond);
u1=x1(:,1)+x1(:,2)+x1(:,3);



figure (3)
plot(t1,u1,'r-',Ndata(:,1),Ndata(:,2),'o');


%figure (4)
%plot(t1,log10(u1),'r-',tdat,ydat,'o');

end