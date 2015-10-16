function []=HostParamFit

%
%loadstruct=load('Data002ST.mat');
%Ndata=loadstruct.Data002ST;
Ndata=[12,8.8250*10^7;24,6.9750*10^7;48,5.8840*10^7]
    %;72,6.4949*10^7]; % [time, data;time,data;time,data]

data=[Ndata];
tdat=data(:,1)
ydat=log10(data(:,2));
ydat1=data(:,2)
%tdat=round(tdat);
length(tdat)
length(ydat)
phi=10^(-7); %cells/(ml*hr)
S0=8.3*10^8;
V0=S0*0.04;
delta=45
%mu=1/24
d=0.0886
gamma=0

r=0.3390%0.001
K=8.9470e+08%10^(10)

%%
function f = funsys(t,X,p) 
f=zeros(1,1);
mu=p(1);


%S1 =X(1);
%D1=X(2);
%I1=X(3);
%V1=X(4);

N1=X(1);

%f(1)= -phi*S1*V1*(1+delta)+gamma*(D1+I1)+r*S1*(1-(S1+D1+I1)/K);
%f(2)=phi*S1*V1*delta-gamma*D1;
%f(3)=phi*S1*V1-gamma*I1-mu*I1;
%f(4)=-phi*S1*V1-d*V1;
r0=r+mu;
K0=K*(r0/r);
%g(0)=0;
%g(12)=0.4;
%g(24)=.95;
f(1)=r0*N1*(1-(0.04*t)*((delta+1)/delta))*(1-(N1/K0))-mu*N1*(1-(0.04*t));

end

 
function[z] = Solve_mybird(param,tdat)
   
      
    init_cond = [S0];
    %tspan=min(tdat):max(tdat);
    [T,x] = ode15s(@(t,X) funsys(t,X,param),tdat,init_cond);
    N1=x(:,1);
    %z1=P()
    z=log10(N1);

end
paramguess=[0.002];

initi_cond =[S0];
tdat1=0:0.1:48;
[t1,x1] = ode15s(@(t,X) funsys(t,X,paramguess),tdat1,initi_cond);
u1=x1(:,1);

figure (2)
plot(t1,u1,'r-');

%%
[paramfit,resnorm] = lsqcurvefit(@Solve_mybird,paramguess,tdat,ydat,[0.001],[1]);
guess=paramfit(1)
guess
resnorm

initi_cond = [S0];
tdat1=0:0.1:48;
[t1,x1] = ode15s(@(t,X) funsys(t,X,paramfit),tdat1,initi_cond);
u1=x1(:,1);



figure (3)
plot(t1,u1,'r-',Ndata(:,1),Ndata(:,2),'o');


figure (4)
plot(t1,log10(u1),'r-',tdat,ydat,'o');

end