function []=WHostParamFit8

%%
loadstruct=load('Data001.mat');
Sdata=loadstruct.Data001;

data=[Sdata];
tdat=data(:,1)
ydat=data(:,2)
%tdat=round(tdat);
length(tdat)
length(ydat)


%%
function f = funsys(t,X,p) 
f=zeros(1,1);
r=p(1);
K=p(2);

S1 =X(1);

f(1)= r*S1*(1-(S1/K));
end

 
function[z] = Solve_mybird(param,tdat)
   
      
    init_cond = [10^7];
    %tspan=min(tdat):max(tdat);
    [T,x] = ode15s(@(t,X) funsys(t,X,param),tdat,init_cond);
    P=x(:,1);
    %z1=P()
    z=P;
    
end
paramguess=[0.00339,(10^8)*(9.4155)];

initi_cond =[10^7];
tdat1=0:0.1:75;
[t1,x1] = ode15s(@(t,X) funsys(t,X,paramguess),tdat1,initi_cond);
u1=x1(:,1);

figure (2)
plot(t1,u1,'r-');

%%
[paramfit,resnorm] = lsqcurvefit(@Solve_mybird,paramguess,tdat,ydat,[0.0001,ydat(end)],[10,10^10]);
paramfit(1)
paramfit(2)
resnorm
ydat(1)
ydat(2)

initi_cond = [10^7]
tdat1=0:0.1:75;
[t1,x1] = ode15s(@(t,X) funsys(t,X,paramfit),tdat1,initi_cond);
u1=x1(:,1);



figure (3)
plot(t1,u1,'r-',Sdata(:,1),Sdata(:,2),'o');

end