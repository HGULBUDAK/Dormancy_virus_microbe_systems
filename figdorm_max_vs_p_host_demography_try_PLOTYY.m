function []=DYNvsEpsilon1  
phi=10^(-9); %cells/(ml*hr)
S0=8.3*10^8;
V0=S0*0.04;
delta=45;
mu=1/24
d=0.0866

r=0.3390%0.001
K=8.9470e+08%10^(10)

pvec=linspace(0,1/24,20); %1/24 approximately 0.042
Dvec=zeros(length(pvec),length(0:0.1:48));
Nvec=zeros(length(pvec),length(0:0.1:48));




    function f = funsys(t,X,p)
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

for i=1:length(pvec)
        p=pvec(i);

[T,Xv]=ode45(@(t,X)funsys(t,X,p),0:0.1:48,[S0;0;0;V0]); 
S=Xv(:,1);
D=Xv(:,2);
I=Xv(:,3);
V=Xv(:,4);
T;
TL=length(T);

Dvec(i,:)=D
Nvec(i,:)=S+D+I

[maxfrac(i),inx(i)]=max(Dvec(i,:)./Nvec(i,:))
T1(i)=T(inx(i))

end
T(2)

figure (1)
plotyy(pvec,maxfrac,pvec,T1)


end


