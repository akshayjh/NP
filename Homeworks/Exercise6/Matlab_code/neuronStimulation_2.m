clear all;
close all;
clc

% Eigene Variablen
Strom = 2; %mA
Phase = 2; % 1:Mono-phase 2:Bi-phase

C=diag(-2*ones(100,1),0)+diag(ones(99,1),1)+diag(ones(99,1),-1);
C(1,1)=-1; C(100,100)=-1;

%Constants set
Cm=1; % Membrane Capcitance uF/cm^2
dt=0.025; % Time Step ms
t=0:dt:30; %Time Array ms

ENa=115; % mv Na reversal potential
EK=-12; % mv K reversal potential
El=10.6; % mv Leakage reversal potential
gbarNa=120; % mS/cm^2 Na conductance
gbarK=36; % mS/cm^2 K conductance
gbarl=0.3; % mS/cm^2 Leakage conductance

k=3^(0.1*(6.3-6.3));

comp=100; %compartments

m=zeros(comp,length(t));
n=m;
h=m;

V=zeros(comp,length(t)); %in this case also start value V(:,1)=0;

IK=m;
Il=m;
INa=m;

gNa=m;
gK=m;
gl=m;

B=m;
X=m;
Ex1=m;

L=0.5*10^-4; %length compartment
r=1.5*10^-4;   %Radius Axon
roh=0.01; %spezifischer Widerstand
A=r^2*pi;
Ra=roh*(L/A);

%lineares Gleichungssystem A*x=B;
A=eye(100)-dt/(Cm*Ra)*C;

% für x=m
am=@(V)(2.5-0.1*V)./(exp(2.5-0.1*V)-1);
bm=@(V)4*exp(-V/18);

% für x=n
an=@(V)(0.1-0.01*V)./(exp(1-0.1*V)-1);
bn=@(V)0.125*exp(-V/80);

% für x=h
ah=@(V)0.07*exp(-V/20);
bh=@(V)1./(exp(3-0.1*V)+1);

m(:,1)=am(V(:,1))./(am(V(:,1))+ bm(V(:,1))); % Initial m-value
n(:,1)=an(V(:,1))./(an(V(:,1))+ bn(V(:,1))); % Initial n-value
h(:,1)=ah(V(:,1))./(ah(V(:,1))+ bh(V(:,1))); % Initial h-value

I = zeros(length(t)-1, 1);
% Externes potential berechnen
for i=1:length(t)-1
    for c=1:comp
        
        if i>=length(1:5/dt) && i<=length(1:6/dt)
            if Phase == 1
                Ex1(c,i)=(dt/(Cm*Ra)*300*(Strom)/(4*pi)*1/sqrt(0.001^2+((c-50)*0.5*10^-4)^2));
                I(i) = Strom;
            end
            if Phase == 2
                Ex1(c,i)=(dt/(Cm*Ra)*300*(-1*Strom)/(4*pi)*1/sqrt(0.001^2+((c-50)*0.5*10^-4)^2));
                I(i) = -Strom;
            end
        elseif i>length(1:6/dt) && i<=length(1:7/dt)
            if Phase>1
                Ex1(c,i)=(dt/(Cm*Ra)*300*(+1*Strom)/(4*pi)*1/sqrt(0.001^2+((c-50)*0.5*10^-4)^2));
                I(i) = Strom;
            end
        else
            Ex1(c,i)=0;
        end
    end
end

Ex=C*Ex1(:,:);
Ex(1,:)=Ex(2,:);
Ex(100,:)=Ex(99,:);


for i=1:length(t)-1
    for c=1:comp
        
        m(c,i+1)=m(c,i)*exp(-k*(am(V(c,i))+bm(V(c,i)))*dt)-am(V(c,i))/(am(V(c,i))+bm(V(c,i)))*(exp(-k*(am(V(c,i))+bm(V(c,i)))*dt)-1);
        n(c,i+1)=n(c,i)*exp(-k*(an(V(c,i))+bn(V(c,i)))*dt)-an(V(c,i))/(an(V(c,i))+bn(V(c,i)))*(exp(-k*(an(V(c,i))+bn(V(c,i)))*dt)-1);
        h(c,i+1)=h(c,i)*exp(-k*(ah(V(c,i))+bh(V(c,i)))*dt)-ah(V(c,i))/(ah(V(c,i))+bh(V(c,i)))*(exp(-k*(ah(V(c,i))+bh(V(c,i)))*dt)-1);
        
        gNa(c,i)=gbarNa*(m(c,i)^3)*h(c,i);
        gK(c,i)=gbarK*(n(c,i)^4);
        gl(c,i)=gbarl;
        
        INa(c,i)=gNa(c,i)*(V(c,i)-ENa);
        IK(c,i)=gK(c,i)*(V(c,i)-EK);
        Il(c,i)=gl(c,i)*(V(c,i)-El);
        
        B(c,i)=V(c,i)+(dt)*((1/Cm)*(-(INa(c,i)+IK(c,i)+Il(c,i))))+Ex(c,i);
        
    end
    V(:,i+1)=A\B(:,i); %lineares GLeichungssystem
    
end

V=V-70;
h=pcolor(V);
c = colorbar
set(h,'edgecolor','none')
ylabel('Compartment #');
xlabel('time in \itms')
ylabel(c, 'Action potential in mV')
if(Phase==1)
    title(['AP, mono-phasic current pulse with I=', num2str(Strom), 'mA']);
elseif(Phase==2)
    title(['AP, bi-phasic current pulse with I=', num2str(Strom), 'mA']);
else
    title('Error!');
end
grid on
set(gca,'FontSize',20);
% imagesc(V)
% [X,Y]=meshgrid(1:100,t);
% plot3(X,Y,V);
% xlabel('compartments c/1')
% ylabel('time t/ms')
% zlabel('voltage V/mV')
%quiver