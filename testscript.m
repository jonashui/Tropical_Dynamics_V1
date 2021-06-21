close all
clear
clc

syms x y eps

%1
%  alpha=-2;
%  xm=y-exp(alpha/eps)*x^2;
%  ym=exp(-1/eps)*x;
%  uint=[0 4];
%  vint=[0 4];

%2
% beta=1;
% xm=y-x^2;
% ym=exp(-1/eps)*(exp(beta/eps)-x);
% uint=[-2 2];
% vint=[-2 2];

%3 autocad
a=0.5;
xm=exp((a-1)/eps)+exp(-1/eps)*(-x-x*y^2);
ym=-y+x+x*y^2;
uint=[-5 5];
vint=[-5 5];


%4 Kristian-Leg
% mu = 0.15;
% rho = 2.5; 
% beta = 0.5;
% xm=(exp(1/eps)+x^2*(y+1)^2)*(mu/rho)-1/rho*(x^2*(y+1)^2);
% ym=exp(beta/eps)*(x^2*(y+1)^2)-y*(exp(1/eps)+x^2*(y+1)^2);
% uint=[-2 2];
% vint=[-2 2];

% 5 Templator
% beta=-1/2;
% k_u=exp(-1/eps);
% K=exp(beta/eps);
% r=0.5;
% k_T=1;
% q=1;
% xm = (r-k_u*x^2-k_T*x^2*y)*(K+y);
% ym = (k_u*x^2+k_T*x^2*y)*(K+y)-q*y;
% uint=[-2 2];
% vint=[-3 2];

%6 Michalis-menten
% g1=2;
% gm1=2;
% g2=4;
% g0=2;
% k1=exp(g1/eps);
% km1=exp(gm1/eps);
% k2=exp(g2/eps);
% e0=exp(g0/eps);
% xm=-k1*e0*x+k1*x*y+km1*y;
% ym=k1*e0*x-k1*x*y-km1*y-k2*y;
% uint=[-4 4];
% vint=[-4 4];


%7 Michaelis-menten 2
% g1=4;
% gm1=2;
% g2=7;
% g0=2;
% gk=gm1-g1;
% gl=g2-g1;
% xm=-exp(g0/eps)*x+y*x+exp(gk/eps)*y;
% ym=exp(g0/eps)*x-x*y-exp(gk/eps)*y-exp(gl/eps)*y;
% uint=[-5 5];
% vint=[-5 5];

% tic
% tropicalCurves(xm,ym);
% toc

[xd,yd]=extractData(xm,ym);
% xd.signs
% xd.alphas
% yd.signs
% yd.alphas

%%%
tic
[A, L]=tropicalCurves2(xm,ym, uint, vint);
toc


%%
syms u v
epsval=0.1;
f=matlabFunction([subs(eps*xm/x, {eps,x,y}, {epsval, exp(u/epsval), exp(v/epsval)});
    subs(eps*ym/y, {eps,x,y}, {epsval, exp(u/epsval), exp(v/epsval)})]);

tspan=[0 10000000];
y0=[-0.5,0];
[t,xp]=ode23s(@(t, xp) f(xp(1), xp(2)), tspan, y0);
%[t,xp]=myRK4System(@(t, xp) f(xp(1), xp(2)), tspan, y0, 1000000);

plot(xp(:,1),xp(:,2),'LineWidth', 2)







