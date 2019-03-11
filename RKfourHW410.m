%Constants:
%%
y0= [0 10^(-5)];
y_0= [10^(-5) 0 ];
tspan= [-0.6 0];

% Using the fourth-order explicit Runge-Kutta solver RK4 to solve the Schrodinger equation on the domain [−x0, x0], where we take x0 = 0.6 nm and set E = 1 eV for now. Initial conditions Psi(-x0) = 0 and Psi'(-x0) = epsilon = 1E-5, and psi(0) is non-zero on entire domain.


%[y,t]=RK4(@(x,y) f(x,y,1), tspan, y0);
%plot(t,y);

[y,t]=RKfa(@(x,y) f(x,y,5.83), tspan, y0);
figure(1)
plot(t,y)

figure(2)
hold on;
title("Schrodinger Equation solution on [-0.6,0.6] domain, E = 22.6")
xlabel("time")
ylabel("Psi")

[y2,t2]=RKfa(@(x,y) f(x,y,22.6), tspan, y0);
plot(t2,y2)
hold off;


% In a single figure, plotting the six corresponding eigenstates, as well as the potential V (x), on the domain [−x0, x0]. Align the zero baseline of each wavefunction with its energy eigenvalue, and use the MATLAB function trapz to normalize appropriately

[a, b] = RKfa(@(x,y) f(x,y,5.8896), tspan, y0);
b1=[ flipud(b(:,1)) ; -b(:,1)]; 
a1=[ flipud(a(:,1)) ; -a(:,1)]; 
a1t= (a1/sqrt(trapz(a1.^2)));



[c,d]=RKfa(@(x,y) f(x,y, 5.8533), tspan, y_0);
d2=[ flipud(d(:,1)) ; -d(:,1)]; 
c2=[ flipud(c(:,1)) ; c(:,1)]; 
c2t= (c2/sqrt(trapz(c2.^2)));


[e,ff]=RKfa(@(x,y) f(x,y, 16.6691 ), tspan, y0);
ff=[ flipud(ff(:,1)) ; -ff(:,1)]; 
e3=[ flipud(e(:,1)) ; e(:,1)]; 

e3t= (e3/sqrt(trapz(e3.^2)));


[g,h]=RKfa(@(x,y) f(x,y,15.3591), tspan, y_0);
h4=[ flipud(h(:,1)) ; -h(:,1)]; 
g4=[ flipud(g(:,1)) ; g(:,1)]; 
g4t= (y3/sqrt(trapz(g4.^2)));

[i,j]=RKfa(@(x,y) f(x,y, 27.6577), tspan, y0);
j5=[ flipud(j(:,1)) ; -j(:,1)]; 
i5=[ flipud(i(:,1)) ; i(:,1)]; 
i5t= (i5/sqrt(trapz(i5.^2)));


[k,l]=RKfa(@(x,y) f(x,y, 22.5648), tspan, y_0);
l6=[ flipud(l(:,1)) ; -l(:,1)]; 
k6=[ flipud(k(:,1)) ; k(:,1)]; 
k6t= (k6/sqrt(trapz(k6.^2)));

disp(trapz(k6t.^2))

figure(3)
plot(b1,a1t + 5.8896);
hold on;
plot(d2,c2t + 5.8533);
hold on;
plot(e3,ff3t + 16.6691);
hold on;
plot(h4,g4t+15.3591);
hold on;
plot(j5,i5t+27.6577);
hold on;
plot(l6,k6t+22.5648);

ylabel('Psi');
xlabel('x');
title('Psi(x) vs x (nm), Energy Eigenvectors')
%xlim([-0.6  1]);


figure(4)
plot(b1,a1t + 5.8896);
xlabel('x')
ylabel('Psi')
title('Psi(x) vs x (nm) for E = 5.889')



function PsiPri= f(x, y,E)
    %E= 5.8896;
    alpha= 500;
    beta=3500;
    hbar= 0.076199682;

    PsiPri= [y(2) (-2/hbar)*(E-(-alpha*x^2+beta*x^4+(alpha^2/(4*beta))))*y(1)];
end


% Explicit Fourth-Order Runge-Kutta with constant step size

function [y,t] = RKfa(fun,tspan,y0)

% f: function evaluating the time derivatives; y'(t) = f(t,y)
% Example f=@(t,y) y, uses function handle to solve y'=y;
% Note that f must accept TWO arguments: t and y(t). For systems of N first
% order ODEs, f and y should both output row vectors with N components.
% y0: row vector of initial conditions
% tspan = [t0 tf] is the time interval over which we solve the ODE

%% Parameters

h = 1e-2;
t0 = tspan(1); tf = tspan(2);

%% Initilization

iter = round((tf-t0)/h);        % number of time steps
y = zeros(iter+1,length(y0));   % preallocating
y(1,:) = y0;                    % y(t0) = y0

for i=1:iter
    k1 = feval(fun, t0 + (i-1)*h,       y(i,:)            );
    k2 = feval(fun, t0 + (i-1)*h + h/2, y(i,:) + (h/2)*k1 );
    k3 = feval(fun, t0 + (i-1)*h + h/2, y(i,:) + (h/2)*k2 );
    k4 = feval(fun, t0 + (i-1)*h + h,   y(i,:) + h*k3     );
    y(i+1,:) = y(i,:) + (h/6)*(k1+2*k2+2*k3+k4);
end

t = t0 + h*(0:iter)';

end


% My own RK4 Method:

function [t2, tt2, F2] = rk4()
close all

%Constants:
m = 1;
v = v;
%A = A;
%v = 1;
g = 1;

% Functions

%fT = @(F, theta, t) -1/v * (m*fF(F, theta, t) + g*sin(theta));
%fF = @(F, theta, t) (-1/m) * (v*F + g*sin(theta));

% initial conditions
tt(1) = 1;
F(1) = 0;
t(1) = 0;


% step size

% Still need to optimize dt
dt = 0.001;
tfinal=3000;
N = ceil(tfinal/dt);

% Update loop
for n = 1:N
    % Update time
    t(n+1) = t(n) + dt;
    % Update F and T
    k1F = f(F(n)          , tt(n)           ,        t(n));
    k2F, k2T = f(F(n)+ k1F*dt/2, tt(n) + k1F*dt/2, t(n) + dt/2);
    k3F, k3T = f(F(n)+ k2F*dt/2, tt(n) + k2T*dt/2, t(n) + dt/2);
    k4F, k4T  = f(F(n)+ k3F*dt/2, tt(n) + k3T*dt/2, t(n) + dt);
    %phi = (1/6)*(k1+2*k2 + 2*k3 + k4);
    %y4(n+1)= y4(n) + phi*dt
    %n=n+1
    F(n+1) = F(n) + dt/6 * (k1F + 2*k2F + 2*k3F + k4F);
    tt(n+1) = tt(n) + dt/6 * (k1T + 2*k2T + 2*k3T + k4T);
end

%figure(1);
%plot(t,theta,'b-','Linewidth',2);
%figure(2);
%plot(theta, F, 'b-','Linewidth',2);
%plot(F, tt)
t2 = t;
tt2 = tt;
F2 = F;

function ydot= f(x, y,E)
    %E= 5.8896;
    alpha= 500;
    beta=3500;
    hbar= 0.076199682;

    ydot= [y(2) (-2/hbar)*(E-(-alpha*x^2+beta*x^4+(alpha^2/(4*beta))))*y(1)];
end


end
