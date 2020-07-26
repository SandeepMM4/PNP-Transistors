function term_project_MM655()




close all;

meu=0.460e-3; % mobility of holes in m2/Vs 
k=1.38e-23;   % boltzman constant in m2kg/s2k
temp=300;     % temperature in kelvin
e=1.67e-19;   % electron charge in coulomb
D= meu*k*temp*10^4/e; % diffusivity of holes in the base of the PNP transistor in cm2/s

% Grid setup
dx=1; % micronmeter               
dt=1;  % ms
L=50;  % micrometer
T=50;  % ms

%time array
t= 0:dt:T;
%x array
x=0:dx:L;

% conc of holes along the length of the base
u(1:length(x))=0.0;
% spatial derivative
d2udx2(1:length(x))=0.0;
dudx(1:length(x))=0.0;
% time derivative
dudt(1:length(x))=0.0;
% analytical concentration of holes along the length of the base
u_analytical(1:length(x))=0.0;
% error belween actual and analytical concentration of holes along the
% length of the base
error(1:length(x))=0.0;

%boundary condition
u(1)=0;
u(L)=0;

%constants
eta=10; 
K=1;
         
% Loop over x and t and update the conc. using the diffusion equation
for it=1:length(t)
    for ix=2:(length(x)-1)
        % spatial derivative d2udx2
        d2udx2(ix)=(u(ix+1)-2*u(ix)+u(ix-1))/(dx^2);
    end
        %spatial derivative dudx
    for ix=2:(length(x)-1)
        dudx(ix)= (u(ix+1)-u(ix-1))/2*dx;
        dudt(ix)=D*(d2udx2(ix))-(D*eta/L)*(dudx(ix));
    end
    % updating concentration u(x,t)
    for ix=2:length(x)-1
        % initial condition for time=0 sec
        if it==1
            u(ix) = (K*L/D)*(1 - exp(-eta*(1 - ix/L)))/eta;
        elseif it>1 && it<=50
            % concentration of holes other than initial condition
            u(ix)=u(ix)+dudt(ix)*dt;
        end
    
    end
end


figure(1); 
plot(x,u,'linewidth',5)
title('Numerical result for holes concentration');
xlabel('x(micrometer)');
ylabel('u(x,t) /cm3');
xlim auto;
ylim auto;


% for analytical solution
%considering the given data for calculating the analytical concentration
ul=15;% equilibrium concentration of holes, u1=n^2/N , N is donor conc. and n is intrinsic carrier concentration
q=1.67e-19;  % charges in coulomb
v=3e-7;      % emitter-base voltage (v)

%Calculate and plot the analytical conc. profile at time t = 50 s
for ix=1:length(x)
    u_analytical(ix)=ul*((exp(q*v/k*temp)-1)*(1-ix/L));
end

figure(2);
plot(x,u_analytical,'linewidth',3);
title('Analytical Solution of Concentration of Holes in Base of PNP Transistor');
xlabel('x (micrometer)');
ylabel('u(x,t) /cm3');
xlim auto;
ylim auto;

% Now calculate the difference between the numerical and analytical concentration
for ix=1:length(x)
    error(ix)=abs(u(ix)-u_analytical(ix));
end

% Plot the error profile
figure(3);
plot(x,error,'linewidth',3)
title('Compare between Actual and Analytical');
xlabel('x');
ylabel('error');
xlim auto;
ylim auto;
hold off;

%Effect of dt with corresponding changing dx

% Grid setup
dx1=0.9;% micronmeter  
dx2=0.4;% micronmeter  

dt1=0.1;% ms
dt2=0.05;% ms

L=50;  % micrometer
T=50;  % ms

%time array
t1= 0:dt1:T;
t2= 0:dt2:T;

%x array
x1=0:dx1:L;
x2=0:dx2:L;

% conc of holes along the length of the base
u1(1:length(x1))=0.0;
u2(1:length(x2))=0.0;

% spatial derivative
d2udx2(1:length(x1))=0.0;
d2udx2(1:length(x2))=0.0;

dudx(1:length(x1))=0.0;
dudx(1:length(x2))=0.0;

% time derivative
dudt(1:length(x1))=0.0;
dudt(1:length(x2))=0.0;

%boundary condition
u1(1)=0;
u1(L)=0;

u2(1)=0;
u2(L)=0;

% Loop over x and t and update the conc. using the diffusion equation
% case:1 for dx1 and dt1
for it=1:length(t1)
    for ix=2:(length(x1)-1)
        % spatial derivative d2udx2
        d2udx2(ix)=(u1(ix+1)-2*u1(ix)+u1(ix-1))/(dx1^2);
    end
        %spatial derivative dudx
    for ix=2:(length(x1)-1)
        dudx(ix)= (u1(ix+1)-u1(ix-1))/2*dx1;
        dudt(ix)=D*(d2udx2(ix))-(D*eta/L)*(dudx(ix));
    end
    % updating concentration u1(x,t)
    for ix=2:length(x1)-1
        % initial condition for time=0 sec
        if it==1
            u1(ix) = (K*L/D)*(1 - exp(-eta*(1 - ix/L)))/eta;
        elseif it>1 && it<=50
            % concentration of holes other than initial condition
            u1(ix)=u1(ix)+dudt(ix)*dt1;
        end
    
    end
end

% case:2 for dx2 and dt2

for it=1:length(t2)
    for ix=2:(length(x2)-1)
        % spatial derivative d2udx2
        d2udx2(ix)=(u2(ix+1)-2*u2(ix)+u2(ix-1))/(dx2^2);
    end
        %spatial derivative dudx
    for ix=2:(length(x2)-1)
        dudx(ix)= (u2(ix+1)-u2(ix-1))/2*dx2;
        dudt(ix)=D*(d2udx2(ix))-(D*eta/L)*(dudx(ix));
    end
    % updating concentration u2(x,t)
    for ix=2:length(x2)-1
        % initial condition for time=0 sec
        if it==1
            u2(ix) = (K*L/D)*(1 - exp(-eta*(1 - ix/L)))/eta;
        elseif it>1 && it<=50
            % concentration of holes other than initial condition
            u2(ix)=u2(ix)+dudt(ix)*dt2;
        end
    
    end
end

figure(4);
a1=subplot(3,1,1);
a2=subplot(3,1,2);
a3=subplot(3,1,3);

plot(a1,x,u,'linewidth',5);
title(a1,'Numerical Solution of Concentration of Holes with dt=1 ms,dx=1 micrometer')
xlabel(a1,'x(micrometer)');
ylabel(a1,'u(x,t) /cm3');
xlim auto;
ylim auto;
hold on;
plot(a2,x1,u1,'linewidth',5);
title(a2,'Numerical Solution of Concentration of Holes with dt=0.1 ms,dx=0.9 micrometer')
xlabel(a2,'x1(micrometer)');
ylabel(a2,'u1(x,t) /cm3');
xlim auto;
ylim auto;
plot(a3,x2,u2,'linewidth',5);
title(a3,'Numerical Solution of Concentration of Holes with dt=0.05 ms,dx=0.4 micrometer')
xlabel(a3,'x2(micrometer)');
ylabel(a3,'u2(x,t) /cm3');
xlim auto;
ylim auto;


