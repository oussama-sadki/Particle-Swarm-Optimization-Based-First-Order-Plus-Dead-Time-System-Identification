clear all
close all
clc
data = load('step_test.csv');
Q1 = data(:,1); % Input
T1 = data(:,2); % Output
% Initial condition
time = 0:1:length(Q1)-1;
ym = zeros(length(Q1),1);
ym0 =  T1(1);
ym(1) = ym0;
um = Q1;
u = zeros(length(Q1),1);
u0 = um(1);
u(1) = u0;
n = length(ym);
yp = T1;
%% Optimization
% PSO Constants Parameters
N = 1000; % Particle
iter = 50;
var = 3;
w = 0.9;
c = 1;
% Search space
a = 0.0;
b = 150;
% Optimization step
c_cf = 0;
% Initialization
for m=1:N
    for j=1:var
        v(m,j) = 0;
        x(m,j) = a+rand*(b-a);
        xp(m,j) = x(m,j);
    end
    kp = x(m,1);
    theta = round(x(m,2));
    tau = x(m,3);
    % Simulate the model
    for i=1:n-1
          if (time(i)-theta) <= 0
              u(i) = u0;
          else
              u(i) = um(time(i)-theta);
          end
        dydt = (-(ym(i)-ym0)+kp*(u(i)-u0))/tau;
        ym(i+1) = ym(i)+dydt*(time(i+1)-time(i));
    end 
    % ISE Objective function
    ff1 = 0;
    sim1 = size(ym);
    for m1=1:sim1
        ff1 = ff1+(ym(m1)-yp(m1))^2;
    end
    ISE(m) = ff1;
end
% Find the best solution
[Best_performance,location] = min(ISE);
fg = Best_performance;
xg(1) = x(location,1);
xg(2) = x(location,2);
xg(3) = x(location,3);
k = 0;
while(k<iter)
    k = k+1;
    for m=1:N
        for j=1:var
            v(m,j) = w*v(m,j)+c*rand*(xp(m,j)-x(m,j))+c*rand*(xg(j)-x(m,j));
            x(m,j) = x(m,j)+v(m,j);
        end
        % Check bound
        for j=1:var
            if x(m,j)<a
                x(m,j) = a;
            end
            if x(m,j)>b
                x(m,j) = b;
            end
        end
        % Model Parameters
        kp = x(m,1);
        theta = round(x(m,2));
        tau = x(m,3);
        % Simultate model
        for i=1:n-1
             if (time(i)-theta) <= 0
                 u(i) = u0;
             else
                 u(i) = um(time(i)-theta);
             end
            dydt = (-(ym(i)-ym0)+kp*(u(i)-u0))/tau;
            ym(i+1) = ym(i)+dydt*(time(i+1)-time(i));
        end
        % ISE Objective function
        ff1 = 0;
        sim1 = length(ym);
        for m1=1:sim1
            ff1 = ff1+(ym(m1)-yp(m1))^2;
        end
        ISEp(m) = ff1;
        % Compare local
        if ISEp(m)<ISE(m)
            ISE(m) = ISEp(m);
            xp(m,1) = x(m,1);
            xp(m,2) = x(m,2);
            xp(m,3) = x(m,3);
        end
    end
    [B_fg,location] = min(ISE);
    % Compare global
    if B_fg<fg
    fg = B_fg;
    xg(1) = xp(location,1);
    xg(2) = xp(location,2);
    xg(3) = xp(location,3);
    end
    c_cf = c_cf+1;
    best_cf_ac(c_cf) = fg;
    kp = xg(1);
    theta = round(xg(2));
    tau = xg(3);
    for i=1:n-1
             if (time(i)-theta) <= 0
                 u(i) = u0;
             else
                 u(i) = um(time(i)-theta);
             end
            dydt = (-(ym(i)-ym0)+kp*(u(i)-u0))/tau;
            ym(i+1) = ym(i)+dydt*(time(i+1)-time(i));
    end
    figure(2)
    plot(time,yp,'b',time,ym,'r--','LineWidth',1.5);grid minor
    xlabel('time(s)');
    ylabel('Temperature °C')
    legend('Data','Model');
    title('Data and System Model fitting')
    pause(0.1)
    k
end
figure(3)
plot(time,u,'k','LineWidth',1);grid minor
xlabel('time(s)');
ylabel('Heater Level');
title('Input');
figure(4);
t_cf = 1:c_cf;
plot(t_cf,best_cf_ac,'r--','LineWidth',2)
xlabel('Iteration')
ylabel('Cost Function')
grid minor
legend('ISE')
title('ISE with each iteration')
min_ISE = fg
kp = xg(1)
theta = xg(2)
tau = xg(3)
params = [kp,theta,tau];
csvwrite('params.csv',params);  