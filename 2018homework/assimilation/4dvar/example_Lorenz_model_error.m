clearvars;

sigma = 10;
rho = 28;
beta = 8/3;

T = 5;
dt = 1.0/300;

x = 1.508870;
y = -1.531271;
z = 25.46091;

for it = 1:round(T/dt)
    [x, y, z] = Lorenz_model(x,y,z,dt,'rk3');
end

x1 = x(1:10:end); y1 = y(1:10:end); z1 = z(1:10:end);
x2 = x(1:5:end);  y2 = y(1:5:end);  z2 = z(1:5:end);

x1_misfit = (x1(3:end)-x1(1:end-2))/(2*dt*10) - sigma* (y1(2:end-1)-x1(2:end-1));
x2_misfit = (x2(3:end)-x2(1:end-2))/(2*dt*5)  - sigma* (y2(2:end-1)-x2(2:end-1));
x3_misfit = 1/6*diff(x,3)/(dt^3)*(dt*10)^2;

figure; hold on;
plot(0:dt:T, x, 'k');
plot(dt*10:dt*10:T-dt*10, x1_misfit, 'b');
plot(dt*5:dt*5:T-dt*5,    x2_misfit, 'r');
plot(1.5*dt:dt:T-1.5*dt,  x3_misfit, 'g');
legend({'True solution','Misfit, dt=0.033','Misfit, dt=0.016','Misfit, dt=0.033,Res'},'FontSize',20);


y2_misfit = (y2(3:end)-y2(1:end-2))/(2*dt*5) - (rho*x2(2:end-1) - y2(2:end-1) - x2(2:end-1).*z2(2:end-1));
z2_misfit = (z2(3:end)-z2(1:end-2))/(2*dt*5) - (x2(2:end-1).*y2(2:end-1) - beta*z2(2:end-1));

Cqq = [x2_misfit'; y2_misfit'; z2_misfit'] * [x2_misfit y2_misfit z2_misfit] / length(x2_misfit);
disp(Cqq);

N = round(0.4/(dt*5)) + 1;
x_acor = zeros(N,1);
y_acor = zeros(N,1);
z_acor = zeros(N,1);
for i = 0:(N-1)
    x_acor(i+1) = mean(x2_misfit(1:end-i).*x2_misfit(1+i:end)) / mean(x2_misfit.^2);
    y_acor(i+1) = mean(y2_misfit(1:end-i).*y2_misfit(1+i:end)) / mean(y2_misfit.^2);
    z_acor(i+1) = mean(z2_misfit(1:end-i).*z2_misfit(1+i:end)) / mean(z2_misfit.^2);
end

figure; hold on;
plot(dt*5*(0:N-1), 0*(0:N-1), 'k');
p1 = plot(dt*5*(0:N-1), x_acor);
p2 = plot(dt*5*(0:N-1), y_acor);
p3 = plot(dt*5*(0:N-1), z_acor);

legend([p1 p2 p3], {'x component','y component','z component'},'FontSize',20);