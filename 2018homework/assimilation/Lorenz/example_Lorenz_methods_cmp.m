clear all;

dt = 0.0001;
T = 2;

x1 = 0.01;
y1 = 0;
z1 = 0;

x2 = 0.01;
y2 = 0;
z2 = 0;

x3 = 0.01;
y3 = 0;
z3 = 0;

for t = 1:ceil(T/dt)
    [x1, y1, z1] = Lorenz_model(x1,y1,z1,dt,'exp');
    [x2, y2, z2] = Lorenz_model(x2,y2,z2,dt,'lf');
    [x3, y3, z3] = Lorenz_model(x3,y3,z3,dt,'rk3');
    
    if mod(t,10) == 0
        subplot(2,2,1);  cla; hold on;
        plot(x1, y1, 'b'); plot(x2, y2, 'r'); plot(x3, y3, 'g');
        scatter(x1(end),y1(end),50,'ob','filled');
        scatter(x2(end),y2(end),50,'or','filled');
        scatter(x3(end),y3(end),50,'og','filled');
        title('x-y'); xlabel('x'); ylabel('y');
        
        subplot(2,2,2);  cla; hold on;
        plot(y1, z1, 'b'); plot(y2, z2, 'r'); plot(y3, z3, 'g');
        scatter(y1(end),z1(end),50,'ob','filled');
        scatter(y2(end),z2(end),50,'or','filled');
        scatter(y3(end),z3(end),50,'og','filled');
        title('y-z'); xlabel('y'); ylabel('z');
        
        subplot(2,2,3);  cla; hold on;
        plot(z1, x1, 'b'); plot(z2, x2, 'r'); plot(z3, x3, 'g');
        scatter(z1(end),x1(end),50,'ob','filled');
        scatter(z2(end),x2(end),50,'or','filled');
        scatter(z3(end),x3(end),50,'og','filled');
        title('z-x'); xlabel('z'); ylabel('x');
        
        pause(0.02);
    end
end