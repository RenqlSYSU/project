clear all;

dt = 0.01;
T = 100;

x1 = -0.01;
y1 = 0;
z1 = 0;

x2 = 0.01;
y2 = 0;
z2 = 0;

l2e = [];
for t = 1:ceil(T/dt)
    [x1, y1, z1] = Lorenz_model(x1,y1,z1,dt,'rk3');
    for tau = 1:ceil(dt/0.0001)
        [x2, y2, z2] = Lorenz_model(x2,y2,z2,0.0001,'rk3');
    end
    
    if mod(t,10) == 0
        subplot(2,2,1);  cla; hold on;
        plot(x1, y1, 'b'); plot(x2, y2, 'r');
        scatter(x1(end),y1(end),50,'ob','filled');
        scatter(x2(end),y2(end),50,'or','filled');
        title('x-y'); xlabel('x'); ylabel('y');
        
        subplot(2,2,2);  cla; hold on;
        plot(y1, z1, 'b'); plot(y2, z2, 'r');
        scatter(y1(end),z1(end),50,'ob','filled');
        scatter(y2(end),z2(end),50,'or','filled');
        title('y-z'); xlabel('y'); ylabel('z');
        
        subplot(2,2,3);  cla; hold on;
        plot(z1, x1, 'b'); plot(z2, x2, 'r');
        scatter(z1(end),x1(end),50,'ob','filled');
        scatter(z2(end),x2(end),50,'or','filled');
        title('z-x'); xlabel('z'); ylabel('x');
        
        l2e = [l2e sqrt((x1(end)-x2(end))^2+(y1(end)-y2(end))^2+(z1(end)-z2(end))^2)];
        subplot(2,2,4); cla;
        plot(l2e);
        title([num2str(t) '/' num2str(ceil(T/dt))]);
        
        pause(0.01);
    end
end