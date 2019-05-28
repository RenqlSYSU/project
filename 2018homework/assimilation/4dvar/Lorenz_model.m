function [x,y,z] = Lorenz_model(xp,yp,zp,dt,tp)

sigma = 10;
rho = 28;
beta = 8/3;

if strcmp(tp,'exp')
    x = xp(end) + dt * sigma * (yp(end) - xp(end));
    y = yp(end) + dt * (rho * xp(end) - yp(end) - xp(end) * zp(end));
    z = zp(end) + dt * (xp(end) * yp(end) - beta * zp(end));
elseif strcmp(tp,'lf')
    if (numel(xp) == 1)
        x = xp + dt * sigma * (yp - xp);
        y = yp + dt * (rho * xp - yp - xp * zp);
        z = zp + dt * (xp * yp - beta * zp);
    else
        x = xp(end-1) + 2*dt * sigma * (yp(end) - xp(end));
        y = yp(end-1) + 2*dt * (rho * xp(end) - yp(end) - xp(end) * zp(end));
        z = zp(end-1) + 2*dt * (xp(end) * yp(end) - beta * zp(end));
    end
elseif strcmp(tp,'rk3')
    k1 = [sigma*(yp(end)-xp(end)) rho*xp(end)-yp(end)-xp(end)*zp(end) xp(end)*yp(end)-beta*zp(end)];
    xx = [xp(end) yp(end) zp(end)] + 0.5*dt * k1;
    k2 = [sigma*(xx(2)-xx(1)) rho*xx(1)-xx(2)-xx(1)*xx(3) xx(1)*xx(2)-beta*xx(3)];
    xx = [xp(end) yp(end) zp(end)] - dt * k1 + 2*dt * k2;
    k3 = [sigma*(xx(2)-xx(1)) rho*xx(1)-xx(2)-xx(1)*xx(3) xx(1)*xx(2)-beta*xx(3)];
    xx = [xp(end) yp(end) zp(end)] + dt/6 * (k1 + 4*k2 + k3);
    x = xx(1);
    y = xx(2);
    z = xx(3);
end

x = [xp; x];
y = [yp; y];
z = [zp; z];