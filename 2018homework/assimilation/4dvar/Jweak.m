function [Jtot, Jdyn, Jinit, Jobs, Jsmth] = Jweak (dt, x, y, z, lobs, obs, x0, y0, z0, ...
    Wqq, Waa, Wee)

sigma = 10;
rho = 28;
beta = 8/3;

gamma = 0.0001;

q = (x(3:end)-x(1:end-2))/(2*dt) - sigma*(y(2:end-1)-x(2:end-1));
q = [q (y(3:end)-y(1:end-2))/(2*dt)-(rho*x(2:end-1)-y(2:end-1)-x(2:end-1).*z(2:end-1))];
q = [q (z(3:end)-z(1:end-2))/(2*dt)-(x(2:end-1).*y(2:end-1)-beta*z(2:end-1))];
Jdyn = dt * sum(q(:,1).^2*Wqq(1,1) + 2*q(:,1).*q(:,2)*Wqq(1,2) + 2*q(:,1).*q(:,3)*Wqq(1,3) ...
    + q(:,2).^2*Wqq(2,2) + 2*q(:,2).*q(:,3)*Wqq(2,3) + q(:,3).^2*Wqq(3,3));

Jinit = [x(1)-x0 y(1)-y0 z(1)-z0] * Waa * [x(1)-x0 y(1)-y0 z(1)-z0]';

Jobs = (x(lobs)-obs(:,1))' * Wee(1,1) * (x(lobs)-obs(:,1)) ...
    +  (y(lobs)-obs(:,2))' * Wee(2,2) * (y(lobs)-obs(:,2)) ...
    +  (z(lobs)-obs(:,3))' * Wee(3,3) * (z(lobs)-obs(:,3));

eta = (x(3:end)+x(1:end-2)-2*x(2:end-1))/(dt^2);
eta = [eta (y(3:end)+y(1:end-2)-2*y(2:end-1))/(dt^2)];
eta = [eta (z(3:end)+z(1:end-2)-2*z(2:end-1))/(dt^2)];
Jsmth = dt * gamma * sum(eta(:,1).^2 + eta(:,2).^2 + eta(:,3).^2);

Jtot = Jdyn + Jinit + Jobs + Jsmth;