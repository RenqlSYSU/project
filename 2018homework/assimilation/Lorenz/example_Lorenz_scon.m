clearvars;

T = 2; % 2,4,8
Waa = 2;
Wee = 20*eye(3);

xset = -25:0.01:25;

xt = 1.5;
yt = 0;
zt = 25;
dt = 0.01;
for t = 1:ceil(T/dt)
    [xt, yt, zt] = Lorenz_model(xt, yt, zt, dt, 'rk3');
end

dobs = 0.5;
tobs = dobs:dobs:T;
nobs = length(tobs);
d = zeros(3,nobs);
for iobs = 1:nobs
    it = min(round(tobs(iobs)/dt),ceil(T/dt));
    d(:,iobs) = [xt(it); yt(it); zt(it)] + mvnrnd(zeros(1,3), Wee^(-1))';
end

Jset = zeros(numel(xset),1);
J1 = Jset;
J2 = Jset;
for ix = 1:numel(xset)
    x = xset(ix);
    y = yt(1);
    z = zt(1);
    for t = 1:ceil(T/dt)
        [x, y, z] = Lorenz_model(x, y, z, dt, 'rk3');
    end
    
    J1(ix) = (x(1)-xt(1))*Waa*(x(1)-xt(1));
    J2(ix) = 0;
    for iobs = 1:nobs
        it = min(round(tobs(iobs)/dt),ceil(T/dt));
        J2(ix) = J2(ix) + (d(:,iobs)-[x(it); y(it); z(it)])' * Wee * (d(:,iobs)-[x(it); y(it); z(it)]);
    end
    Jset(ix) = J1(ix) + J2(ix);
end

figure;
subplot(3,1,1); hold on;
plot(xset,J1,'k','Linewidth',1);
plot([1.5 1.5],ylim, 'r');
title('J1');
[~, xloc] = min(J1);
disp(['Mininum of J1 : ' num2str(xset(xloc))]);

subplot(3,1,2); hold on;
plot(xset,J2,'k','Linewidth',1);
plot([1.5 1.5],ylim, 'r');
title('J2');
[~, xloc] = min(J2);
disp(['Mininum of J2 : ' num2str(xset(xloc))]);

subplot(3,1,3); hold on;
plot(xset,Jset,'k','Linewidth',1);
plot([1.5 1.5],ylim, 'r');
title('J');
[~, xloc] = min(Jset);
disp(['Mininum of J : ' num2str(xset(xloc))]);