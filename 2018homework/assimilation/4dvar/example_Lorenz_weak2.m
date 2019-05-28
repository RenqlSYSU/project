
if 1
    
    T = 20;
    dt = 1.0/60;
    t = (0:dt:T)';
    nt = length(t);
    
    xt = 1.508870;
    yt = -1.531271;
    zt = 25.46091;
    for it = 1:round(T/(dt/10))
        [xt, yt, zt] = Lorenz_model(xt,yt,zt,dt/10,'rk3');
    end
    xt = xt(1:10:end);
    yt = yt(1:10:end);
    zt = zt(1:10:end);
    
    var_obs = 2;
    x0 = xt(1) + random('Normal',0,var_obs, 1,1);
    y0 = yt(1) + random('Normal',0,var_obs, 1,1);
    z0 = zt(1) + random('Normal',0,var_obs, 1,1);
    
    dobs = round(0.25/dt);
    lobs = 1+dobs:dobs:nt;
    nobs = numel(lobs);
    obs = [xt(lobs) yt(lobs) zt(lobs)] + random('Normal', 0,var_obs, nobs,3);
    
    Wqq = [0.1491 0.1505 0.0007; 0.1505 0.9048 0.0014; 0.0007 0.0014 0.9180]^(-1);
    Waa = 1/var_obs*eye(3);
    Wee = 1/var_obs*eye(3);
    
    x = zeros(nt,1);
    y = zeros(nt,1);
    z = ones(nt,1)*23;
end

figure;
% subplot(2,2,1);  plot(t, xt, 'k'); hold on; title('x');
% subplot(2,2,2);  plot(t, yt, 'k'); hold on; title('y');
% subplot(2,2,3);  plot(t, zt, 'k'); hold on; title('z');

r = 0.002;
delta = 0.002;

niter = 200;
J = zeros(niter,1);
J1 = zeros(niter,1);
J2 = zeros(niter,1);
J3 = zeros(niter,1);
J4 = zeros(niter,1);
[J(1), J1(1), J2(1), J3(1), J4(1)] = Jweak2(dt, [x; y; z], lobs, obs, x0, y0, z0, Wqq, Waa, Wee);

options = optimoptions(@simulannealbnd,'MaxIterations',2000);

for iter = 2:niter
    
    xyz = [x; y; z];
    
    fun = @(var) Jweak2(dt, var, lobs, obs, x0, y0, z0, Wqq, Waa, Wee);
    xyz = simulannealbnd(fun,xyz,[],[],options);
    x = xyz(1:nt);
    y = xyz(nt+1:2*nt);
    z = xyz(2*nt+1:3*nt);    
    
    [J(iter),J1(iter),J2(iter),J3(iter),J4(iter)] ...
        = Jweak2(dt, xyz, lobs, obs, x0, y0, z0, Wqq, Waa, Wee);

    subplot(2,2,1); cla; plot(t, xt, 'k'); hold on; title('x'); plot(t, x, 'r'); scatter(t(lobs),obs(:,1),'o','filled');
    subplot(2,2,2); cla; plot(t, yt, 'k'); hold on; title('y'); plot(t, y, 'b'); scatter(t(lobs),obs(:,2),'o','filled');
    subplot(2,2,3); cla; plot(t, zt, 'k'); hold on; title('z'); plot(t, z, 'g'); scatter(t(lobs),obs(:,3),'o','filled');
    subplot(2,2,4); cla; hold on;
    plot(1:iter, J(1:iter));   plot(1:iter, J1(1:iter)); 
    plot(1:iter, J3(1:iter));  plot(1:iter, J4(1:iter)); 
    legend({'Total', 'Dyn', 'Obs', 'Smooth'});
    xlim([1 niter]); ylim([0 J(1)*1.2]);
    pause(0.1);
end