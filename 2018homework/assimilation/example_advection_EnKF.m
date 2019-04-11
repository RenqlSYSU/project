clear all;

nx = 1000;
dx = 1;
x = 0:dx:(nx-1);

Cxx = ones(nx, nx);
for i = 1:nx
    for j = i:nx
        dij = min(x(j)-x(i), x(i)+nx*dx-x(j));
        Cxx(i,j) = Cxx(i,j) * exp(-dij^2/20^2);
        Cxx(j,i) = Cxx(i,j);
    end 
end
figure; imagesc(Cxx); colorbar; pause;

f0 = mvnrnd(zeros(1,nx), Cxx)';
figure; plot(f0); colorbar; pause;

c = 1;
dt = 1;
nt = 1000;

ft = f0; % True value

G = eye(nx, nx) * (1-c*dt/dx) + diag(ones(nx-1,1), -1) * (c*dt/dx);
G(1,nx) = c*dt/dx;

nset = 50;
faset = repmat(f0, [1 nset]) + mvnrnd(zeros(nset,nx), Cxx)';
% Cqq = zeros(nx, nx);
Cqq = Cxx/100;

Coo = 0.1 * eye(4);
M = zeros(4, nx);
M(1,100) = 1;
M(2,350) = 1;
M(3,600) = 1;
M(4,850) = 1;

for it = 1:nt
    ft = G * ft;
    
    d = M * ft + mvnrnd(zeros(1,4), Coo)';
    fset = G * faset + mvnrnd(zeros(nset,nx), Cqq)';
    fmean = mean(fset, 2);
    Cff = 1/(nset-1) * (fset - repmat(fmean,[1 nset])) * (fset - repmat(fmean,[1 nset]))';
    K = Cff * M' / (M * Cff * M' + Coo);
    Caa = (eye(nx) - K * M) * Cff;
    dset = repmat(d,[1 nset]) + mvnrnd(zeros(nset,4), Coo)';
    faset = fset + K * (dset - M * fset);
    
    clf;
    subplot(2,1,1);
    hold on;
    plot(ft, 'r');
    plot(mean(faset,2), 'b');
    scatter([100 350 600 850], d, 50, 'filled','dk');
    subplot(2,1,2);
    plot(diag(Caa));
    title(['t = ' num2str(it)]);
    pause(0.1);
end
