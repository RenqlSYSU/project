%Kalman Filtering, take upwind scheme for example
%assume nx=1000 grids and nt=1000 iterations,
clear all;

nx = 1000;  %Number of grids
dx = 1;     %grid spacing
x = 0:dx:(nx-1); %create an array of 1*nx

%-------The analysis error at initial time is given by people--------
Cxx = ones(nx, nx); %Create array of all one
for i = 1:nx
    for j = i:nx
        dij = min(x(j)-x(i), x(i)+nx*dx-x(j));
        Cxx(i,j) = Cxx(i,j) * exp(-dij^2/20^2);
        Cxx(j,i) = Cxx(i,j);
    end 
end
figure; imagesc(Cxx); colorbar;
title('C_{aa}^0');

%-------create initial value using initial analysis error-------------------
f0 = mvnrnd(zeros(1,nx), Cxx); %returns a matrix R of random vectors chosen from the 
%same multivariate normal distribution, with mean vector first number and covariance matrix second number. 
%zeros means Create array of all zero
figure; plot(f0); colorbar; 
title('Initial value');
pause;

c = 1;
dt = 1;    %time step
nt = 2000; %number of iteration

ft = f0; % True value

%---create predictor array-----------------
G = eye(nx, nx) * (1-c*dt/dx) + diag(ones(nx-1,1), -1) * (c*dt/dx); 
%Create diagonal matrix (the first number must be vector) or get diagonal elements of matrix (the first number is matrix)
G(1,nx) = c*dt/dx;

fa0 = f0 + mvnrnd(zeros(1,nx), Cxx)';
fa = fa0; %analysis value
Caa = Cxx; %analytical error
Cqq = zeros(nx, nx); %predictor error
% Cqq = Cxx/10;

%-----obseration-------------------------------
Coo = 0.1 * eye(4);%The observation error
M = zeros(4, nx);  %The observation operator
M(1,100) = 1;
M(2,350) = 1;
M(3,600) = 1;
M(4,850) = 1;

dt_obs = 1;
for it = 1:nt
    ft = G * ft; %true value
    
    ff = G * fa; %predicted value
    
    if mod(it,dt_obs) == 0 %Remainder after division
        d = M * ft + mvnrnd(zeros(1,4), Coo)';
        Cff = G * Caa * G' + Cqq; %prediction error
        K = Cff * M' / (M * Cff * M' + Coo);
        Caa = (eye(nx) - K * M) * Cff;
        fa = ff + K * (d - M * ff);
    else
        fa = ff;
    end
    
    subplot(2,1,1); cla;%subplot(m,n,p) divides the current figure into an m-by-n grid and creates axes in the position specified by p.
    hold on;
    plot(ft, 'r');
    plot(fa, 'b');
    if mod(it,dt_obs) == 0
        scatter([100 350 600 850], d, 50, 'filled','dk');
    end
    title(['t = ' num2str(it)]);
    
    if mod(it,dt_obs) == 0
        subplot(2,1,2); cla;  %Clear axes
        plot(diag(Caa));
        pause(0.1);
    else
        pause(0.1);
    end
end
