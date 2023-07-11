clc
clear all
close all

%% Defining the problem domain
n = 51; % Number of points
dom_length = 1; % its the same for y direction 
h = dom_length/(n-1); %space steps = 0.02
x = 0:h:dom_length; %X domain cordinates
y = 0:h:dom_length; %Y domain cordinates
dt = 0.0001; %Time steps
Re = 100; %Reynolds number
delta = 4.5; %Artificial compressibility

%% Initializing the variables
%Final collocated variables
u_final(n,n)=0;
v_final(n,n)=0;
p_final(n,n)=1;
u_final(1,:) = 1;

%Staggered variables
u(n+1,n)=0;
v(n,n+1)=0;
p(n+1,n+1)=1;
u(1,:)=2;

u_new(n+1,n)=0;
v_new(n,n+1)=0;
p_new(n+1,n+1)=1;
u_new(1,:)=2;

%% Solving the governing equations
% error = 1;
% iterations = 0;
% error_req = 1e-6; %final allowed error 
% figure(1); %for error monitoring

% while error > error_req
    % x-momentum eq.
    for k=1:100000
    for i = 2:n
        for j = 2:n - 1
            pressure = -(p(i,j+1) - p(i,j))/h;
            diffusion = (1/Re)*((u(i+1,j) - 2*u(i,j) + u(i-1,j))/(h*h) + (u(i,j+1) - 2*u(i,j) + u(i,j-1))/(h*h));
            advection_x = ((0.5*(u(i,j)+u(i,j+1)))^2 - (0.5*(u(i,j)+u(i,j-1)))^2)/h;
            advection_y = ((0.25*(u(i,j)+u(i-1,j))*(v(i-1,j)+v(i-1,j+1))) - (0.25*(u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))))/h;
            u_new(i,j) = u(i,j) + dt*(diffusion - advection_x - advection_y + pressure);
        end
    end
    
    % x-momentum eq. - Boundary condition for u velocity
    u_new(1,:) = 2 - u_new(2,:);% in order to have u=1 on the wall
    u_new(n + 1,:) = -u_new(n,:);% in order to have u=0 on the wall
    u_new(2:n,1) = 0;
    u_new(2:n,n) = 0;    
    
    % y-momentum eq.
    for i = 2:n - 1
        for j = 2:n
            pressure = -(p(i,j) - p(i+1,j))/h;
            diffusion = (1/Re)*((v(i+1,j) - 2*v(i,j) + v(i-1,j))/(h*h) + (v(i,j+1) - 2*v(i,j) + v(i,j-1))/(h*h));
            advection_y = ((0.5*(v(i,j)+v(i-1,j)))^2 - (0.5*(v(i,j)+v(i+1,j)))^2)/h;
            advection_x = ((0.25*(u(i,j)+u(i+1,j))*(v(i,j)+v(i,j+1))) - (0.25*(u(i,j-1)+u(i+1,j-1))*(v(i,j)+v(i,j-1))))/h;
            v_new(i,j) = v(i,j) + dt*(diffusion - advection_x - advection_y + pressure);
        end
    end
    
    % y-momentum eq. - Boundary condition for v velocity
    v_new(:,1) = -v_new(:,2); % in order to have v=0 on the wall
    v_new(:,n + 1) = -v_new(:,n); % in order to have v=0 on the wall
    v_new(1,2:n) = 0;
    v_new(n,2:n) = 0;
    
    % Continuity eq.
    for i = 2:n
        for j = 2:n
            p_new(i,j) = p(i,j) - delta*dt*(u(i,j) - u(i,j-1) + v(i-1,j) - v(i,j))/h;
        end
    end
    
    % Continuity eq. - Boundary condition for pressure
    p_new(1,:) = p_new(2,:);
    p_new(n + 1,:) = p_new(n,:);
    p_new(:,1) = p_new(:,2);
    p_new(:,n + 1) = p_new(:,n);
    
    % Continuity error 
    error = 0;
    for i = 2:n - 1
        for j = 2:n - 1
            error = error + abs((u_new(i,j) - u_new(i,j-1) + v_new(i-1,j) - v_new(i,j))/h);
        end
    end
    
    % Error monitoring after every 100 timesteps
%     if(rem(iterations, 100)) == 0
%        figure(1);
%        semilogy(iterations, error, '-ko')
%        hold on
%        xlabel('Iterations')
%        ylabel('Residual Error')
%     end
    
    
    u = u_new;
    v = v_new;
    p = p_new;
%     iterations = iterations + 1;
end

% final values of velocity for our physical points
for i = 1:n
    for j = 1:n
        u_final(i,j) = 0.5*(u(i,j) + u(i+1,j));
        v_final(i,j) = 0.5*(v(i,j) + v(i,j+1));
        p_final(i,j) = 0.25*(p(i,j) + p(i,j+1) + p(i+1,j) + p(i+1,j+1));
    end
end

%% Centreline u validation - Comparison with benchmark
figure(11);
plot(u_final(:,(n+1)/2),1-y, 'LineWidth', 1)

tab_data = xlsread('Re=1000-Ghia', 'A2:B18');
y_ghia = tab_data(:,1);
u_ghia = tab_data(:,2);

figure(11); hold on
plot(u_ghia, y_ghia, 'o', 'LineWidth', 1)
xlabel('u')
ylabel('y')
legend('Numerical solve', 'Ghia etal 1982', 'location', 'southeast')
% legend('Numerical solve', 'location', 'southeast')

figure(33);
plot( x , u_final((n+1)/2,:) , 'LineWidth', 1)
tab_data = xlsread('Re=100-Ghia', 'A2:B18');
x_ghia = tab_data(:,1);
v_ghia = tab_data(:,2);

figure(33); hold on
plot(x_ghia, v_ghia , 'o', 'LineWidth', 1)
xlabel('x')
ylabel('v')
% legend('Numerical solve', 'location', 'southeast')
legend('Numerical solve', 'Ghia etal 1982', 'location', 'southeast')
% %% Contour and vector visuals for 2D domane.
x_dom = ((1:n)-1).*h;
y_dom = 1-((1:n)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);
figure(21);
contourf(X,Y,v_final, 21, 'LineStyle', 'none')
colorbar
colormap('jet')
xlabel('x')
ylabel('y')

figure(22);
hold on
quiver(X, Y, u_final, v_final, 5, 'k')