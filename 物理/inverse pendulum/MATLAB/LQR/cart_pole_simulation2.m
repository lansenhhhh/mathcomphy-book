function cart_pole_simulation2
    % Parameters
    mt = 1.0; % Mass of the cart
    mp = 0.2; % Mass of the pendulum
    I = 0.03; % Moment of inertia of the pendulum
    l = 0.5; % Length of the pendulum
    b = 0.1; % Damping coefficient
    g = 9.81; % Acceleration due to gravity

    % Initial conditions
    theta0 = 0.1+pi; % Initial angle
    dtheta0 = 0; % Initial angular velocity
    x0 = 0; % Initial position
    dx0 = 0; % Initial velocity

    init = [x0; dx0; theta0; dtheta0]; % Vector of initial conditions

    % Time span
    tspan = [0 200];
    [T, Y] = ode45(@(t, y) nonlinear(t, y, mt, mp, I, l, b, g), tspan, init);
    figure;
    subplot(2,1,1);
    plot(T, Y(:,1), 'r-');
    title('Nonlinear System: Cart Position (x)');
    xlabel('Time (s)');
    ylabel('Position (m)');

    subplot(2,1,2);
    plot(T, Y(:,3), 'b-'); 
    title('Nonlinear System: Pendulum Angle (theta)');
    xlabel('Time (s)');
    ylabel('Angle (rad)');

    % Solve ODE for linear system
    [T_linear, Y_linear] = ode45(@(t, y) linear(t, y, mt, mp, I, l, b, g), tspan, init);
    figure;
    subplot(2,1,1);
    plot(T_linear, Y_linear(:,1), 'r-'); 
    title('Linear System: Cart Position (x)');
    xlabel('Time (s)');
    ylabel('Position (m)');

    subplot(2,1,2);
    plot(T_linear, Y_linear(:,3)-pi, 'b-'); 
    title('Linear System: Pendulum Angle (theta)');
    xlabel('Time (s)');
    ylabel('Angle (rad)');
end

function dy = nonlinear(t, y, mt, mp, I, l, b, g)
    x = y(1);
    dx = y(2);
    theta = y(3);
    dtheta = y(4);
    u = 0.01*cos(0.5*t);

    M = [mt+mp, mp*l*cos(theta); 
         mp*l*cos(theta), I+mp*l^2];
    C = [b*dx - mp*l*dtheta^2*sin(theta);
         mp*g*l*sin(theta)];
    R = [u; 0];

    ac = M \ (R - C);

    dy = [dx; ac(1); dtheta; ac(2)];
end

function dy = linear(t, y, mt, mp, I, l, b, g)
    x = y(1);
    dx = y(2);
    theta = y(3);
    dtheta = y(4);

    % Linearize around theta = 0
    M = [mt+mp, mp*l; mp*l, I+mp*l^2];
    C = [b*dx; mp*g*l*theta];
    R = [0; 0];

    ac = M \ (-C + R);

    dy = [dx; ac(1); dtheta; ac(2)];
end
