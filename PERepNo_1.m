function RepNo = PERepNo(paramset)

    % Reproduction number
    beta = paramset(1);
    sigma = paramset(2); 
    alpha = paramset(3);
    rho = paramset(4);
    phi = paramset(5);
    tau = paramset(6);
    epsilon = paramset(7);
    kappa = paramset(8);
    lambda = paramset(9);
    delta = paramset(10);
    gamma = paramset(11);
    theta = paramset(12);
    mu = paramset(13);
    zeta = paramset(14);

    % Basic Rep. No equation:
    W = alpha+zeta; %A
    X = gamma+mu+zeta;  %B
    Y = gamma+zeta-gamma*tau+phi*tau;  %C

    RepNo = alpha*beta*(gamma*rho+phi*tau+rho*zeta-gamma*rho*tau ...
        -sigma*(rho-1)*X)/(W*X*Y);

end