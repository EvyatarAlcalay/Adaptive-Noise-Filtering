
function znext = adaptivepredict(zvec)
    Z = zeros(length(zvec), 1);
    for i = 1 : length(zvec)
        Z(i) = zvec(i);
    end
    
    % Initialization of parameters
    bestNRdB = -inf;
    N = length(Z);
    znext = 0;
    
    % Run RLS for 3 Times for different L's and lambda's
    L_values_RLS = [10 20 30];
    lambda_vals = [0.01, 0.5, 0.999];
    for i = 1:length(L_values_RLS)
        delta = 1;
        disp(lambda_vals(i));
        [znext, bestNRdB] = RLS(Z, L_values_RLS(i), lambda_vals(i), delta, bestNRdB, znext);
    end

    mu_vals = [0.85 9 0.01];
    L_values_LMS = [5 10 15];
    for i = 1:length(L_values_LMS)
        disp(mu_vals(i));
        w_n = zeros(1, L_values_LMS(i));
        [znext, bestNRdB] = LMS(Z, w_n, mu_vals(i), N, L_values_LMS(i), bestNRdB, znext);
    end
    fprintf(["final Best estimator is:"]);
    disp(znext)
    fprintf(["Best NRdB value is:"]);
    disp(bestNRdB);
end

% Helper Functions:

function [znext, bestNRdB] = RLS(X, L, lambda, delta, bestNRdB, znext)
    e = zeros(length(X),1);
    N = length(X);
    w_n = zeros(L, 1);
    P = eye(L) / delta;
    for i = L + 1 : N
        Y_n = X(i-L:i-1);
        X_hat = w_n.' * Y_n;
        e(i) = X(i) - X_hat;
        K = ((P * Y_n) / lambda) / (1 + (1 / lambda) * Y_n.' * P * (1 / lambda) * Y_n);
        w_n = w_n + K * e(i);
        P = P / lambda - (K * Y_n.' * P) / lambda;
    end
    currNRdB = 10 * log10(mean(X.^2)/ mean(e.^2));
    % finds the best NRdb and updates znext if a better NRdB val was found
    if currNRdB > bestNRdB
        bestNRdB = currNRdB;
        znext = w_n.'* Y_n;
        disp(znext);
    end
end


function [znext, bestNRdB] = LMS(X, w_n, mu, N, L, bestNRdB, znext)
    e = zeros(N,1);
    for n = 1:N - L
        % get Y_n
        Y_n = X(n:n+L-1);

        % Get the current estimator
        X_n_hat = w_n * flip(Y_n);

        % Calculate output and error
        e(n) = X(n+L) - X_n_hat;
        
        % Calculate W_n+1
        w_n = w_n + mu.*flip(Y_n') * e(n);
    end
    currNRdB = 10 * log10(mean(X.^2)/ mean(e.^2));
    % Finds the best NRdb and updates znext if a better value
    % was found
    if currNRdB > bestNRdB
        bestNRdB = currNRdB;
        znext = w_n.'* Y_n;
        disp(znext);
    end
end
