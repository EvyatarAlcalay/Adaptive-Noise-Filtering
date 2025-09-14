
%% The code for section 4,5 in Q1

% Define parameters
alpha = 0.9;
Sigma_N_2 = 1;
Sigma_N = sqrt(Sigma_N_2);


%Generates the Z signal
fs = 48000; % sampling frequency in Hz.
T = 10; % duration of the signal in seconds
N = T*fs; % number of samples
G_n = randn(N, 1);
N_n = randn(N, 1);
% Filter the noise to create the WSS signal
b = [1]; % numerator coefficients of the filter
a = [1, -alpha]; % denominator coefficients of the filter
Z_wave = filter(b, a, G_n);
Z_n = Z_wave + N_n;

% Calculate empirical mean and second moment
empiricalMean = mean(Z_n);
empiricalSecondMoment = mean(Z_n.^2);

% Compare to theoretical values
theoreticalMean = 0; % theoretical mean of WSS process is 0
theoreticalSecondMoment = 1/(1-alpha^2) + Sigma_N; % theoretical second moment of WSS process
fprintf('Empirical mean: %.4f\nTheoretical mean: %.4f\n', empiricalMean, theoreticalMean);
fprintf('Empirical second moment: %.4f\nTheoretical second moment: %.4f\n', empiricalSecondMoment, theoreticalSecondMoment);

% 4.2
beta = sqrt(1/(2*empiricalSecondMoment + Sigma_N_2));
fprintf('beta:');
disp(beta);
estimatedBeta = 0.2945;
fprintf('Estimated beta is:');
disp(estimatedBeta)

% Calculate the optimal filters coefficients of orders for every L from 1 to 5
for L= 1:5
    % Calculate the R matrix and w*
    [w_star, ~, ~] = generate_R_p_w_star(alpha, Sigma_N_2, L);
    
    % generate Z_hat
    b = [0 ; w_star];
    a = [1];
    Z_hat = filter(b,a,Z_n);
    
    %calculate the predicition error signal
    e_n = Z_n - Z_hat;

    %play the sound
%     sound(beta*Z_n,fs);
%     sound(beta*e_n,fs);

    %calculate the avarage estimation error
    avg_estimate_err = mean(e_n.^2);
    %calculate the empirical mean of Z form order 2
    avg_Z = mean(Z_n.^2);
    %calculate NRBD
    NRdB = 10 * (log10(avg_Z/avg_estimate_err));
    fprintf(["L=%d, NRdB is:%f, Average Estimation Error:%f\n"], L,NRdB,avg_estimate_err);
end

%% 2A: calculation of the R matrix and its eignables values
% Define parameters
alpha = 0.9;
Sigma_N_2 = 0.5;
L = 4;

% 2B: here we start the steepest decent
[w_star, R, p] = generate_R_p_w_star(alpha, Sigma_N_2, L);

% Calculate eigenvalues of R
[V, D] = eigs(R);
lambda_max = max(diag(D)); % Largest eigenvalue of R

% Define step sizes to try
iter = 100;
mu_values = [0.001, 0.01, 0.1, 0.2];
mu_colores = ["green", "red", "blue", "cyan"];
iter_vec = linspace(0, iter, iter + 1);

% Iterate the steepest descent algorithm for each value of mu
figure;
hold on;
for i = 1:length(mu_values)
    % Initialize weight vector w
    w_n = zeros(L, 1);
    mu = mu_values(i);
    color = mu_colores(i);
    for j = 1:100
        % Calculate gradient of cost function
        grad = p - R*w_n;
        % Update weight vector
        w_n = w_n + mu.*grad;
        % Calculate error norm
        cn_norm = power(norm(w_n - w_star),2);
        w_star_norm = power(norm(w_star),2);
        % Calculate error in dB scale
        error_dB = 10 * log10(cn_norm / w_star_norm);
        % Plot error versus iteration
        plot(j,error_dB,'.', Color=color);
        xlabel('Iteration');
        ylabel('Error (dB)');
        ylim([-100, 20]);
        xlim([0,100]);
%         legend('mu=' + mu);
    end
end


% Explain results
% The error should decrease as the iteration number increases, but the rate of
% convergence will depend on the value of mu. If mu is too small, the algorithm
% will converge slowly. If mu is too large, the algorithm may not converge or
% may converge to a suboptimal solution. The maximal value of mu for which the
% error decreases can be related to the bound 2/lambda_max. If mu is larger than
% this bound, the algorithm may not converge.

%% Q3

% Set the number of samples

F_s = 48 * 10^3;
T = 10;
N = F_s*T;
samplesNum = linspace(0, N, N + 1);
Sigma_N_2 = 0.5;
alpha = 0.9;

% Generate WSS Process
X = generateProcessWSS(N, alpha, Sigma_N);

% Define parameters
L_values = [1, 2, 4];
mu_values = [0.01, 0.001, 0.0001];

% Create a cell array to store w_star values
w_star_values = cell(length(L_values), 1);

% Calculate autocorrelation for given L
for i = 1:length(L_values)
    L = L_values(i);
    [w_star, ~, ~] = generate_R_p_w_star(alpha, Sigma_N_2, L);
    % Store w_star in cell array
    w_star_values{i} = w_star;
end

% Generate process
figure(Name="LMS Error dB");
title("LMS Error dB");

% Iterate over L and mu values
for i = 1:length(L_values)
    L = L_values(i);
    w_n = zeros(L, 1);
    lin_space_vec = zeros(N,1);
    for j = 1:length(mu_values)
        hold on;
        % Initialize filter coefficients
        mu = mu_values(j);
        [w_n, lin_space_vec,e] = LMS_Q3(X', w_n, w_star_values{i}, mu, N, L, lin_space_vec);
        subplot(3,3,(i-1)*3 + j)
        plot(linspace(0,N + 1, N),lin_space_vec);
        xlabel('Iteration');
        ylabel('Relative Wn error [dB]');
        legend("L = " + L + "  , \mu = " + mu);
        NRdB =  10*log10(mean(X.^2) / mean(e.^2));
        title(sprintf('L = %d, mu = %.4f, NRdB = %.4f', L, mu, NRdB));
        grid("minor")
    end
end

%% Q4.1
F_s = 48 * 10^3;
T = 10;
N = F_s*T;
Sigma_N = sqrt(0.5);
L = 2;
alpha = 0.9;
sigma_squared = 0.5;
lambda = 0.99;

% Generate the WSS process
X = generateProcessWSS(N, alpha, Sigma_N);

% Create a cell array to store w_star values
w_star_values = cell(length(L_values), 1)';
disp(w_star_values)

% Calculate autocorrelation for l from 1 to 2

[w_star, R, p] = generate_R_p_w_star(alpha, Sigma_N, L);


e = zeros(1, length(X));
w_n = [0; 0];
errorDB = zeros(N,1);
for delta = [0.1, 1, 10, 100]
    P = eye(L)/delta;
    figure;
    [w_n, errorDB, e] = RLS(X, w_star, lambda, e);
    NRdb_LMS = 10 * log10(mean(X.^2)/ mean(e.^2));
    plot(linspace(0,N-1,N), errorDB);
    title([" Delta = "+ delta + " and NRdb = " + NRdb_LMS]);
end


%% 4.2

lambda = 0.99;
delta = 0.01;
MaxNRdb = 0;
MaxDelta = 0;
for j = 1 : 20
    P = eye(L)/delta;
    w_n = [0; 0];
    [w_n, errorDB, e] = RLS(X, w_star, lambda, e);
    NRdb_LMS = 10 * log10(mean(X.^2)/ mean(e.^2));
    if NRdb_LMS > MaxNRdb
        MaxNRdb = NRdb_LMS;
        MaxDelta = delta;
    end
    delta = delta * 10;
end

disp(MaxDelta);

%% Question 5

% 5(a)

sigma_n = sqrt(0.5);
Sigma_N = 0.5;
alpha = 0.9;
% sound_file = ["vacuumcleaner.wav" "city.wav" "cafe.wav" "airplane.wav"];

for sound_file = ["vacuumcleaner.wav" "city.wav" "cafe.wav" "airplane.wav"];
    
    [signal, Fs] = audioread(sound_file);
    L = 1;
    w_n = [1];
    R = [1/(1-alpha^2)];
    P = [alpha/(1-alpha^2)];
    w_star = inv(R) * P;
    N = length(signal);
    lin_space_vec = zeros(N,1);

    lambda = 0.99;

    % Trivial case LMS
    mu = 0.01;
    [~, ~ , e] = LMS_Q3(signal, w_n, w_star, mu, N, L, lin_space_vec);
    NRdbTestLms = 10 * log10((mean(signal.^2))/mean(e.^2));
    fprintf(['The trivial LMS NRDB for ' + sound_file + ' is ' + ...
        NRdbTestLms + '\n']);

    % Trivial case RLS

    w_n = [1];
    e = zeros(length(signal),1);
    [~, ~, e] = RLS(signal, w_star, lambda, e);
    NRdbTestRls = 10 * log10(mean(signal.^2)/ mean(e.^2));
    fprintf(['The trivial RLS NRDB for ' + sound_file + ' is ' + ...
    NRdbTestLms + '\n']);
end
%% 5(b) + (c) + (d)

sigma_n = sqrt(0.5);
Sigma_N = 0.5;
alpha = 0.9;
M = 10000;

for sound_file = ["vacuumcleaner.wav"];
    
    [signal, Fs] = audioread(sound_file);
    L = 5;
    [w_star, R, p] = generate_R_p_w_star(alpha, Sigma_N, L);
    N = length(signal);
    lin_space_vec = zeros(N,1);

    % LMS plot
    mu = 15;
    w_n = zeros(L, 1);
    [~, e] = LMS(signal, w_n, mu, N, L);
    NRdb_LMS = 10 * log10((mean(signal.^2))/mean(e.^2));
    figure;
    subplot(1,2,1);
    plot(10*log10(movvar(signal,M)) + 1)
    hold on
    plot(10*log10(movvar(e,M)))
    xlabel("Iteration");
    ylabel('Insantaneous Power [dB]');
    title(sprintf(['%s, LMS, L=%d, mu=%f, NR[dB]=%f'],sound_file, L, mu, NRdb_LMS))
    legend('Original Sound', 'Prediction Error', 'Location' , 'Best')

     % RLS plot

    L = 5;
    lambda = 0.9;
    delta = 1;
    P = eye(L)/delta;
%     w_n = zeros(L,1);
    e = zeros(N,1);
    [w_n, errorDB, e] = RLS(signal, w_star, lambda, e);
    NRdb_RLS = 10 * log10(mean(signal.^2)/ mean(e.^2));
    subplot(1,2,2);
    plot(10*log10(movvar(signal,M)) + 1)
    hold on
    plot(10*log10(movvar(e,M)))
    xlabel("Iteration");
    ylabel('Insantaneous Power [dB]');
    title(sprintf(['Sound: %s, LMS, L=%delta, lambda=%f, NR[dB]=%f'],sound_file, L, lambda, NRdb_RLS))
    legend('Original Sound', 'Prediction Error', 'Location' , 'Best')
end

for sound_file = ["city.wav" "cafe.wav" "airplane.wav"];
    
    [signal, Fs] = audioread(sound_file);
    L = 5;
    [w_star, R, p] = generate_R_p_w_star(alpha, Sigma_N, L);
    N = length(signal);
    lin_space_vec = zeros(N,1);

    % LMS plot
    mu = 5;
    w_n = zeros(L, 1);
    [~, e] = LMS(signal, w_n, mu, N, L);
    NRdb_LMS = 10 * log10((mean(signal.^2))/mean(e.^2));
    figure;
    subplot(1,2,1);
    plot(10*log10(movvar(signal,M)) + 1)
    hold on
    plot(10*log10(movvar(e,M)))
    xlabel("Iteration");
    ylabel('Insantaneous Power [dB]');
    title(sprintf(['%s, LMS, L=%d, mu=%f, NR[dB]=%f'],sound_file, L, mu, NRdb_LMS))
    legend('Original Sound', 'Prediction Error', 'Location' , 'Best')

     % RLS plot

    L = 20;
    lambda = 0.99;
    delta = 1;
    P = eye(L)/delta;
    e = zeros(N,1);
    [w_n, errorDB, e] = RLS(signal, w_star, lambda, e);
    NRdb_RLS = 10 * log10(mean(signal.^2)/ mean(e.^2));
    subplot(1,2,2);
    plot(10*log10(movvar(signal,M)) + 1)
    hold on
    plot(10*log10(movvar(e,M)))
    xlabel("Iteration");
    ylabel('Insantaneous Power [dB]');
    title(sprintf(['Sound: %s, LMS, L=%delta, lambda=%f, NR[dB]=%f'],sound_file, L, lambda, NRdb_RLS))
    legend('Original Sound', 'Prediction Error', 'Location' , 'Best')
end

for sound_file = ["cafe.wav"];
    
    [signal, Fs] = audioread(sound_file);
    L = 5;
    [w_star, R, p] = generate_R_p_w_star(alpha, Sigma_N, L);
    N = length(signal);
    lin_space_vec = zeros(N,1);

    % LMS plot
    mu = 22.932;
    w_n = zeros(L, 1);
    [~, e] = LMS(signal, w_n, mu, N, L);
    NRdb_LMS = 10 * log10((mean(signal.^2))/mean(e.^2));
    figure;
    subplot(1,2,1);
    plot(10*log10(movvar(signal,M)) + 1)
    hold on
    plot(10*log10(movvar(e,M)))
    xlabel("Iteration");
    ylabel('Insantaneous Power [dB]');
    title(sprintf(['%s, LMS, L=%d, mu=%f, NR[dB]=%f'],sound_file, L, mu, NRdb_LMS))
    legend('Original Sound', 'Prediction Error', 'Location' , 'Best')

     % RLS plot

    L = 27;
    lambda = 0.99;
    delta = 10;
    P = eye(L)/delta;
    e = zeros(N,1);
    [w_n, errorDB, e] = RLS(signal, w_star, lambda, e);
    NRdb_RLS = 10 * log10(mean(signal.^2)/ mean(e.^2));
    subplot(1,2,2);
    plot(10*log10(movvar(signal,M)) + 1)
    hold on
    plot(10*log10(movvar(e,M)))
    xlabel("Iteration");
    ylabel('Insantaneous Power [dB]');
    title(sprintf(['Sound: %s, LMS, L=%delta, lambda=%f, NR[dB]=%f'],sound_file, L, lambda, NRdb_RLS))
    legend('Original Sound', 'Prediction Error', 'Location' , 'Best')
end

for sound_file = ["airplane.wav"];
    
    [signal, Fs] = audioread(sound_file);
    L = 5;
    [w_star, R, p] = generate_R_p_w_star(alpha, Sigma_N, L);
    N = length(signal);
    lin_space_vec = zeros(N,1);

    % LMS plot
    mu = 1;
    w_n = zeros(L, 1);
    [~, e] = LMS(signal, w_n, mu, N, L);
    NRdb_LMS = 10 * log10((mean(signal.^2))/mean(e.^2));
    figure;
    subplot(1,2,1);
    plot(10*log10(movvar(signal,M)) + 1)
    hold on
    plot(10*log10(movvar(e,M)))
    xlabel("Iteration");
    ylabel('Insantaneous Power [dB]');
    title(sprintf(['%s, LMS, L=%d, mu=%f, NR[dB]=%f'],sound_file, L, mu, NRdb_LMS))
    legend('Original Sound', 'Prediction Error', 'Location' , 'Best')

     % RLS plot

    L = 27;
    lambda = 0.99;
    delta = 10;
    P = eye(L)/delta;
    e = zeros(N,1);
    [w_n, errorDB, e] = RLS(signal, w_star, lambda, e);
    NRdb_RLS = 10 * log10(mean(signal.^2)/ mean(e.^2));
    subplot(1,2,2);
    plot(10*log10(movvar(signal,M)) + 1)
    hold on
    plot(10*log10(movvar(e,M)))
    xlabel("Iteration");
    ylabel('Insantaneous Power [dB]');
    title(sprintf(['Sound: %s, LMS, L=%delta, lambda=%f, NR[dB]=%f'],sound_file, L, lambda, NRdb_RLS))
    legend('Original Sound', 'Prediction Error', 'Location' , 'Best')
end

%% Tests for Q5: Finding a good match for lambda, mu, delta:
for L = [2,5,10,15,20,25,30]
    for sound_file = ["vacuumcleaner.wav" "city.wav" "cafe.wav" "airplane.wav"];
        
        [signal, Fs] = audioread(sound_file);
        Sigma_N = 0.5;
        alpha = 0.9;
        [w_star, R, p] = generate_R_p_w_star(alpha, Sigma_N, L);
        [V, D] = eigs(R);
        lambda_max = max(diag(D)); % Largest eigenvalue of R
        disp(lambda_max);
        N = length(signal);
        M = 10000;
    
        % LMS plot
        for mu = [10^-4, 10^-3, 10^-2, 10^-1,0.5, 1, 5 ,10, 15, 20, 22, 22.9323]
            w_n = zeros(L, 1);
            disp(mu);
            [~, e] = LMS(signal, w_n, mu, N, L);
            NRdb_LMS = 10 * log10((mean(signal.^2))/mean(e.^2));
            figure;
            grid on;
            plot(10*log10(movvar(signal,M)) + 1)
            hold on
            plot(10*log10(movvar(e,M)))
            xlabel("Iteration");
            ylabel('Insantaneous Power [dB]');
            title(sprintf(['%s, LMS, L=%d, mu=%f, NR[dB]=%f'],sound_file, L, mu, NRdb_LMS))
            legend('Original Sound', 'Prediction Error', 'Location' , 'Best')
            grid on;
        end
    
         % RLS plot
    %      delta = [100, 10, 1, 10^-1, 10^-2, 10^-3, 10^-4]
        for lambda = [0.99, 10^-1, 10^-2, 10^-3, 10^-4, 10^-5, 10^-6, 0.5, 10, 50]
            for delta = [1]
                fprintf(['lambda=%f, delta=%f'],lambda, delta);
                P = eye(L)/delta;
                e = zeros(N,1);
                [w_n, errorDB, e] = RLS(signal, w_star, lambda, e);
                NRdb_RLS = 10 * log10(mean(signal.^2)/ mean(e.^2));
                figure;
                plot(10*log10(movvar(signal,M)) + 1)
                hold on
                plot(10*log10(movvar(e,M)))
                xlabel("Iteration");
                ylabel('Insantaneous Power [dB]');
                title(sprintf(['%s, RLS, L=%d, delta=%f, lambda=%f, NR[dB]=%f'],sound_file, L, delta, lambda, NRdb_RLS))
                legend('Original Sound', 'Prediction Error', 'Location' , 'Best')
                grid on;
            end
        end
    end
end

%% Helper Functions:

function X = generateProcessWSS(N, alpha, Sigma_N)
    % Initialize the signal arrays
    X = zeros(N,1);
    signal = zeros(N,1);

    % Set the initial values of X and Y
    X(1) = randn(1)/sqrt(1-alpha^2);

    % Generate the signal
    for n = 2:N

        % Generate the white noise process Nn
        Nn = randn(1)*Sigma_N;

        % Generate the process Gn
        Gn = randn(1);

        % Generate the process zn
        signal(n) = alpha*signal(n-1) + Gn;

        % Generate the desired signal Xn
        Xn = signal(n) + Nn;

        % Update the signal arrays
        X(n) = Xn;
    end
end

function [w_star, R, p] = generate_R_p_w_star(alpha, Sigma_N, L)
    vec = [];
    for l = 0:L-1
        if l==0
            rz = 1/(1 - alpha^2) + Sigma_N;
        else
            rz = (alpha^abs(l))/(1 - alpha^2);
        end
        vec = vertcat(vec, rz);
    end
    R = toeplitz(vec);
    
    p = [];
    for l = 1:L
        p_i = (alpha^abs(l))/(1 - alpha^2);
        p = vertcat(p, p_i);
    end
    w_star = inv(R)*p;
end

function [w_n, lin_space_vec, e] = LMS_Q3(X, w_n, w_star, mu, N, L, lin_space_vec)
    e = zeros(N,1);
    for n = 1:N - L
        % get Y_n
        Y_n = X(n:n+L-1);

        % Get the current estimator
        X_n_hat = flip(Y_n)*w_n;

        % Calculate output and error
        e(n) = X(n+L) - X_n_hat;
        % Calculate W_n+1
        w_n = w_n + mu.*flip(Y_n') * e(n);

        % Calculate C_n^2
        c_n_norm = power(norm(w_n - w_star),2);
        w_star_norm = power(norm(w_star),2);
        c_n_squared = 10*log10(c_n_norm / w_star_norm);
        lin_space_vec(n) = c_n_squared;
    end
end

function [w_n, e] = LMS(X, w_n, mu, N, L)
    e = zeros(N,1);
    for n = 1:N - L
        % get Y_n
        Y_n = X(n:n+L-1);

        % Get the current estimator
        X_n_hat = flip(Y_n')*w_n;

        % Calculate output and error
        e(n) = X(n+L) - X_n_hat;
        % Calculate W_n+1
        w_n = w_n + mu.*flip(Y_n) * e(n);
    end
end

function [w_n, errorDB, e] = RLS(X, w_star, lambda, e)
    N = length(X);
    L = length(w_star);
    w_n = zeros(L, 1);
    P = eye(L) / lambda;
    errorDB = zeros(N, 1);
    for i = L + 1 : N
        Y = X(i-L:i-1);
        X_hat = w_n.' * Y;
        e(i) = X(i) - X_hat;
        K = ((P * Y) / lambda) / (1 + (1 / lambda) * Y.' * P * (1 / lambda) * Y);
        w_n = w_n + K * e(i);
        P = P / lambda - (K * Y.' * P) / lambda;
        errorDB(i) = 10*log10((norm(w_star - w_n))^2/(norm(w_star))^2);
    end
end

function [w_n, lin_space_vec, e] = LMS_Q3_NEW(X, w_n, w_star, mu, N, L, lin_space_vec)
    e = zeros(N,1);
    for n = L+1:N
        % get Y_n
        Y_n = X(n-1:-1:n-L,1);

        % Get the current estimator
        X_n_hat = w_n'*Y_n;

        % Calculate output and error
        e(n) = X(n) - X_n_hat;
        % Calculate W_n+1
        w_n = w_n + mu.*Y_n.*e(n);

        % Calculate C_n^2
        c_n_norm = norm(w_n - w_star)^2;
        w_star_norm = norm(w_star)^2;
        c_n_squared = 10*log10(c_n_norm / w_star_norm);
        lin_space_vec(n) = c_n_squared;
    end
end
