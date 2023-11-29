phi = 1;         % Autoregressive parameter
sigma = 1;        % Standard deviation of white noise
num_samples = 100; % Number of samples

% Generate white noise
epsilon = sigma * randn(1, num_samples);

% Generate AR(1) process
X = filter(1, [1, -phi], epsilon);
plot(X)

%%
X(1) = epsilon(1);
for I = 2:100
    X(I) = X(I-1) + epsilon(I);
end
plot(X)