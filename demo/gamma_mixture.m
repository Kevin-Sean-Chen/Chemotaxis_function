% EM algorithm for mixture of gamma distributions,
% which is a model for displacement (speed vectors) in worm movements

%% load some displacement data
x = alldis(1:100000)  /33/(5/14);
mask = alltrials(1:100000);
mask(find(x>6)) = NaN;
figure
hist(x,100)
N = length(x);

%%% https://minerva.it.manchester.ac.uk/~saralees/gammapaper.pdf

%%
% Initialize the parameters of the mixture of gamma distributions
n_mixtures = 3;  % Number of mixture components
% alpha_est = [1 7 ];  % Initial shape parameters
% beta_est = [.5 .2 ];  % Initial rate parameters
alpha_est = [1. 5 6];
beta_est = [.1 .1 .1];
lamb_est = ones(1,n_mixtures)/n_mixtures;  % Initial component probabilities
z_nk = zeros(N, n_mixtures);  % latent z
% EM algorithm
maxIter = 100;  % Maximum number of iterations
log_likelihood = zeros(maxIter, 1);
for iter = 1:maxIter
    % Expectation step (E-step)
    for k = 1:n_mixtures
        z_nk(:, k) = lamb_est(k) * gampdf(x, alpha_est(k), beta_est(k));
    end
    z_nk = z_nk ./ sum(z_nk, 2);
    
    % Maximization step (M-step)
    lamb_est = sum(z_nk,1) / N;
    
%     for k = 1:n_mixtures
%         
%         %%% numerical method
%         lfun = @(xx)nLL_weighted_gamma(xx, z_nk(:,k), x, mask);  % objective function
%         opts = optimset('display','iter');
%         LB = zeros(1,2);
%         UB = ones(1,2)*300;%[inf inf inf inf];
%         prs0 = rand(1,2);%[alpha_est(k), beta_est(k)];%, alpha_est(2), beta_est(2)];
%         [param_est,fval,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon(lfun,prs0,[],[],[],[],LB,UB,[],opts);  % constrained optimization
%         alpha_est(k) = param_est(1);
%         beta_est(k) = param_est(2);
%     end    
        %%% closed-form method(?)
    for k = 1:n_mixtures
        sum_z = sum(z_nk(:,k));
        sum_zx = sum(z_nk(:,k)'.*x);
        sum_zxlx = sum(z_nk(:,k)'.*x.*log(x));
        sum_zlx = sum(z_nk(:,k)'.*log(x));
        alpha_est(k) = sum_z*sum_zx / (sum_z*sum_zxlx - sum_zlx*sum_zx);
        beta_est(k) = (sum_z*sum_zxlx - sum_zlx*sum_zx) / sum_z^2;
    end
    
%     alpha_est = [param_est(1), param_est(3)];
%     beta_est = [param_est(2), param_est(4)];
    
    % Calculate the log-likelihood
    log_likelihood(iter) = nansum(nansum(log(z_nk),1));  %sum(log(sum(z_nk, 2))); %
    
    % Check convergence
    if iter > 1 && abs(log_likelihood(iter) - log_likelihood(iter-1)) < 1e-7
        break;
    end
    disp(['iteration',num2str(iter)])
end

% Plot the data and the estimated mixture of gamma distributions
figure
histogram(x, 'Normalization', 'pdf')
hold on

x_vals = linspace(min(x), max(x), 1000);
y_vals = zeros(numComponents, length(x_vals));
for k = 1:n_mixtures
    y_vals(k, :) = lamb_est(k) * gampdf(x_vals, alpha_est(k), beta_est(k));
    plot(x_vals, y_vals(k, :), 'LineWidth', 2); hold on
end
plot(x_vals, sum(y_vals,1))

% Customize the plot
xlabel('x')
ylabel('Probability Density')
title('Mixture of Gamma Distributions Estimation')
legend('Data', 'Component 1', 'Component 2')

hold off


%%
% alpha_est = [.95 3.5];
% beta_est = [.8 .5];
% lamb_est = [.4, .6];
% figure
% histogram(x, 'Normalization', 'pdf')
% hold on
% 
% x_vals = linspace(min(x), max(x), 1000);
% y_vals = zeros(numComponents, length(x_vals));
% for k = 1:numComponents
%     y_vals(k, :) = lamb_est(k) * gampdf(x_vals, alpha_est(k), beta_est(k));
%     plot(x_vals, y_vals(k, :), 'LineWidth', 2); hold on
%     plot(x_vals, sum(y_vals,1))
% end
% 
% % Customize the plot
% xlabel('x')
% ylabel('Probability Density')
% title('Mixture of Gamma Distributions Estimation')
% legend('Data', 'Component 1', 'Component 2')
% 
% hold off

%% nLL function
%Negative Log-likelihood for chemotaxis with kernels
function [NLL] = nLL_weighted_gamma(THETA, alpha, dis, mask)

    %%% loading parameters
    k1 = THETA(1);
    theta1 = THETA(2);
%     k2 = THETA(3);
%     theta2 = THETA(4);
    
    %%% log-likelihood
    N = length(dis);
    dis = mask(1:end).*dis;  % remove masked elements
    L1 = (k1-1)*nansum(alpha(:)'.*log(dis)) - nansum(alpha(:)'.*dis)/theta1 - nansum(alpha(:)*k1*log(theta1) + alpha(:)*log(gamma(k1)));
%     L2 = (k2-1)*nansum(alpha(:,2)'.*log(dis)) - nansum(alpha(:,2)'.*dis)/theta2 - 1*(N*k2*log(theta2) + N*log(gamma(k2)));
    %%% log_likelihood = sum((k-1)*log(data) - data/theta - k*log(theta) -gammaln(k));  %nansum(alpha(:,2))
    
    marginalP = L1;% + L2;
    NLL = -marginalP;
end