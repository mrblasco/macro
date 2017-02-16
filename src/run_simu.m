function [] = run_simu(params_filename, shocks_filename, output_filename)

% params_filename = './economies_two.csv';
% shocks_filename = './shocks_system.csv';

params = csvread(params_filename, 1);
shocks = csvread(shocks_filename);
fout = fopen(output_filename, 'w');
fout_log = fopen( sprintf('%s.log', output_filename), 'w');
fprintf(fout, '# Economy_id, simulation_id, period, C, Z, K, EulerResidual\n');
fprintf(fout_log, '# Economy_id, score\n');

N = 50;
T = 2000;
K0 = 1.0;
Z0 = 1.0;
numTestCases = size(params, 1);
fprintf('numTestCases: %d\nSettings: N = %d, T = %d, K0 = %f, Z0 = %f\n', numTestCases, N, T, K0, Z0);
shocks = shocks(:, 1:N);

tic
total_score = 0;
for test = 1:numTestCases
    fprintf('Case %d / %d\n', test, numTestCases);
    
    Beta = params(test, 1);
    Eta = params(test, 2);
    Alpha = params(test, 3);
    Delta = params(test, 4);
    Rho = params(test, 5);
    Sigma = params(test, 6);
    
    C = zeros(T+1, N);
    K = K0 * ones(T+1, N);
    W = zeros(T+1, N);
    Z = Z0 * ones(T+1, N);
    eer = zeros(T+1, N);
    
    hasError = false;
    test_score = 0;
    
    SetEconomyParameters(Beta, Eta, Alpha, Delta, Rho, Sigma, N, T);
    
    for t = 1:T
        ti = t+1;
        Z(ti, :) = exp( Rho * log( Z(ti-1, :) ) + Sigma * shocks(ti-1, :) );
    end
    
    for t = 0:T
        ti = t+1;
        
        W(ti, :) = (1-Delta) * K(ti, :) + Z(ti, :) .* K(ti, :).^Alpha;
        
        if (t ~= T)
            C(ti, :) = ConsumptionDecisionRule(K(ti,:), Z(ti,:));
            
            for n = 1:N
                if (C(ti, n) <= 0)
                    fprintf('[ERROR] C[%d][%d] cannot be negative.\n', t, n-1);
                    test_score = -1;
                    hasError = true;
                    break;
                end
                if (C(ti, n) > W(ti, n))
                    fprintf('[ERROR] C[%d][%d] cannot exceed wage.\n', t, n-1);
                    test_score = -2;
                    hasError = true;
                    break;
                end
                if (isinf(C(ti, n)))
                    fprintf('[ERROR] C[%d][%d] must be a real number.\n', t, n-1);
                    test_score = -3;
                    hasError = true;
                    break;
                end
            end
            
            K(ti+1, :) = W(ti, :) - C(ti, :);
        else
            C(ti, :) = W(ti, :);
        end
        
        if hasError
            break
        end
        
        if t > 0
            eer(ti-1, :) = Beta * ((C(ti-1,:) ./ C(ti,:)).^Eta) .* (1 - Delta + Z(ti,:) .* Alpha .* (K(ti,:).^(Alpha - 1))) - 1;
        end
        
    end
    
    if hasError == false
        test_score = 0;
        for t = 0:T-1
            ti = t+1;
            test_score = test_score + abs( sum(eer(ti,:)) );
        end
        test_score = 1000000 / (1 + test_score / (N*T));
        total_score = total_score + test_score;
    end
    
    fprintf('\tscore: %.4f\n', test_score);
    
    if hasError == false
        fprintf('\tSaving...\n');
        for ni = 1:N
            for t = 0:T-1
                ti = t+1;
                fprintf(fout, '%d,%d,%d,%.16g,%.16g,%.16g,%.16g\n', test-1, ni-1, t, C(ti,ni), Z(ti,ni), K(ti,ni), eer(ti,ni));
            end
            fprintf(fout, '%d,%d,%d,%.16g,%.16g,%.16g,\n', test-1, ni-1, t, C(ti,ni), Z(ti,ni), K(ti,ni));
        end
        
    end
    fprintf(fout_log, '%d,\t%.4f\n', test-1, test_score);
end
total_score = total_score / numTestCases;
fprintf('Total score: %.4f\n', total_score);
toc

fclose(fout);
fclose(fout_log);

end
