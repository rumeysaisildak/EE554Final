%Assignment - 5

clear all;
clc;

snr_db = (-20 : 2 : 20);
Nt = 4;
P = 1;
beta = 0;
iterations = 5e5;
SecrecyCapacity = zeros(3, length(snr_db));
for i = 1 : 3
    for j = 1 : length(snr_db)
        fprintf("w = %d, SNR dB = %d has just started.\n", i, snr_db(j))
        N0 = 1 / 10^(snr_db(j) / 10);

        Cs_sum = 0;
        for k = 1 : iterations
            h = randn(Nt, 1) + 1i * randn(Nt, 1);
            h = h ./ sqrt(2);

            g = randn(Nt, 1) + 1i * randn(Nt, 1);
            g = g ./ sqrt(2);
            
            Ph = eye(Nt) - h * h' / norm(h)^2;
            Sigma = beta * Ph;
            Pz = trace(Sigma) / Nt;
            
            z = randn(Nt, 1) + 1i * randn(Nt, 1);
            z = z * sqrt(Pz) ./ sqrt(2);

            if i == 1
                w = ones(Nt, 1);
            elseif i == 2
                w = g;
            else
                w = h;
            end
            
            Cb = log2(1 + (abs(h' * w)^2 * P) / (abs(h' * z)^2 + N0));
            Ce = log2(1 + (abs(g' * w)^2 * P) / (abs(g' * z)^2 + N0));
            
            Cs_instant = Cb - Ce;
            if Cs_instant >= 0
                Cs_sum = Cs_sum + Cs_instant;
            end
        end
        Cs_ort = Cs_sum / iterations;
        SecrecyCapacity(i, j) = Cs_ort;
        fprintf("w = %d, SNR dB = %d has just finished.\n\n", i, snr_db(j))
    end
end

semilogy(snr_db, SecrecyCapacity(1, :), '-^', 'LineWidth', 3, ...
                                              'Color', [0.4660 0.6740 0.1880], ...
                                              'MarkerEdgeColor', 'k', ...
                                              'MarkerFaceColor', '#77DD77', ...
                                              'MarkerSize', 12);
hold on

semilogy(snr_db, SecrecyCapacity(2, :), '-h', 'LineWidth', 3, ...
                                              'Color', [0.9290 0.6940 0.1250], ...
                                              'MarkerEdgeColor', 'k', ...
                                              'MarkerFaceColor', [0.4940 0.1840 0.5560], ...
                                              'MarkerSize', 18);
hold on

semilogy(snr_db, SecrecyCapacity(3, :), '-p', 'LineWidth', 3, ...
                                              'Color', [0.6350 0.0780 0.1840], ...
                                              'MarkerEdgeColor', 'k', ...
                                              'MarkerFaceColor', [0.3010 0.7450 0.9330], ...
                                              'MarkerSize', 20);
xlabel('1 / \sigma^{2}');
ylabel('Secrecy Capacity');
title('Secrecy Capacity for w = [1 1 1 1]^{T}, w = g, and w = h');
legend("w = [1 1 1 1]^{T}","w = g","w = h");
grid;                                          




