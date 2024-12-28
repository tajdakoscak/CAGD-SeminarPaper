function C_table = C_table(m, k, l)
    % Compute the C-table using Algorithm 2.1
    % Input:
    %   m - degree of the Bézier curve
    %   k - derivative degree at t = 0
    %   l - derivative degree at t = 1
    % Output:
    %   C_table - coefficients matrix

    % Initialize C-table with zeros
    C_table = zeros(m +1, m + 1);

    % Helper functions for A(u) and B(u)
    A = @(u) (u - m) * (u - k) * (u + k + 2) / (u + 1);
    B = @(u) u * (u - m - l - 2) * (u - m + l) / (u - m - 1);

    j = k + 1;
    for h = m - l - 1:-1:k + 1
        if h == m - l - 1
            C_table(j, h) = (1 / nchoosek(m, k + 1)) * ...
                            (1 / nchoosek(m, l + 1)) * ...
                            (-1)^(m - k - l - 2) * factorial(m + k + l + 3) / ...
                            (factorial(m - k - l - 2) * factorial(2 * k + 2) * factorial(2 * l + 2));
        else
            C_table(j, h) = ((h - m) * (h - k) * (h + k + 3)) / ...
                            ((h + 1) * ((h - m)^2 - (l + 1)^2)) * ...
                            C_table(j, h + 1);
        end
    end

    for j = k + 1:m - l - 2
        for h = k + 1:m - l - 1
            term1 = 0; term2 = 0; term3 = 0; term4 = 0;

            if h >= k + 1 && h <= m - l - 1
                term1 = 2 * (j - h) * (j + h - m) * C_table(j, h);
            end
            if h > k + 1
                term2 = B(h) * C_table(j, h - 1);
            end
            if h < m - l - 1
                term3 = A(h) * C_table(j, h + 1);
            end
            if j > k + 1
                term4 = -B(j) * C_table(j - 1, h);
            end

            C_table(j + 1, h) = (1 / A(j)) * (term1 + term2 + term3 + term4);
        end
    end
end
