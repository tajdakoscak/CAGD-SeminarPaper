function C_table = C_table(m, k, l)
    % Input:
    %   m - stopnja Bezierjeve krivulje
    %   k - stopnja odvoda pri t = 0
    %   l - stopnja odvoda pri t = 1
    % Output:
    %   C_table - do indeksov m-l-1; drugje je enaka 0


    C_table = zeros(m -l-1, m -l-1);

    A = @(u) (u - m) * (u - k) * (u + k + 2) / (u + 1);
    B = @(u) u * (u - m - l - 2) * (u - m + l) / (u - m - 1);


    % Compute the first row (j = k+1)
    j = k + 1;
    for h = m - l - 1:-1:k + 1
        if h == m - l - 1
            C_table(j, h) = (1 / nchoosek(m, k + 1)) * ...
                            (1 / nchoosek(m, l + 1)) * ...
                            (-1)^(m - k - l - 2) * factorial(m + k + l + 3) / ...
                            ((m - k - l - 2) * factorial(2 * k + 2) * factorial(2 * l + 2));
        else
            C_table(j, h) = ((h - m) * (h - k) * (h + k + 3)) / ...
                            ((h + 1) * ((h - m)^2 - (l + 1)^2)) * ...
                            C_table(j, h + 1);
        end
    end

    for j = k+1:m-l-2
        for h = k+1: m-k-1
            if h == m-k-1
                C_table(j+1, h) = 1/A(j) * (2*(j-h)*(j+h-m) * C_table(j,h) + ...
                B(h) * C_table(j,h-1) - ...
                B(j) * C_table(j-1,h))

            else

            C_table(j+1, h) = 1/A(j) * (2*(j-h)*(j+h-m) * C_table(j,h) + ...
            B(h) * C_table(j,h-1) + ...
            A(h) * C_table(j,h+1)- ...
            B(j) * C_table(j-1,h))
            end
        end

    end
