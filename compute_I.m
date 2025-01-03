function I = compute_I(N, M, x, y)
    I = 0;
    for i = 0:N
        for j = 0:M
            I = I + nchoosek(N, i) * nchoosek(M, j) * x(i+1) * y(j+1) / nchoosek(N+M, i+j);
        end
    end
    I = I / (N + M + 1);
end