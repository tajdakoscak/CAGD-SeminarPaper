function d = delta(B, j, k)
% forward  difference
    if j == 0
        d = B(k+1, :);
    else
        d = delta(B, j-1, k+1) - delta(B, j-1, k);
    end
end
