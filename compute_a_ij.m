function a_ij = compute_a_ij(h_i, m_i, j)
    % compute_a_ij calculates the coefficient a_ij
    % Input:
    %   h_i - scalar parameter
    %   m_i - degree of the curve
    %   j - index for the coefficient
    % Output:
    %   a_ij - computed coefficient
    
    % Calculate the numerator h_i^j
    numerator = h_i^j;

    a_ij = numerator / shifted_factorial(m_i-j+1,j);
    
end
