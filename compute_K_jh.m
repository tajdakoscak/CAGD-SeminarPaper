function K_jh = compute_K_jh(j, h, m, k, l)
t1 = nchoosek(m,h) / nchoosek(m,j);
t2 = (-1)^(j-k-1) * shifted_factorial(k+1-h,m-k-l-1) * shifted_factorial(k+2, h) * shifted_factorial(l+2,m-h);
t3 = (j-h) * factorial(j-k-1) * factorial(m-l-j-1) * shifted_factorial(k+2,j) * shifted_factorial(l+2,m-j);
K_jh = t1 * t2 / t3;
end
