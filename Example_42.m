clear all
format long


% Parameters for Example 4.2
s = 2;
n = [8, 12];
m = [6, 7];

P1 = [
    0.313, 0.520;
    0.198, 0.493;
    0.245, 0.412;
    0.346, 0.446;
    0.466, 0.528;
    0.397, 0.518;
    0.301, 0.553;
    0.296, 0.473;
    0.299, 0.418
];

P2 = [
    0.299, 0.418;
    0.301, 0.380;
    0.323, 0.342;
    0.328, 0.294;
    0.256, 0.252;
    0.230, 0.272;
    0.173, 0.323;
    0.237, 0.427;
    0.294, 0.278;
    0.350, 0.327;
    0.400, 0.267;
    0.417, 0.296;
    0.396, 0.323
];

r = [1, 3, 1];

t = [0, 0.49, 1];

H(1) = t(2)-t(1);
H(2) = t(3)-t(2);

H = [t(2)-t(1), t(3)-t(2)];

k = [r(1), r(2)];

l = [r(2), r(3)];


% Step 1: Compute C-tables
C1 = C_table(m(1), r(1), r(2));
C2 = C_table(m(2), r(2), r(3));

% Step 2: Compute E1, E2

% Define the variables
syms kappa_x kappa_y lambda1_x lambda1_y lambda2_x lambda2_y lambda3_x lambda3_y

assume(kappa_x, 'real');
assume(kappa_y, 'real');
assume(lambda1_x, 'real');
assume(lambda1_y, 'real');
assume(lambda2_x, 'real');
assume(lambda2_y, 'real');
assume(lambda3_x, 'real');
assume(lambda3_y, 'real');


kappa = [kappa_x; kappa_y];
lambda1 = [lambda1_x; lambda1_y];
lambda2 = [lambda2_x; lambda2_y];
lambda3 = [lambda3_x; lambda3_y];

% this are all the variables
KAPPA = [kappa';lambda1'; lambda2'; lambda3'];

Q1 = sym(zeros(m(1) + 1, 2));

%  eqations (3.3) za Q1
for h = 0:r(2)
    if h == 0
        Q1(m(1)-h+1,:) = KAPPA(1,:);
    else
        suma = 0;
        for j = 1:h
            suma = suma + (-1)^j * nchoosek(h,j) * compute_a_ij(H(1), m(1), j) *KAPPA(1+j,:);
        end
        Q1(m(1)-h+1,:) = KAPPA(1,:) + suma;
    end
end

Q2 = sym(zeros(m(2) + 1, 2));

%  eqations (3.3) za Q2
for h = 0:r(2)
    if h == 0
        Q2(h+1,:) = KAPPA(1,:);
    else
        suma = 0;
        for j = 1:h
            suma = suma + nchoosek(h,j) * compute_a_ij(H(2), m(2), j) *KAPPA(1+j,:);
        end
        Q2(h+1,:) = KAPPA(1,:) + suma;
    end
end

%  eqations (3.1) & (3.2) za Q1 in Q2

for j = 0:r(1)
    if j == 0
        Q1(j+1,:) = nchoosek(n(1),j)./ nchoosek(m(1),j) .* P1(1,:);
    else
        suma2 = 0;
        for h = 0:j-1
            suma2 = suma2 + (-1)^(j+h) * nchoosek(j,h) * Q1(h+1,:);
        end
        Q1(j+1,:) = nchoosek(n(1),j) / nchoosek(m(1),j) * delta(P1, j, 0) - suma2;
    end
end


for j = 0:r(3)
    if j == 0
        Q2(m(2)+1-j,:) = (-1)^j * nchoosek(n(2),j)/nchoosek(m(2),j) * delta(P2, j, n(2)-j);
    else
        suma = 0;
        for h = 1:j
            suma = suma + (-1)^(h) * nchoosek(j,h) .* Q2(m(2) + 1 -j+h,:);
        end
        Q2(m(2)+1-j,:) = (-1)^j * nchoosek(n(2),j)/nchoosek(m(2),j) * delta(P2, j, n(2)-j) - suma;
    end
end


% Lema 2.3 za Q1(3), Q2(4), Q2(5)
for j = k(1)+1:m(1)-l(1)-1
    suma1 = sym(zeros(1,2));
    for h = k(1)+1:m(1)-l(1)-1
        sumaa = 0;
        for i = 0:n(1)
            sumaa = sumaa + nchoosek(n(1),i) / nchoosek(n(1)+m(1), i+h) .* P1(i+1,:);
        end
        vh = 1/(n(1)+m(1)+1) * nchoosek(m(1),h) *(sumaa);
        suma1 = suma1 + C1(j,h).* vh;
    end

    suma2 = sym(zeros(1,2));
    combined_range = [0:k(1), m(1)-l(1):m(1)];
    for h = combined_range
        KJH = compute_K_jh(j, h, m(1), k(1), l(1));
        suma2 = suma2 + KJH.* Q1(h+1,:);
    end
    Q1(j+1,:) = suma1-suma2;
end


for j = k(2)+1:m(2)-l(2)-1
    suma1 = sym(zeros(1,2));
    for h = k(2)+1:m(2)-l(2)-1
        sumaa = 0;
        for i = 0:n(2)
            sumaa = sumaa + nchoosek(n(2),i) / nchoosek(n(2)+m(2), i+h) .* P2(i+1,:);
        end
        vh = 1/(n(2)+m(2)+1) * nchoosek(m(2),h) *(sumaa);
        suma1 = suma1 + C2(j,h).* vh;
    end

    suma2 = sym(zeros(1,2));
    combined_range = [0:k(2), m(2)-l(2):m(2)];
    for h = combined_range
        KJH = compute_K_jh(j, h, m(2), k(2), l(2));
        suma2 = suma2 + KJH.* Q2(h+1,:);
    end
    Q2(j+1,:) = suma1-suma2;
end


p1_1 = P1(:,1); %dim = 1 (h = 1, i = 1)
p2_1= P1(:,2); % dim = 2 (h = 2, i = 1)

p1_2 = P2(:,1); %dim = 1 (h = 1, i = 2)
p2_2= P2(:,2); % dim = 2 (h = 2, i = 2)


q1_1 = Q1(:,1); %dim = 1 (h = 1, i = 1)
q2_1= Q1(:,2); % dim = 2 (h = 2, i = 1)

q1_2 = Q2(:,1); %dim = 1 (h = 1, i = 2)
q2_2= Q2(:,2); % dim = 2 (h = 2, i = 2)


suma = compute_I(n(1),n(1),p1_1,p1_1) + compute_I(m(1),m(1),q1_1,q1_1) - 2.* compute_I(n(1),m(1),p1_1,q1_1);
E1 = suma + compute_I(n(1),n(1),p2_1,p2_1) + compute_I(m(1),m(1),q2_1,q2_1) - 2.* compute_I(n(1),m(1),p2_1,q2_1);
E1 = H(1) * E1;

suma = compute_I(n(2),n(2),p1_2,p1_2) + compute_I(m(2),m(2),q1_2,q1_2) - 2.* compute_I(n(2),m(2),p1_2,q1_2);
E2 = suma + compute_I(n(2),n(2),p2_2,p2_2) + compute_I(m(2),m(2),q2_2,q2_2) - 2.* compute_I(n(2),m(2),p2_2,q2_2);
E2 = H(2) * E2;

E = E1 + E2;

% Step 3: Resi sistem linearnih enacb (3.6) in (3.7)
grad_kappa = gradient(E, [kappa_x, kappa_y]);
grad_lambda1 = gradient(E, [lambda1_x, lambda1_y]);
grad_lambda2 = gradient(E, [lambda2_x, lambda2_y]);
grad_lambda3 = gradient(E, [lambda3_x, lambda3_y]);

eqs = [grad_kappa; grad_lambda1; grad_lambda2; grad_lambda3];

sol = solve(eqs == 0, [kappa_x, kappa_y, lambda1_x, lambda1_y, lambda2_x, lambda2_y, lambda3_x, lambda3_y]);


kappa_opt = [sol.kappa_x, sol.kappa_y];
lambda1_opt = [sol.lambda1_x, sol.lambda1_y];
lambda2_opt = [sol.lambda2_x, sol.lambda2_y];
lambda3_opt = [sol.lambda3_x, sol.lambda3_y];


% Step 4, Step 5, Step 6 in Step 7 (samo vstavi optimalne lambda in kappa) 
% DOBIŠ Q1, Q2
Q1_opt = subs(Q1, [kappa, lambda1, lambda2, lambda3], ...
              [kappa_opt', lambda1_opt', lambda2_opt', lambda3_opt']);
Q2_opt = subs(Q2, [kappa, lambda1, lambda2, lambda3], ...
              [kappa_opt', lambda1_opt', lambda2_opt', lambda3_opt']);

% Računanje napake:
q1_1 = double(Q1_opt(:,1)); % Dim 1 of Q1_opt
q2_1 = double(Q1_opt(:,2)); % Dim 2 of Q1_opt
q1_2 = double(Q2_opt(:,1)); % Dim 1 of Q2_opt
q2_2 = double(Q2_opt(:,2)); % Dim 2 of Q2_opt

E1_part1 = compute_I(n(1), n(1), p1_1, p1_1) + compute_I(m(1), m(1), q1_1, q1_1) - 2 * compute_I(n(1), m(1), p1_1, q1_1);
E1_part2 = compute_I(n(1), n(1), p2_1, p2_1) + compute_I(m(1), m(1), q2_1, q2_1) - 2 * compute_I(n(1), m(1), p2_1, q2_1);
E1 = H(1) * (E1_part1 + E1_part2);

E2_part1 = compute_I(n(2), n(2), p1_2, p1_2) + compute_I(m(2), m(2), q1_2, q1_2) - 2 * compute_I(n(2), m(2), p1_2, q1_2);
E2_part2 = compute_I(n(2), n(2), p2_2, p2_2) + compute_I(m(2), m(2), q2_2, q2_2) - 2 * compute_I(n(2), m(2), p2_2, q2_2);
E2 = H(2) * (E2_part1 + E2_part2);

format short e

fprintf('Optimized E1 - Algorithm 3.1: %.2e\n', double(E1));
fprintf('Optimized E2 - Algorithm 3.1: %.2e\n', double(E2));
fprintf('Optimized E - Algorithm 3.1: %.2e\n', double(E1+E2));

t1_values = linspace(0, 1, 501);
t2_values = linspace(0, 1, 501);

P1_curve = bezier(P1, t1_values);
Q1_curve = bezier(Q1_opt, t1_values);
distances_E1 = sqrt(sum((P1_curve - Q1_curve).^2, 2));
E1_infinity = max(distances_E1);

P2_curve = bezier(P2, t2_values);
Q2_curve = bezier(Q2_opt, t2_values);
distances_E2 = sqrt(sum((P2_curve - Q2_curve).^2, 2));
E2_infinity = max(distances_E2);

fprintf('E1 infinity - Algorithm 3.1: %.2e\n', E1_infinity);
fprintf('E2 infinity - Algorithm 3.1: %.2e\n', E2_infinity);
fprintf('E infinity - Algorithm 3.1: %.2e\n', max(E2_infinity,E1_infinity));


%%%%%% REMARK 3.2.
kappa_x = P2(1,1);
kappa_y= P2(1,2);

% Define the variables
syms lambda1_x lambda1_y lambda2_x lambda2_y lambda3_x lambda3_y t

assume(lambda1_x, 'real');
assume(lambda1_y, 'real');
assume(lambda2_x, 'real');
assume(lambda2_y, 'real');
assume(lambda3_x, 'real');
assume(lambda3_y, 'real');

kappa = [kappa_x; kappa_y];
lambda1 = [lambda1_x; lambda1_y];
lambda2 = [lambda2_x; lambda2_y];
lambda3 = [lambda3_x; lambda3_y];

% this are all the variables
KAPPA = [kappa';lambda1'; lambda2'; lambda3'];

Q1_R = sym(zeros(m(1) + 1, 2));

%  eqations (3.3) za Q1
for h = 0:r(2)
    if h == 0
        Q1_R(m(1)-h+1,:) = KAPPA(1,:);
    else
        suma = 0;
        for j = 1:h
            suma = suma + (-1)^j * nchoosek(h,j) * compute_a_ij(H(1), m(1), j) *KAPPA(1+j,:);
        end
        Q1_R(m(1)-h+1,:) = KAPPA(1,:) + suma;
    end
end

Q2_R = sym(zeros(m(2) + 1, 2));

%  eqations (3.3) za Q2
for h = 0:r(2)
    if h == 0
        Q2_R(h+1,:) = KAPPA(1,:);
    else
        suma = 0;
        for j = 1:h
            suma = suma + nchoosek(h,j) * compute_a_ij(H(2), m(2), j) *KAPPA(1+j,:);
        end
        Q2_R(h+1,:) = KAPPA(1,:) + suma;
    end
end

%  eqations (3.1) & (3.2) za Q1 in Q2

for j = 0:r(1)
    if j == 0
        Q1_R(j+1,:) = nchoosek(n(1),j)./ nchoosek(m(1),j) .* P1(1,:);
    else
        suma2 = 0;
        for h = 0:j-1
            suma2 = suma2 + (-1)^(j+h) * nchoosek(j,h) * Q1_R(h+1,:);
        end
        Q1_R(j+1,:) = nchoosek(n(1),j) / nchoosek(m(1),j) * delta(P1, j, 0) - suma2;
    end
end


for j = 0:r(3)
    if j == 0
        Q2_R(m(2)+1-j,:) = (-1)^j * nchoosek(n(2),j)/nchoosek(m(2),j) * delta(P2, j, n(2)-j);
    else
        suma = 0;
        for h = 1:j
            suma = suma + (-1)^(h) * nchoosek(j,h) .* Q2_R(m(2) + 1 -j+h,:);
        end
        Q2_R(m(2)+1-j,:) = (-1)^j * nchoosek(n(2),j)/nchoosek(m(2),j) * delta(P2, j, n(2)-j) - suma;
    end
end


% Lema 2.3 za Q1(3), Q2(4), Q2(5)
for j = k(1)+1:m(1)-l(1)-1
    suma1 = sym(zeros(1,2));
    for h = k(1)+1:m(1)-l(1)-1
        sumaa = 0;
        for i = 0:n(1)
            sumaa = sumaa + nchoosek(n(1),i) / nchoosek(n(1)+m(1), i+h) .* P1(i+1,:);
        end
        vh = 1/(n(1)+m(1)+1) * nchoosek(m(1),h) *(sumaa);
        suma1 = suma1 + C1(j,h).* vh;
    end

    suma2 = sym(zeros(1,2));
    combined_range = [0:k(1), m(1)-l(1):m(1)];
    for h = combined_range
        KJH = compute_K_jh(j, h, m(1), k(1), l(1));
        suma2 = suma2 + KJH.* Q1_R(h+1,:);
    end
    Q1_R(j+1,:) = suma1-suma2;
end


for j = k(2)+1:m(2)-l(2)-1
    suma1 = sym(zeros(1,2));
    for h = k(2)+1:m(2)-l(2)-1
        sumaa = 0;
        for i = 0:n(2)
            sumaa = sumaa + nchoosek(n(2),i) / nchoosek(n(2)+m(2), i+h) .* P2(i+1,:);
        end
        vh = 1/(n(2)+m(2)+1) * nchoosek(m(2),h) *(sumaa);
        suma1 = suma1 + C2(j,h).* vh;
    end

    suma2 = sym(zeros(1,2));
    combined_range = [0:k(2), m(2)-l(2):m(2)];
    for h = combined_range
        KJH = compute_K_jh(j, h, m(2), k(2), l(2));
        suma2 = suma2 + KJH.* Q2_R(h+1,:);
    end
    Q2_R(j+1,:) = suma1-suma2;
end


p1_1 = P1(:,1); %dim = 1 (h = 1, i = 1)
p2_1= P1(:,2); % dim = 2 (h = 2, i = 1)

p1_2 = P2(:,1); %dim = 1 (h = 1, i = 2)
p2_2= P2(:,2); % dim = 2 (h = 2, i = 2)


q1_1_R = Q1_R(:,1); %dim = 1 (h = 1, i = 1)
q2_1_R= Q1_R(:,2); % dim = 2 (h = 2, i = 1)

q1_2_R = Q2_R(:,1); %dim = 1 (h = 1, i = 2)
q2_2_R= Q2_R(:,2); % dim = 2 (h = 2, i = 2)


suma = compute_I(n(1),n(1),p1_1,p1_1) + compute_I(m(1),m(1),q1_1_R,q1_1_R) - 2.* compute_I(n(1),m(1),p1_1,q1_1_R);
E1_R = suma + compute_I(n(1),n(1),p2_1,p2_1) + compute_I(m(1),m(1),q2_1_R,q2_1_R) - 2.* compute_I(n(1),m(1),p2_1,q2_1_R);
E1_R = H(1) * E1_R;

suma = compute_I(n(2),n(2),p1_2,p1_2) + compute_I(m(2),m(2),q1_2_R,q1_2_R) - 2.* compute_I(n(2),m(2),p1_2,q1_2_R);
E2_R = suma + compute_I(n(2),n(2),p2_2,p2_2) + compute_I(m(2),m(2),q2_2_R,q2_2_R) - 2.* compute_I(n(2),m(2),p2_2,q2_2_R);
E2_R = H(2) * E2_R;

E_R = E1_R + E2_R;

grad_lambda1 = gradient(E_R, [lambda1_x, lambda1_y]);
grad_lambda2 = gradient(E_R, [lambda2_x, lambda2_y]);
grad_lambda3 = gradient(E_R, [lambda3_x, lambda3_y]);

eqs = [grad_lambda1; grad_lambda2; grad_lambda3];

sol_R = solve(eqs == 0, [lambda1_x, lambda1_y, lambda2_x, lambda2_y, lambda3_x, lambda3_y]);

kappa_opt = [kappa_x, kappa_y];
lambda1_opt = [sol_R.lambda1_x, sol_R.lambda1_y];
lambda2_opt = [sol_R.lambda2_x, sol_R.lambda2_y];
lambda3_opt = [sol_R.lambda3_x, sol_R.lambda3_y];

% Step 4, Step 5, Step 6 in Step 7 (samo vstavi optimalne lambda in kappa) 
% DOBIŠ Q1, Q2
Q1_opt_R = subs(Q1_R, [lambda1, lambda2, lambda3], ...
              [lambda1_opt', lambda2_opt', lambda3_opt']);
Q2_opt_R = subs(Q2_R, [lambda1, lambda2, lambda3], ...
              [lambda1_opt', lambda2_opt', lambda3_opt']);

% Računanje napake:
q1_1_R = double(Q1_opt_R(:,1)); % Dim 1 of Q1_opt
q2_1_R = double(Q1_opt_R(:,2)); % Dim 2 of Q1_opt
q1_2_R = double(Q2_opt_R(:,1)); % Dim 1 of Q2_opt
q2_2_R = double(Q2_opt_R(:,2)); % Dim 2 of Q2_opt

E1_part1 = compute_I(n(1), n(1), p1_1, p1_1) + compute_I(m(1), m(1), q1_1_R, q1_1_R) - 2 * compute_I(n(1), m(1), p1_1, q1_1_R);
E1_part2 = compute_I(n(1), n(1), p2_1, p2_1) + compute_I(m(1), m(1), q2_1_R, q2_1_R) - 2 * compute_I(n(1), m(1), p2_1, q2_1_R);
E1_R = H(1) * (E1_part1 + E1_part2);

E2_part1 = compute_I(n(2), n(2), p1_2, p1_2) + compute_I(m(2), m(2), q1_2_R, q1_2_R) - 2 * compute_I(n(2), m(2), p1_2, q1_2_R);
E2_part2 = compute_I(n(2), n(2), p2_2, p2_2) + compute_I(m(2), m(2), q2_2_R, q2_2_R) - 2 * compute_I(n(2), m(2), p2_2, q2_2_R);
E2_R = H(2) * (E2_part1 + E2_part2);

format short e

fprintf('Optimized E1 - Remark 3.2: %.2e\n', double(E1_R));
fprintf('Optimized E2 - Remark 3.2: %.2e\n', double(E2_R));
fprintf('Optimized E - Remark 3.2: %.2e\n', double(E1_R+E2_R));


t1_values = linspace(0, 1, 501);
t2_values = linspace(0, 1, 501);

Q1_curve_R = bezier(Q1_opt_R, t1_values);
distances_E1 = sqrt(sum((P1_curve - Q1_curve_R).^2, 2));
E1_infinity_R = max(distances_E1);

Q2_curve_R = bezier(Q2_opt_R, t2_values);
distances_E2 = sqrt(sum((P2_curve - Q2_curve_R).^2, 2));
E2_infinity_R = max(distances_E2);

fprintf('E1 infinity - Remark 3.2: %.2e\n', E1_infinity_R);
fprintf('E2 infinity - Remark 3.2: %.2e\n', E2_infinity_R);
fprintf('E infinity - Remark 3.2: %.2e\n', max(E2_infinity_R, E1_infinity_R));


disp('Q1 Table (Algorithm 3.1):');
disp(double(Q1_opt));

disp('Q2 Table (Algorithm 3.1):');
disp(double(Q2_opt));

disp('Q1 Table (Remark 3.2):');
disp(double(Q1_opt_R));
disp('Q2 Table (Remark 3.2):');
disp(double(Q2_opt_R));




figure;
hold on;


plot(Q1_curve(:,1), Q1_curve(:,2), 'r--', 'LineWidth', 2); % Reduced Segment 1 (red dashed)
plot(Q2_curve(:,1), Q2_curve(:,2), 'r--', 'LineWidth', 2); % Reduced Segment 2 (red dashed)

plot(Q1_curve_R(:,1), Q1_curve_R(:,2), 'g:', 'LineWidth', 2); % Reduced Segment 1 (green dashed)
plot(Q2_curve_R(:,1), Q2_curve_R(:,2), 'g:', 'LineWidth', 2); % Reduced Segment 2 (green dashed)

plot(P1_curve(:,1), P1_curve(:,2), 'b-', 'LineWidth', 2); % Original Segment 1 (blue solid line)
plot(P2_curve(:,1), P2_curve(:,2), 'b-', 'LineWidth', 2); % Original Segment 2 (blue solid line) 

plot(P1(:,1), P1(:,2), 'bo--', 'MarkerFaceColor', 'b', 'MarkerSize', 1); % Original Segment 1 Control Points
plot(P2(:,1), P2(:,2), 'bo--', 'MarkerFaceColor', 'b', 'MarkerSize', 1); % Original Segment 2 Control Points

plot(Q1_opt(:,1), Q1_opt(:,2), 'ro--', 'MarkerFaceColor', 'r', 'MarkerSize', 1); % Reduced Segment 1 Control Points
plot(Q2_opt(:,1), Q2_opt(:,2), 'ro--', 'MarkerFaceColor', 'r', 'MarkerSize', 1); % Reduced Segment 2 Control Points

plot(Q1_opt_R(:,1), Q1_opt_R(:,2), 'go--', 'MarkerFaceColor', 'g', 'MarkerSize', 1); % Reduced Segment 1 Control Points
plot(Q2_opt_R(:,1), Q2_opt_R(:,2), 'go--', 'MarkerFaceColor', 'g', 'MarkerSize', 1); % Reduced Segment 2 Control Points


scatter(P2(1,1), P2(1,2), 50, 'k', 'filled');


legend('Reduced Segment', '', ...
       'Reduced Segment', '', ...
       'Original Segment', '', ...
       'Original Control Points', '', ...
       'Reduced Control Points ', 'Reduced Control Points 2');
title('Original and Reduced Bézier Curves with Control Points');

hold off;


figure;
hold on;



plot(Q1_curve(:,1), Q1_curve(:,2), 'r--', 'LineWidth', 2); % Reduced Segment 1 (red dashed)
plot(Q2_curve(:,1), Q2_curve(:,2), 'r--', 'LineWidth', 2); % Reduced Segment 2 (red dashed)

plot(Q1_curve_R(:,1), Q1_curve_R(:,2), 'g:', 'LineWidth', 2); % Reduced Segment 1 (green dashed)
plot(Q2_curve_R(:,1), Q2_curve_R(:,2), 'g:', 'LineWidth', 2); % Reduced Segment 2 (green dashed)

plot(P1_curve(:,1), P1_curve(:,2), 'b-', 'LineWidth', 2); % Original Segment 1 (blue solid line)
plot(P2_curve(:,1), P2_curve(:,2), 'b-', 'LineWidth', 2); % Original Segment 2 (blue solid line) 

scatter(P2(1,1), P2(1,2), 50, 'k', 'filled');

legend('Original Segment', '', ...
       'Reduced Segment - Algorithm 3.1', '', 'Reduced Segment - Remark 3.2');
title('Original and Reduced Bézier Curves');
hold off;

