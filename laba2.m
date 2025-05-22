clear, clc

n = 10;

x = zeros(n,1);
y=zeros(n,1);

for s = 1:15
    c = 10 ^ s;
    x0 = 1 + rand(n,1);
    
    % Собственные значения
    lam1 = 1;
    lam2 = c;
    lamn = 1 + (c- 1) * rand;
    v(1) = 1;
    v(2) = c;
    v(3:n) = lamn;
    v = v(randperm(n));
    D = diag(v);  
    E = eye(n);

    % Ортогональная матрица Q
    w = rand(n, 1);
    Q = E - 2 * w * transpose(w) / (norm(w)) ^ 2;

    % Построение A = Q^T * B * Q
    L = ones(n);
    B = D + triu(L) - diag(L);
    A = Q ^ (-1) * B * Q;

    % Правая часть
    b = A * x0;
    
    % Решение через LDRЙ
    [L, D, R] = ldr_factorization(A);
    LDR = {L, D, R};
    x_solved = solve_ldr(LDR, b);

    % Вычисление ошибок
    veps(s) = norm(x_solved - x0);
    vepsn(s) = norm(A * x_solved - b);
    vc(s) = c;
end 

% Построение графика
figure Name 'Зависимость норм ошибки и невязки от обусловленности'
loglog(vc,veps);
hold on
loglog(vc,vepsn);
legend('Ошибки', 'Невязки');
xlabel('Обусловленность ');
ylabel('Нормы');
grid on;
hold off
    
function [L, D, R] = ldr_factorization(A)
    n = size(A, 1);
    L = eye(n); 
    R = eye(n);
    D = zeros(n);
    
    for m = 1:n
        D(m, m) = A(m, m) - sum(L(m, 1:m-1) .* diag(D(1:m-1, 1:m-1))' .* R(1:m-1, m)');
        for i = m+1:n
            R(m, i) = (A(m, i) - sum(L(m, 1:m-1) .* diag(D(1:m-1, 1:m-1))' .* R(1:m-1, i)')) / D(m, m);
            L(i, m) = (A(i, m) - sum(L(i, 1:m-1) .* diag(D(1:m-1, 1:m-1))' .* R(1:m-1, m)')) / D(m, m);
        end
    end
end

function x = solve_ldr(LDR, b)
    L = LDR{1}; 
    D = LDR{2}; 
    R = LDR{3}; 
    
    % Прямая подстановка для решения LZ = B
    z = forward_substitution(L, b);

    % Диагональная подстановка для решения DY = Z
    y = diagonal_substitution(D, z);

    % Обратная подстановка для решения RX = Y
    x = backward_substitution(R, y);
end

% Прямая подстановка
function x = forward_substitution(A, b)
    n = length(b);
    x = zeros(n, 1);
    for i = 1:n
        x(i) = (b(i) - sum(A(i, 1:i-1) .* x(1:i-1)')) / A(i, i);
    end
end

% Диагональная подстановка
function x = diagonal_substitution(A, b)
    n = length(b);
    x = zeros(n, 1);
    for i = 1:n
        x(i) = b(i) / A(i, i);
    end
end

% Обратная подстановка
function x = backward_substitution(A, b)
    n = length(b);
    x = zeros(n, 1);
    for i = n:-1:1
        x(i) = (b(i) - sum(A(i, i+1:n) .* x(i+1:n)')) / A(i, i);
    end
end

