clear; clc;
n = 11;

% Построение симметричной матрицы
E = eye(n);
w = rand(n, 1); 
lam = linspace(200, 1, n)';
lammax = max(lam);
D = diag(lam);

C = triu(rand(n));  % Верхняя треугольная матрица
C(1:n+1:end) = 0;    % Обнуляем диагональ
W = C + D;

Q = E - 2 * (w * w') / (norm(w)^2);  % Ортогональная матрица Хаусхолдера
A = Q' * W * Q;  % Симметричная матрица

% Метод степенных итераций
normir = zeros(1, 8);  % ошибка по λ
nev = zeros(1, 8);     % невязка
e = zeros(1, 8);       % точности
iter = zeros(1, 8);    % итерации

for k = 1:8
    x_prev = ones(n, 1); 
    x_prev = x_prev / norm(x_prev);
    lam_prev = 1;
    eps = 0.1^k;
    it = 0;

    while true
        x_next = A * x_prev;
        lam_next = (x_next' * x_prev) / (x_prev' * x_prev);
        
        if abs(lam_next - lam_prev) < eps
            break;
        end

        x_next = x_next / norm(x_next);  % Евклидова нормировка
        x_prev = x_next;
        lam_prev = lam_next;
        it = it + 1;
    end

    normir(k) = abs(lam_next - lammax);  % Ошибка по собственному значению
    nev(k) = norm(A * x_next - x_next * lam_next);  % Невязка
    e(k) = eps;
    iter(k) = it;
end

disp("Последняя собственная невязка:");
disp(nev(end));
disp("Последняя ошибка в собственном значении:");
disp(normir(end));

% Графики
figure('Name', 'Норма фактической ошибки и норма невязки от точности');
loglog(e, normir); hold on;
loglog(e, nev);
legend('Ошибка', 'Невязка', 'Location', 'northwest');
xlabel('Точность');
ylabel('Норма');
grid on; grid minor;
title('Сходимость метода степенных итераций');
hold off;

figure('Name', 'Число итераций от точности');
semilogx(e, iter);
xlabel('Точность');
ylabel('Количество итераций');
grid on; grid minor;
title('Итерации метода степенных итераций');
