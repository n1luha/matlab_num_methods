clear, clc

n = 15; % размерность
d = 2;  % заданный определитель

D = diag(linspace(50, 131, n));
E = eye(n);
w = rand(n, 1);
Q = E - 2 * w * transpose(w) / (norm(w)) ^ 2;
now_det = det(D);
alpha = (d / now_det) ^ (1 / n);
D = D * alpha;

% Сборка симметричной матрицы
A = Q ^ (-1) * D * Q;

% Истинное решение
xt = rand(n, 1);
b = A*xt;

% Подготовка итерационного метода
alpha = 1 / max(diag(D));
C = max(abs(1 - max(abs(diag(D)))), abs(1 - min(abs(diag(D)))));

for t = 1:6
    eps = 0.1^t;
    it = 0;
    x0 = b;
    x1 = A*x0 + b;
    while norm(x1 - x0) > ((1 - C)*eps)/C
        x0 = x1;
        x1 = x0 - alpha*A*x0 + alpha*b;
        it = it + 1;
    end
    xp = x1;

    mist(t) = norm(xp - xt);
    nev(t) = norm(A*xp - b);
    e(t) = eps;
    ited(t) = it;
end

disp("xt = ");
disp(xt);
disp("x = ");
disp(xp);

% Графики
figure('Name', 'Норма фактической ошибки и норма невязки от точности');
loglog(e, mist); hold on;
loglog(e, nev);
legend('фактическая ошибка', 'невязка', 'Location', 'northwest');
xlabel('Точность');
ylabel('Норма');
grid on; grid minor;
title('Анализ точности итерационного метода');
hold off;

figure('Name', 'Число итераций от заданной точности');
semilogx(e, ited);
xlabel('Точность');
ylabel('Количество итераций');
grid on; grid minor;
title('Сходимость итерационного метода');
