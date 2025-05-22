clear; clc;

%% МПД алгебраическое уравнение
f = @(x) 2*x.^4 + 8*x.^3 + 8*x.^2 - 1;
a0 = 0.1;
b0 = 0.6;
n1 = fzero(f, [a0, b0]);

for k = 1:6
    a = a0; b = b0;
    my_eps = 0.1^k;
    fa = f(a);
    it = 0;
    while (b - a) > my_eps
        it = it + 1;
        c = (a + b)/2;
        if f(c) * fa < 0
            b = c;
        else
            a = c;
        end
    end
    wc(k) = abs(c - n1);
    we(k) = my_eps;
    wi(k) = it;
end

disp("МПД алгебраическое:");
disp(["Корень:" c, "Итераций:" it]);

%% Метод хорд — алгебраическое уравнение
f = @(x) 2*x.^4 + 8*x.^3 + 8*x.^2 - 1;
a0 = 0.1;
b0 = 0.6;
df = @(x) 8*x.^3 + 24*x.^2 + 16*x;
M1 = df(b0);
m1 = df(a0);
n2 = fzero(f, [a0, b0]);

for k = 1:6
    my_eps = 0.1^k;
    a = b0;
    c = a0;
    it = 0;
    while ((M1 - m1)/M1) * abs(a - c) >= my_eps
        a = c;
        c = a - (f(a) * (b0 - a)) / (f(b0) - f(a));
        it = it + 1;
    end
    xc(k) = abs(c - n2);
    xe(k) = my_eps;
    xi(k) = it;
end

disp("Метод хорд — алгебраическое:");
disp(["Корень:" c, "Итераций:" it]);

%% МПД трансцендентное уравнение
f = @(x) (x - 3)*cos(x) - 1;
a0 = -5;
b0 = -4;
n3 = fzero(f, [a0, b0]);

for k = 1:6
    a = a0; b = b0;
    my_eps = 0.1^k;
    fa = f(a);
    it = 0;
    while (b - a) > my_eps
        it = it + 1;
        c = (a + b)/2;
        if f(c) * fa < 0
            b = c;
        else
            a = c;
        end
    end
    yc(k) = abs(c - n3);
    ye(k) = my_eps;
    yi(k) = it;
end

disp("МПД трансцендентное:");
disp(["Корень:" c, "Итераций:" it]);

%% Метод хорд — трансцендентное уравнение
f = @(x) (x - 3)*cos(x) - 1;
a0 = -5;
b0 = -4;
df = @(x) -x*sin(x) + 3*sin(x) + cos(x);
M1 = df(a0);
m1 = df(b0);
n4 = fzero(f, [a0, b0]);

for k = 1:6
    my_eps = 0.1^k;
    a = b0;
    c = a0;
    it = 0;
    while ((M1 - m1)/M1) * abs(a - c) >= my_eps
        a = c;
        c = a - (f(a) * (b0 - a)) / (f(b0) - f(a));
        it = it + 1;
    end
    zc(k) = abs(c - n4);
    ze(k) = my_eps;
    zi(k) = it;
end

disp("Метод хорд — трансцендентное:");
disp(["Корень:" c, "Итераций:" it]);

%% Построение графиков
figure;
subplot(2,1,1);
title('Зависимость ошибки от точности')
loglog(we, wc, xe, xc, ye, yc, ze, zc, xe, xe);
legend({'МПД алг.', 'МХ алг.', 'МПД трансц.', 'МХ трансц.', 'бисектриса'},'Location','southeast');
xlabel('Точность');
ylabel('Ошибка');
grid on; grid minor;

subplot(2,1,2);
title('Зависимость числа итераций от точности')
semilogx(we, wi, xe, xi, ye, yi, ze, zi);
legend({'МПД алг.', 'МХ алг.', 'МПД трансц.', 'МХ трансц.'},'Location','northeast');
xlabel('Точность');
ylabel('Итерации');
grid on; grid minor;
