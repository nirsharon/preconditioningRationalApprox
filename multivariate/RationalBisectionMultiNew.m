%% Part 1: Rational Bisection Code and Model Construction
function [p, q, g, Err] = RationalBisectionMultiNew(f,T,T1,T2,d1,eps,N1,N2,xx,mm,f1)
    Vn = [];
    Vm = [];

    for i = 1:length(T)
        if rem(T(i), d1) == 0
            x = T1(d1);
            y = T2(T(i)/d1);
        else
            x = T1(rem(T(i), d1));
            y = T2(1 + (T(i) - rem(T(i), d1)) / d1);
        end

        Vn(i,:) = chebyshevBasis2D(x, y, N1, xx, mm);
        Vm(i,:) = chebyshevBasis2D(x, y, N2, xx, mm, true);
    end

    n = size(Vn,2);
    m = size(Vm,2);

    L = 0;
    U = max(abs(f(T)));
    u = 0;

    while U - L > eps
        u = (U + L)/2;
        if checkValue(f, n, m, T, u, Vn, Vm)
            U = u;
        else
            L = u;
        end
    end

    [p, q, ~] = createModel(f, n, m, T, U, Vn, Vm);

    G1 = Vn * p;
    G2 = Vm * q;
    g = G1 ./ (1 + G2);

    % errors or deviations on the selected domain
    Err = zeros(numel(T),1);
    F = f(T);
    for r = 1:numel(T)
        Err(r) = F(r) - g(r);
    end

    g = reshape(g,[length(T1),length(T2)]);
    Err = reshape(Err,[length(T1),length(T2)]);

    %F = f(T);
    %Err = F - g;
    %g = reshape(g, [length(T1), length(T2)]);
    %Err = reshape(Err, [length(T1), length(T2)]);
end

function V = chebyshevBasis2D(x, y, N, xx, mm, skipZero)
    if nargin < 6
        skipZero = false;
    end
    x1 = (x - min(xx)) / (max(xx) - min(xx));
    y1 = (y - min(mm)) / (max(mm) - min(mm));
    V = [];
    for n = (skipZero + 0):N
        for k = 0:n
            V(end+1) = cos(k * acos(x1)) * cos((n - k) * acos(y1));
        end
    end
end

function C = checkValue(f, n, m, T, Lk, Vn, Vm)
    [~, ~, u] = createModel(f, n, m, T, Lk, Vn, Vm);
    C = (u <= 1e-6);
end

function [p, q, Lk] = createModel(f, n, m, T, Lk, Vn, Vm)
    F = f(T);
    F1 = F + Lk;
    F2 = F - Lk;

    A = [Vn -F1'.*Vm -ones(numel(T),1);
        -Vn F2'.*Vm -ones(numel(T),1);
        zeros(numel(T), n) -Vm zeros(numel(T),1);
        zeros(numel(T), n) Vm zeros(numel(T),1);
        zeros(numel(T), n) zeros(numel(T),m) -ones(numel(T),1)];

    D = 0.9 * ones(numel(T), 1);
    D1 = 100 * ones(numel(T), 1);
    b = [F1'; -F2'; D; D1; ones(numel(T),1)];

    obj = [zeros(1, n + m), 1];
    options = optimoptions('linprog', 'Display', 'none', 'MaxIterations', 1e10);
    [w, ~, exitflag, output] = linprog(obj, A, b, [], [], [], [], options);

    if exitflag ~= 1 || isempty(w)
        error('Linear program failed: %s', output.message);
    end

    p = w(1:n);
    q = w(n+1:n+m);
    Lk = w(end);
end