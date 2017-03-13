n = 40; m = 5;
I = eye(n);
fndCounterEx = false;
i = 1;
tol = 1e-10;
while ~fndCounterEx
    V = randn(n,m);
    A = randn(n,n); B = randn(n,n);
    if rank(V) == m && rank(A) == n && rank(B) == n,
        V = orth(V);
        x = randn(m,1);
        x = V*x;
        a = randn; b = randn;
        v1 = (I - V*V') * (a * A + b * B) \ (A*x);
        v2 = (I - V*V') * (a * A + b * B) \ (B*x);
        if norm(v1 + (b/a) * v2) > tol,
            fndCounterEx = true;
            disp('CounterEx Found');
            i
            break;
        else
            i
            i = i+1;
        end
    else
        disp('rank deficient');
    end
end