% This script tests solve_CG.m for correctness
% Uses rmg.m

fprintf('Test 1: Identity matrices ...');
for n = 1:10
    A = eye(n);
    x = rand(n, 1);
    b = A * x;
    assert(norm(solve_CG(A, b) - x) < 1e-8);
end
fprintf('\tpassed\n');

fprintf('Test 2: Random matrices ...')
for n = 1:5
    A = rand(n);
    A = A' * A;
    while cond(A) > 1e5
        A = rand(n);
        A = A' * A;
    end
    x = rand(n, 1);
    b = A * x;
    assert(norm(solve_CG(A, b) - x) < 1e-8);
end
fprintf('\tpassed\n');

fprintf('Test 3: Random complex matrices ...')
for n = 1:5
    [~, A] = rmg(rand(n));
    x = rand(n, 1);
    b = A * x;
    assert(norm(solve_CG(A, b) - x) < 1e-8);
end
fprintf('\tpassed\n');

fprintf('All tests passed!\n');