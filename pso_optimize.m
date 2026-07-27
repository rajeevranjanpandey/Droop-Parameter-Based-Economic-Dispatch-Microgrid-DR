function [gbest, gbest_val, history] = pso_optimize(obj, dim, lb, ub, varargin)
%PSO_OPTIMIZE  Standard velocity-based Particle Swarm Optimization.
%
%   Same calling convention as tfo_optimize.m, so both can be benchmarked
%   on an identical objective/bounds/evaluation-budget for fair comparison.
%
%   [gbest, gbest_val, history] = pso_optimize(obj, dim, lb, ub, ...)
%
%   INPUTS
%     obj  : function handle, cost = obj(x), x is a 1 x dim row vector
%     dim  : number of decision variables
%     lb,ub: scalars or 1 x dim vectors (lower/upper bounds)
%
%   NAME-VALUE OPTIONS
%     'PopSize'       (50)
%     'Iters'         (2000)
%     'C1'            (2)      -- cognitive coefficient
%     'C2'            (2)      -- social coefficient
%     'Wmax'          (0.9)    -- inertia weight, start
%     'Wmin'          (0.4)    -- inertia weight, end (linear decrease)
%     'Seed'          ([])
%     'RecordHistory' (false)
%     'Verbose'       (false)
%
%   OUTPUTS
%     gbest     : 1 x dim best solution found
%     gbest_val : scalar best objective value
%     history   : Iters x 1 vector of gbest_val per iteration (if requested)

p = inputParser;
addParameter(p, 'PopSize', 50);
addParameter(p, 'Iters', 2000);
addParameter(p, 'C1', 2);
addParameter(p, 'C2', 2);
addParameter(p, 'Wmax', 0.9);
addParameter(p, 'Wmin', 0.4);
addParameter(p, 'Seed', []);
addParameter(p, 'RecordHistory', false);
addParameter(p, 'Verbose', false);
parse(p, varargin{:});
o = p.Results;

if ~isempty(o.Seed)
    rng(o.Seed);
end

lb = reshape(lb, 1, []);
ub = reshape(ub, 1, []);
if isscalar(lb), lb = repmat(lb, 1, dim); end
if isscalar(ub), ub = repmat(ub, 1, dim); end

% ---- initialize swarm (bsxfun for pre-R2016b compatibility) ----
X = bsxfun(@plus, lb, bsxfun(@times, rand(o.PopSize, dim), (ub - lb)));
V = 0.1 * bsxfun(@times, rand(o.PopSize, dim), (ub - lb));

fitness = zeros(o.PopSize, 1);
if o.Verbose
    fprintf('PSO: evaluating initial swarm (%d particles)...\n', o.PopSize);
end
for i = 1:o.PopSize
    fitness(i) = obj(X(i, :));
    if o.Verbose
        fprintf('  init particle %d/%d -> cost = %g\n', i, o.PopSize, fitness(i));
    end
end

Pbest = X;
Pbest_val = fitness;
[gbest_val, idx] = min(fitness);
gbest = X(idx, :);

if o.RecordHistory
    history = zeros(o.Iters, 1);
else
    history = [];
end

for t = 1:o.Iters
    W = o.Wmax - (t * (o.Wmax - o.Wmin) / o.Iters);

    r1 = rand(o.PopSize, dim);
    r2 = rand(o.PopSize, dim);

    V = W .* V ...
        + o.C1 .* r1 .* (Pbest - X) ...
        + o.C2 .* r2 .* (bsxfun(@minus, gbest, X));

    X = X + V;
    X = bsxfun(@min, bsxfun(@max, X, lb), ub);   % clamp to bounds

    for i = 1:o.PopSize
        f = obj(X(i, :));
        fitness(i) = f;
        if f < Pbest_val(i)
            Pbest_val(i) = f;
            Pbest(i, :) = X(i, :);
        end
    end

    [best_now, idx] = min(Pbest_val);
    if best_now < gbest_val
        gbest_val = best_now;
        gbest = Pbest(idx, :);
    end

    if o.RecordHistory
        history(t) = gbest_val;
    end

    if o.Verbose
        fprintf('PSO iter %d/%d -> best = %g\n', t, o.Iters, gbest_val);
    end
end

end
