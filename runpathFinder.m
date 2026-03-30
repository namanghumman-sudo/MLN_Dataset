function d = runpathFinder(sessionData, knn, parameters, s0)
% Minimal placeholder for debugging only
% Returns unweighted shortest-path distance from seed s0 on knn graph

n = size(sessionData, 1);

if s0 < 1 || s0 > n
    error('Seed s0 out of range');
end

A = sparse(knn);
A = max(A, A');     % make undirected
A = A ~= 0;         % unweighted adjacency

d = inf(n,1);
d(s0) = 0;

visited = false(n,1);
queue = zeros(n,1);
head = 1;
tail = 1;
queue(tail) = s0;
visited(s0) = true;

while head <= tail
    u = queue(head);
    head = head + 1;

    nbrs = find(A(u,:));
    for v = nbrs
        if ~visited(v)
            visited(v) = true;
            d(v) = d(u) + 1;
            tail = tail + 1;
            queue(tail) = v;
        end
    end
end

% replace disconnected nodes
bad = isinf(d);
if any(bad)
    mx = max(d(~bad));
    if isempty(mx)
        mx = 0;
    end
    d(bad) = mx + 1;
end
end