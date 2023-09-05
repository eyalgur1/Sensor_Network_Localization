A=(net.Matrices.true_distances).*(net.Matrices.true_distances>0);
g = digraph(A);
bins = conncomp(g, 'Type', 'weak');
isConnected = all(bins == 1)
figure()
plot(g)
