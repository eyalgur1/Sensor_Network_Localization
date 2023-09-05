function Network=createNetwork_bench(n,K,m,R,sigma,PP,dd)
%% INFORMATION

% DESCRIPTION: creates the benchmark sensor networks that are available at
%              https://web.stanford.edu/~yyye/Col.html.

% INPUTS:
% * n - dimension of the space
% * K - total number of sensors
% * m - number of anchor sensors
% * R - max radio range as suggested by the network's data (r argument). It
%       has no effect in truncating the network.
% * sigma - Gaussian(0, noise_D) noise level as suggested by the paper
%             https://ieeexplore.ieee.org/abstract/document/7782840/. It
%             has no effect in noising the network.
% * PP - real positions of the sensors, as obtained by the datasets
%        available at https://web.stanford.edu/~yyye/Col.html.
% * dd - Adjacency matrix of the networks with noised distances as obtained 
%        by the datasets available at 
%        https://web.stanford.edu/~yyye/Col.html.

% OUTPUTS:
% * Q_tilde - incidence matrix of non-anchors vs. non-anchors
% * QD_tilde - Q_tilde with noised distances.
% * A_tilde - indicator matrix for non-anchors in non-anchors vs. anchors edges
% * B_tilde - indicator matrix for anchors in non-anchors vs. anchors edges
% * BD_tilde - B_tilde with noised distances
% * a - vector of anchor real locations
% * X_real - real network location matrix
% * x_real - real network location vector
% * P_tilde - Q_tilde'*Q_tilde + A_tilde'*A_tilde
% * IM_net - incidence matrix of the entire network
% * Laplacian_net - Laplacian matrix of the entire network

%% CREATING THE NETWORK

S=PP'; % Location matrix, each row is a sensor 
distance=zeros(K);
for i=1:n
    distance=distance+(S(:,i)*ones(1,K)-ones(K,1)*S(:,i)').^2;
end
distance=sqrt(distance); % returns the distance

dist_in_range=distance.*(dd>0);

num_edges=nnz(dist_in_range)/2;

Network=struct;

Network.K=K;
Network.anchors=m;
if n==2
    Network.positions=(S(:,1)+S(:,2)*1i)';
else
    Network.positions=S';
end
Network.is_anchor=logical([zeros(1,K-m),ones(1,m)]);

IM=zeros(num_edges,K);
IMD=IM;
Nodes = cell(1,K);
cur_l=0;
for i=1:K    
    %create edges for node i
    i_neighbors=find(dist_in_range(i,i+1:K))+i;
    n_i_neighbors=length(i_neighbors);        
    %update the incidence matrix and the incidence matrix with distances
    IM(cur_l+1:cur_l+n_i_neighbors,i)=1;
    IM(cur_l+1:cur_l+n_i_neighbors,i_neighbors)=-eye(n_i_neighbors);    
    IMD(cur_l+1:cur_l+n_i_neighbors,i)=dist_in_range(i,i_neighbors)';
    IMD(cur_l+1:cur_l+n_i_neighbors,i_neighbors)=-diag(dist_in_range(i,i_neighbors)');
    cur_l=cur_l+n_i_neighbors;
    
    %create the node
    Nodes{i}.neighbors=[find(dist_in_range(i,1:i-1)),i_neighbors];
    Nodes{i}.N_neigh=length(Nodes{i}.neighbors);
    Nodes{i}.distance=abs(dist_in_range(i,Nodes{i}.neighbors));
    Nodes{i}.real_distance=abs(distance(i,Nodes{i}.neighbors));
    Nodes{i}.position=Network.positions(i);
    if i>K-m
        Nodes{i}.is_anchor=true;
    else
        Nodes{i}.is_anchor=false;
    end
end
for j=1:K
    for l=1:Nodes{j}.N_neigh      
        %updating the position of current node in neighbor's neighbor list
        neighbor_index=Nodes{j}.neighbors(l);
        Nodes{j}.neigh_pos(l)=find(Nodes{neighbor_index}.neighbors==j);
    end
end

%% OUTPUT

edges_with_anchor=(sum(abs(IM(:,K-m+1:K)),2)==1);
A_tilde=abs(IM(edges_with_anchor,1:K-m));
B_tilde=abs(IM(edges_with_anchor,K-m+1:K));
AD_tilde=abs(IMD(edges_with_anchor,1:K-m));
BD_tilde=abs(IMD(edges_with_anchor,K-m+1:K));
a=reshape(S(K-m+1:K,:)',m*n,1); %locations of anchors
X_real=S';
x_real=reshape(X_real,K*n,1);
edges_of_nona_nona=(sum(abs(IM(:,1:K-m)),2)==2);
Q_tilde=IM(edges_of_nona_nona,1:K-m);
QD_tilde=IMD(edges_of_nona_nona,1:K-m);
P_tilde=Q_tilde'*Q_tilde+A_tilde'*A_tilde;
Laplacian=IM'*IM;

GI = struct;
Matrices = struct;

Network.node = Nodes;
GI.num_anchors = m;
GI.num_non_anchors = K-m;
GI.num_rel_edges = sum(edges_with_anchor)+sum(edges_of_nona_nona);
GI.anchor_indices = K-m+1:K;
GI.R = R;
GI.Noise = sigma;
GI.avg_neigh = mean(diag(Laplacian));

Matrices.Q_tilde = sparse(Q_tilde);
Matrices.QD_tilde = sparse(QD_tilde);
Matrices.A_tilde = sparse(A_tilde); %indicator matrix - rows are edges between non anchors and anchors, columns are non-anchors nodes
Matrices.AD_tilde = sparse(AD_tilde);
Matrices.B_tilde = sparse(B_tilde); %indicator matrix - rows are edges between non anchors and anchors, columns are anchors nodes
Matrices.BD_tilde = sparse(BD_tilde);
Matrices.P_tilde = sparse(P_tilde); 
Matrices.a = a; 
Matrices.X_real = X_real;
Matrices.x_real = x_real;
Matrices.IM_net = sparse(IM);
Matrices.IMD_net = sparse(IMD);
Matrices.Laplacian_net = sparse(Laplacian);
Matrices.true_distances = sparse(dist_in_range);

Network.GI = GI;
Network.Matrices = Matrices;

end