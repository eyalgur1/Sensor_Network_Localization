function out=alg_EML_CVX(net)
%% INFORMATION

% DESCRIPTION: implemnts using CVX the EML relaxatiom method (Problem (12)) described in https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6725647&tag=1 

% INPUTS:
% * net - network structure as obtained by createNetwork.

% OUTPUTS:
% * out - structure containing:
%         location_estimation: final obtained location of all sensors (anchors are last m columns).
%         time: total run time of the cvx progrma.

%% INITIALIZATION

n=size(net.Matrices.X_real,1); 
K=net.K; 
m=net.anchors; 
N=K-m;
QDtilde=net.Matrices.QD_tilde; 
ADtilde=net.Matrices.AD_tilde;
amatrix=reshape(net.Matrices.a,n,m); 
sigma=net.GI.Noise;
[edges_nn,edges_na]=create_edges(net);

%% ALGORITHM

starta=tic;
%cvx_solver sdpt3
cvx_begin sdp
variable Y(N,N) symmetric
variables X(n,N) delta(size(QDtilde,1),1) epsilon(size(ADtilde,1),1)  d(size(QDtilde,1),1) e(size(ADtilde,1),1)
first_term=0;
for i=1:size(QDtilde,1)
    edge=edges_nn(i,1:4);
    first_term=first_term+delta(edge(1))-2*d(edge(1))*edge(4)+edge(4)^2;
end
second_term=0;
for i=1:size(ADtilde,1)
    edge=edges_na(i,1:4);
    second_term=second_term+epsilon(edge(1))-2*e(edge(1))*edge(4)+edge(4)^2;
end
f=(1/sigma)*(first_term+second_term);
minimize f
subject to
for i=1:size(QDtilde,1)
    edge=edges_nn(i,1:3);
    Y(edge(2),edge(2))+Y(edge(3),edge(3))-2*Y(edge(2),edge(3))==delta(edge(1))
    [1 d(edge(1));d(edge(1)) delta(edge(1))]>=0
    delta(edge(1))>=0
    d(edge(1))>=0
end
for i=1:size(ADtilde,1)
    edge=edges_na(i,1:3);
    Y(edge(2),edge(2))-2*(X(:,edge(2))')*amatrix(:,edge(3)-N)+norm(amatrix(:,edge(3)-N))^2==epsilon(edge(1))
    [1 e(edge(1));e(edge(1)) epsilon(edge(1))]>=0
    epsilon(edge(1))>=0
    e(edge(1))>=0
end
for i=1:size(QDtilde,1)
    edge=edges_nn(i,1:3);
    [eye(n) X(:,edge(2)) X(:,edge(3));X(:,edge(2))' Y(edge(2),edge(2)) Y(edge(2),edge(3));X(:,edge(3))' Y(edge(2),edge(3)) Y(edge(3),edge(3))]>=0
end
cvx_end
enda=toc(starta);

%% OUTPUT

out=struct;
out.location_estimation=[X amatrix]; % location output, m last columns are anchors
out.time=enda;
out.execution_date=date;