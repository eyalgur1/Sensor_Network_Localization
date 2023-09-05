function out=alg_ESDPD(net)
%% INFORMATION

% DESCRIPTION: implements ESDPD method accoridng to the paper https://stanford.edu/~boyd/papers/pdf/semidef_relax_sensor_net_loc.pdf
%              and based on their code.

% INPUTS:
% * net - the true network structure

%% INITIALIZATION

PP=net.Matrices.X_real;
K=size(PP,2); % total number of sensors
n=size(net.Matrices.X_real,1);
m=net.anchors;
N=K-m;
a=net.Matrices.a;
amatrix=reshape(a,n,m);
noised_distances=net.Matrices.noised_distances;

%% ALGORITHM

% Input: 
%  * PP: 2xn matrix representing all point true locations on 2D
%        where n is the number of total sensor and anchor points
%  * m: the number of Anchor points; last m columns of PP; it should be at least 3
%  * noised_distances: the sparse (and noisy) distance matrix between PP(:,i) and PP(:,j)
% Output:
%  * Xiter: position estimation for all n-m unknown points made by the method

starta=tic;
Xiter = ESDPD(PP,m,noised_distances);
enda=toc(starta);

%% OUTPUT
out=struct;
out.location_estimation=[Xiter amatrix]; % location output, m last columns are anchors
out.time=enda;
out.total_estimated_parallel_time=enda/N; % this method is fully distributed and parallelizable
out.execution_date=date;

end
