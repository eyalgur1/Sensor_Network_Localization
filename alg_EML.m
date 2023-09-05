function out = alg_EML(net)
%This is my Maximum Likelihood ESDP formulation. Here I use Gaussian noise

%PP: 2xn matrix representing all point on 2D
%n : the number of total points
%m : the number of anchor points; the first m columns of PP
%r : the radio range
%nf: noisy factor
%degree: the limit on the number of edges connected to any sensor point, typically set it to 5-10

%Output
%PP: 2xn matrix representing all point on 2D, the last m columns are anchors
%dd: the uppertriangle distance matrix from PP(:i) to PP(:,j), i<j.
%Xiter: position estimation for all n-m unknown points made by the method

%%

disp('Generating the input file for SeDuMi...')

%--------2020--------%
m=net.anchors;
K=net.K;
N=K-m;
r=net.GI.R;
degree=100; % maximum degree of neighbors for non-anchors (set it high enough so it has no effect, since maximum degree is dictated by the radius)
Xreal=net.Matrices.X_real;
amatrix=reshape(net.Matrices.a,2,m);
PP=zeros(size(Xreal));
PP(:,1:m)=Xreal(:,K-m+1:end); % GenerateDistances takes PP with anchors as the first m columns
PP(:,m+1:end)=Xreal(:,1:K-m);
noised_distances=sparse(triu(net.Matrices.noised_distances)); % GenerateDistances takes dd (measured distances) as an upper trinagular matrix
%--------------------%

%Use Generate distances to generate the distance
[PP, ~, dnoise, Adjacency] = GenerateDistances(PP,m,K,r,0,degree,noised_distances); % set noise to zero so as not to create a different realization

starta=tic;
[i1,i2]=find(Adjacency(1:end-m,1:end-m));
Edges = size(i1,1);
[i3,i4]=find(Adjacency(1:end-m,end-m+1:end));
EdgesA = size(i3,1);

%%

%Define the variables for the SDP solution
%Position:
X = sdpvar(2,K-m, 'full','real');

%Noise
Dij = sdpvar(Edges,1, 'full','real');
dij = sdpvar(Edges,1, 'full','real');
Eik = sdpvar(EdgesA,1, 'full','real');
eik = sdpvar(EdgesA,1, 'full','real');

%Y=X'X as vec(Y) for nonzero element
Yoffdiag = sdpvar(Edges,1,'full','real');
Ydiag = sdpvar(K-m,1,'full','real');

%%
%Set the problem up
constraintset = (Dij(:)>=0) + (Eik(:)>=0);

%Inter-range measurements
ee=1:Edges;
ii = i1(ee); jj = i2(ee); ij=ee;
costfunction = sum(Dij)-2*dij'*diag(dnoise(ii,jj));
M0 = Ydiag(ii)+Ydiag(jj) - 2*Yoffdiag(ij) - Dij(ij);
constraintset = constraintset+ (M0(:)==0);
for ee=1:Edges
    ii = i1(ee); jj = i2(ee); ij=ee;
    Xij = X(1:2, [ii,jj]);
    %Yij = [Ydiag(ii), Yoffdiag(ij); Yoffdiag(ij), Ydiag(jj)];
    constraintset = constraintset+...
        ([eye(2), Xij; Xij', [Ydiag(ii), Yoffdiag(ij); Yoffdiag(ij), Ydiag(jj)]]>=0);
    constraintset = constraintset+...
        ([1, dij(ij); dij(ij), Dij(ij)]>=0);
end

%Anchors
ee=1:EdgesA;
ii = i3(ee); kk = K-m+i4(ee); ik=ee;
costfunction = costfunction+sum(Eik)-2*eik'*diag(dnoise(ii,kk));
M = Ydiag(ii)+diag(PP(1:2,kk)'*PP(1:2,kk)) - diag(2*X(1:2,ii)'*PP(1:2,kk)) - Eik(ik);
constraintset = constraintset+(M(:)==0);
for ee=1:EdgesA
    ik=ee;
    constraintset = constraintset+...
        ([1, eik(ik); eik(ik), Eik(ik)]>=0);
end

%%

%Solution via SeDuMi
disp('Solving using SeDuMi...')

options = sdpsettings('solver','sedumi','sedumi.eps',10^-5);

optimize(constraintset, costfunction,options);
%Check if constraints are OK [Optional]
%checkset(constraintset)

enda=toc(starta);
%%

%Post-processing
%Getting the position variable
Xiter = double(X);

%Some plotting
disp('Plotting...')
myPostProcessingGraphs(PP, Xiter, m);

%% OUTPUT
%--------2020--------%
out=struct;
out.location_estimation=[Xiter amatrix];
out.norm2_diff=sum(sum((Xiter(:,1:N)-Xreal(:,1:N)).^2)); % final norm2_diff of the algorithm
out.time=enda;
out.execution_date=date;
%--------------------%
end