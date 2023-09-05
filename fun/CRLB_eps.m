function CRLB=CRLB_eps(net,sigma,J,bias)
%% INFORMATION

% DESCRIPTION: calculates the CRLB of a network with a given Gaussian noise (sigma) according to the papers:
%              1. https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=1212671
%              2. https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=4753258 

% INPUTS:
% * net - network structure as obtained by createNetwork.

% OUTPUTS:
% * CRLB - value of the bound (take sqaure root if needed).

%% INITIALIZATION

Xreal=net.Matrices.X_real;
N=net.K-net.anchors;
n=size(Xreal,1);
FIM=zeros(n*N); % FIM

%% IMPLEMENTATION

for i=1:N % update the FIM
    for j=net.node{i}.neighbors
        delta=Xreal(:,i)-Xreal(:,j);
        Temp=-(delta*delta')/(norm(delta)^2);
        if j<=N
            FIM(n*(i-1)+1:n*i,n*(j-1)+1:n*j)=Temp;
        end
        FIM(n*(i-1)+1:n*i,n*(i-1)+1:n*i)=FIM(n*(i-1)+1:n*i,n*(i-1)+1:n*i)-Temp;
    end
end
FIM=(FIM+FIM')/2;
FIM=FIM./(sigma^2); % FIM

% check if FIM is invertible
[~,flag] = chol(FIM);
if flag>0
    error('FIM is not PD')
end

invF=FIM\eye(n*N);

%sqrt(trace(invF))/N
sqrt(trace(J*invF*J'+bias*bias'))/N
%sqrt(trace(J*invF*J'+bias*bias'))/N
%sqrt(trace(J*invF*J'-bias*bias'))/N

CRLB=trace(J*invF*J'+bias*bias');