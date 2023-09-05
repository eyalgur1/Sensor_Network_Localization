function CRLB=CRLB(net,sigma)
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
FIM=zeros(2*N); % FIM

%% IMPLEMENTATION

for i=1:N % update the FIM
    for j=net.node{i}.neighbors
        delta=Xreal(:,i)-Xreal(:,j);
        Temp=-(delta*delta')/(norm(delta)^2);
        if j<=N
            FIM(2*i-1:2*i,2*j-1:2*j)=Temp;
        end
        FIM(2*i-1:2*i,2*i-1:2*i)=FIM(2*i-1:2*i,2*i-1:2*i)-Temp;
    end
end
FIM=(FIM+FIM')/2;
FIM=FIM./(sigma^2); % FIM

% check if FIM is invertible
[~,flag] = chol(FIM);
if flag>0
    error('FIM is not PD')
end

invF=FIM\eye(2*N);
CRLB=trace(invF);