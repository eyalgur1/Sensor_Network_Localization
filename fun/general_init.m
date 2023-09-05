function [Qtilde,QDtilde,Atilde,ADtilde,Btilde,BDtilde,Ptilde,a,Xreal,xreal,...
    IMnet,IMDnet,Laplaciannet,K,m,N,n,Ne,Neq,invPtilde,...
    amatrix,F,RMSE,x_const,x_u_coef,v_x_coef,v_const]=general_init(Network)
%% INFORMATION

% DESCRIPTION: returns matrices and functions which are relavent for all
%              algorithms.

% INPUTS:
% * Network - the true network structure

%% OUTPUT

Xreal=Network.Matrices.X_real; % output location matrix
xreal=Network.Matrices.x_real; % output location vector
K=Network.K; % total number of sensors
m=Network.GI.num_anchors; % total number of anchors
N=K-m; % total number of non-anchors
n=size(Network.Matrices.X_real,1); % dimension
Ne=size(Network.Matrices.Q_tilde,1)+size(Network.Matrices.A_tilde,1); % total number of non-anchor edges
Neq=size(Network.Matrices.Q_tilde,1); % total number of non-anchor to non-anchor edges
a=Network.Matrices.a; % anchors real location matrix
Ptilde=Network.Matrices.P_tilde;
invPtilde=inv(Ptilde);
Atilde=Network.Matrices.A_tilde;
Btilde=Network.Matrices.B_tilde;
Qtilde=Network.Matrices.Q_tilde;
QDtilde=Network.Matrices.QD_tilde;
ADtilde=Network.Matrices.AD_tilde;
BDtilde=Network.Matrices.BD_tilde;
IMnet=Network.Matrices.IM_net;
IMDnet=Network.Matrices.IMD_net;
Laplaciannet=Network.Matrices.Laplacian_net;
amatrix=reshape(a,n,m);
x_coef=(amatrix*Btilde'*Atilde);
v_x_coef=[QDtilde;ADtilde]'; % [QD;AD] for the objective function
x_const=x_coef*invPtilde; % (invP*A'*B)*a' for the objective function
x_u_coef=v_x_coef'*invPtilde; % invP*[QD', AD'] for the objective function
v_const=[zeros(2,size(QDtilde,1)),-amatrix*BDtilde']; % [zeros(size(QD,1),1);-BD*a] for the objective function
F=@(x,u)(sum(sum(((x*Ptilde').*x)))-sum(sum(2*u.*(x*v_x_coef+v_const)))-sum(sum(2*x_coef.*x))); % objective function
RMSE=@(x)(sqrt(1/K))*sqrt(sum(sum(([x,amatrix]-Xreal).^2))); % RMSE function

end