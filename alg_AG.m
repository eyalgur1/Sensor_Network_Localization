function out=alg_AG(Network,x,u,max_iter,Lipschitz)
%% INFORMATION

% DESCRIPTION: implements AG as a warm start method for the function
%              x->F(x,u) for a fixed u.

% INPUTS:
% * Network - the true network structure
% * x - starting point
% * u - starting point for objective function values calculations
% * max_iter - maximum number of AG iterations
% * Lipschitz - constant method:
%               - 'c' - conservative 2*K (distributed)
%               - 'l' - liberal 2*(2*max_nonanchor_nonanchor_degree+max_nonanchor_anchor_degree) (distributed)
%               - 't' - tight 2*||P|| (centralized)

%% INITIALIZATION

% setting matrices, objective function, etc.
[Qtilde,QDtilde,Atilde,ADtilde,Btilde,~,Ptilde,a,Xreal,~,...
    ~,~,~,K,~,N,n,~,~,~,...
    ~,F,RMSE,~,~,~,~]=general_init(Network);

A=kron(Atilde,eye(n));
B=kron(Btilde,eye(n));
AD=kron(ADtilde,eye(n));
QD=kron(QDtilde,eye(n));
P2=kron(2*Ptilde,eye(n));
u_matrix=u;
u0_matrix=u_matrix;
u=reshape(u,n*size(u,2),1);

startgrad=tic;
x_const_grad=(-2*(a'*B'*A+u'*[QD;AD]))';
gF=@(x)(P2*x+x_const_grad); % gradient of F
endgrad=toc(startgrad);

% setting the Lipschitz constant of gF
if contains(Lipschitz,'c')
    L=2*K;
elseif contains(Lipschitz,'l')
    L=2*(2*max(diag((Qtilde')*Qtilde))+max(sum(Atilde)));
elseif contains(Lipschitz,'t')
    L=max(eigs(P2));
end

y=reshape(x,size(x,2)*n,1);
s=1; % stepsize initialization
x=y;
fun_v=zeros(1,max_iter); % function values
RMSE_v=zeros(1,max_iter); % RMSE values
norm2_diff=zeros(1,max_iter); % PRMSE values
t=zeros(1,max_iter); % non-parallel time values
iter=0;

%% ALGORITHM

while iter<max_iter
    tic
    x_prev=x;
    s_prev=s;
    
    % AG update rule
    x=y-(1/L)*gF(y);
    s=(1+sqrt(1+4*(s^2)))/2;
    y=x+((s_prev-1)/s)*(x-x_prev);
    
    iter=iter+1;
    
    % iteration information
    t(iter)=toc;
    x_matrix=reshape(x,n,N);
    x_prev_matrix=reshape(x_prev,n,N);
    RMSE_v(iter)=RMSE(x_matrix);
    norm2_diff(iter)=sum(sum((x_matrix(:,1:N)-Xreal(:,1:N)).^2));
    fun_v(iter)=F(x_matrix,u_matrix);
    %if ~mod(iter,100) || iter==1
     %   fprintf('AG warm start:    Iter=%2d    RMSE=%10.10f    F(k)-F(k+1)=%10.10f\n', iter, RMSE_v(iter), F(x_prev_matrix,u0_matrix)-F(x_matrix,u0_matrix))
    %end
end

%% OUTPUT

t(1)=t(1)+endgrad;
tp=t/N; % estimated parallel time computation
fun_v=fun_v-sum((sum(abs(QDtilde),2)/2).^2)-sum(sum(ADtilde,2).^2)-norm(kron(Btilde,eye(n))*a)^2;
out=output(Network,x,t,tp,iter,x_matrix,u_matrix,RMSE_v,fun_v,norm2_diff,0);
out.is_warm_start=true; % whether the method was used as a warm start method