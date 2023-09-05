%% INFORMATION

% DESCRIPTION: implements the AMFC method (AMU with one cluster).

% INPUTS:
% * Network - the true network structure
% * x0 - starting point
% * u0 - starting point
% * STOP - takes {eps,k} eps or k that holds the following arguments:
%          * 0<eps<1 - used in the stopping_criteria function
%          * k>=1 - maximum number of iteration used in the stopping_criteria
%                   function

function out=alg_AMAG(Network,x0,u0,STOP,Lipschitz,r,s)
%% INITIALIZATION

% setting matrices, objective function, etc.
[Qtilde,QDtilde,Atilde,ADtilde,Btilde,~,Ptilde,~,Xreal,~,...
    ~,~,~,K,m,N,n,Ne,~,~,...
    amatrix,F,RMSE,~,~,v_x_coef,v_const]=general_init(Network);

[stopping,RMSE_v,fun_v,t,tp]=stopping_criteria(F,STOP); % setting the stopping criterion
norm2_diff=fun_v;

if length(x0)>n*N % takes only the first n*(N-m) coordinates of x0
    x0 = x0(1:n*N);
end

x0=reshape(x0,n,N);
u0=reshape(u0,n,Ne);
x=x0;
u=u0;
uk=u;
xk=x;
yk=reshape(x,n*N,1);
tk=1;

A=kron(Atilde,eye(n));
B=kron(Btilde,eye(n));
AD=kron(ADtilde,eye(n));
QD=kron(QDtilde,eye(n));
P2=kron(2*Ptilde,eye(n));
u=reshape(u,n*size(u,2),1);
a=reshape(amatrix,n*m,1);
aBA=-2*a'*B'*A;
QDAD=-2*[QD;AD];
gF=@(x,u)(P2*x+(aBA+u'*QDAD)'); % gradient of F

if contains(Lipschitz,'c')
    L=2*K;
elseif contains(Lipschitz,'l')
    L=2*(2*max(diag((Qtilde')*Qtilde))+max(sum(Atilde)));
elseif contains(Lipschitz,'t')
    L=max(eigs(P2));
end

iter=0;
out_iter=0;

%% ALGORITHM

first_iteration=true;
while all(stopping(x, u, xk, uk, iter))==0 || first_iteration
    first_iteration=false;
    
    starta=tic;
    out_iter=out_iter+1;
    x=xk;
    u=uk;
    tk_prev=tk;
    
    x=reshape(x,n*N,1);
    uk=reshape(uk,n*Ne,1);
    
    starti=tic;
    max_inner_iter=s+2^floor(out_iter/r);
    
    inner_iter=0;
    while inner_iter<=max_inner_iter
        inner_iter=inner_iter+1;
        iter=iter+1;
        x=reshape(xk,n*N,1);
        xk=yk-(1/L)*gF(yk,uk);
        tk=(1+sqrt(1+4*tk^2))/2;
        yk=xk+((tk_prev-1)/tk)*(xk-x);
        if inner_iter<=max_inner_iter
            t(iter)=toc(starta); % upper bound on the non-parallelizable running time
            tp(iter)=t(end)/N; % running time estimation
            xk=reshape(xk,n,N);
            uk=reshape(uk,n,Ne);
            RMSE_v(iter)=RMSE(xk);
            fun_v(iter)=F(xk,uk);
            norm2_diff(iter)=sum(sum((xk(:,1:N)-Xreal(:,1:N)).^2));
            xk=reshape(xk,n*N,1);
            uk=reshape(uk,n*Ne,1);
        end
    end
    endi=toc(starti);
    
    xk=reshape(xk,n,N);
    x=reshape(x,n,N);
    uk=reshape(uk,n,Ne);
    
    [uk,~,endb]=u_update(xk,uk,v_x_coef,v_const,Ne,n,K,[]); % update the u block
    %iter=iter+1;
    
    t(iter)=toc(starta); % upper bound on the non-parallelizable running time
    tp(iter)=endi/N+endb; % running time estimation
    RMSE_v(iter)=RMSE(xk);
    fun_v(iter)=F(xk,uk);
    norm2_diff(iter)=sum(sum((xk(:,1:N)-Xreal(:,1:N)).^2));
    if ~mod(iter,100) || iter==1
        fprintf('AAG:    Iter=%2d    RMSE=%10.10f    F(k)-F(k+1)=%10.10f\n', iter, RMSE_v(iter), F(x,u)-F(xk,uk))
    end
end

%% OUTPUT
out=output(Network,xk,t,tp,iter,x0,u0,RMSE_v,fun_v,norm2_diff,1);
disp('AMFC is Done!')