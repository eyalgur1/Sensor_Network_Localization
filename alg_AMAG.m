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
iter=0;
out_iter=0;

%% ALGORITHM

first_iteration=true;
while all(stopping(x, u, xk, uk, iter))==0 || first_iteration
    first_iteration=false;
    out_iter=out_iter+1;
    
    start_keep=tic;
    x=xk;
    u=uk;
    end_keep=toc(start_keep);
    
    AGiter=s+2^(floor(out_iter/r));
    outAG=alg_AG(Network,xk,uk,AGiter,Lipschitz);
    xk=outAG.location_estimation(:,1:N);
    
    fun_v(iter+1:iter+AGiter)=outAG.function_values(2:end);
    RMSE_v(iter+1:iter+AGiter)=outAG.RMSE(2:end);
    norm2_diff(iter+1:iter+AGiter)=outAG.norm2_diff(2:end);
    t(iter+1:iter+AGiter)=outAG.time;
    tp(iter+1:iter+AGiter)=outAG.estimated_parallel_time;
    
    [uk,u_time,~]=u_update(xk,uk,v_x_coef,v_const,Ne,n,K,[]);
    u_max_time=max(u_time(1:N));
    
    fun_v(end)=F(xk,uk);
    RMSE_v(end)=RMSE(xk);
    norm2_diff(end)=sum(sum((xk(:,1:N)-Xreal(:,1:N)).^2));
    tp(end)=tp(end)+u_max_time+end_keep/N;
    t(end)=t(end)+sum(u_time(1:N))+end_keep;
    
    iter=iter+AGiter;
    
    %if ~mod(iter,100) || iter==1
        fprintf('AMAG:    Iter=%2d    RMSE=%10.10f    F(k)-F(k+1)=%10.10f\n', iter, RMSE_v(iter), F(x,u)-F(xk,uk))
    %end
end

%% OUTPUT
out=output(Network,xk,t,tp,iter,x0,u0,RMSE_v,fun_v,norm2_diff,1);
disp('AMAG is Done!')