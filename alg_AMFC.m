function out=alg_AMFC(Network,x0,u0,STOP)
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

%% INITIALIZATION

% setting matrices, objective function, etc.
[~,~,~,~,~,~,~,~,Xreal,~,...
    ~,~,~,K,~,N,n,Ne,~,~,...
    ~,F,RMSE,x_const,x_u_coef,v_x_coef,v_const]=general_init(Network);

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

%% ALGORITHM

first_iteration=true;
while all(stopping(x, u, xk, uk, iter))==0 || first_iteration  
    first_iteration=false;
    
    starta=tic;
    x=xk;
    u=uk;
    
    starti=tic;
    xk=x_const+u*x_u_coef; % invP*([QD', AD']*u + A'*B*a) - update the x block   
    endi=toc(starti);
 
    [uk,~,endb]=u_update(xk,uk,v_x_coef,v_const,Ne,n,K,[]); % update the u block
    iter=iter+1;
    
    t(iter)=toc(starta); % upper bound on the non-parallelizable running time
    tp(iter)=endi+endb; % running time estimation
    RMSE_v(iter)=RMSE(xk);
    fun_v(iter)=F(xk,uk);
    norm2_diff(iter)=sum(sum((xk(:,1:N)-Xreal(:,1:N)).^2));
    if ~mod(iter,100) || iter==1
        fprintf('AMFC:    Iter=%2d    RMSE=%10.10f    F(k)-F(k+1)=%10.10f\n', iter, RMSE_v(iter), F(x,u)-F(xk,uk))  
    end
end

%% OUTPUT
out=output(Network,xk,t,tp,iter,x0,u0,RMSE_v,fun_v,norm2_diff,1);
disp('AMFC is Done!')