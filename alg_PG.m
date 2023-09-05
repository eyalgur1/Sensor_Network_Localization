function out=alg_PG(Network,x0,max_iter,varargin)
%% INFORMATION

% DESCRIPTION: implements the PG method as described in the paper https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7032222

% INPUTS:
% * Network - the true network structure
% * x0 - starting point
% * max_iter - maximum number of iterations
% * varargin - an optional argument for a warm start:
%              * empty - no warm start
%              * SF or AG - SF or AG warm start
%              * k - k>=1 is the maximum number of warm start iterations. If
%                    k is set, then k-l is the maximum number of AMCC iterations
%              * Lipschitz - a string variable that defines the Lipschitz
%                            constant of AG. See AG algorithm.


%% INITIALIZATION

[Qtilde,QDtilde,Atilde,ADtilde,~,~,~,~,Xreal,~,...
    ~,~,~,K,m,N,n,Ne,~,~,...
    amatrix,F,RMSE,~,~,v_x_coef,v_const]=general_init(Network);

if length(x0)>n*N % takes only the first n*(N-m) coordinates of x0
    x0 = x0(1:n*N);
end
x0=reshape(x0,n,N);
xa_init=[x0 amatrix]; % starting point in case of no warm start
delta=diag(Qtilde'*Qtilde); % non-anchor to non-anchor degree
delta_max=max(delta); % maximum non-anchor to non-anchor degree
u=zeros(2,Ne); % for function value update
distances=[sum(abs(QDtilde),2)/2;sum(ADtilde,2)]; % measured (noisy) distances
L=2*delta_max+max(sum(Atilde))+2; % upper bound on the Lipschitz constant set to 2*(2*max_nonanchor_nonanchor_degree+max_nonanchor_anchor_degree)+2
b=(L-delta(1:K-m)-sum(Atilde)')/L; % vector of constants for PG updates
fun_v=zeros(1,max_iter); % function values
RMSE_v=zeros(1,max_iter); % RMSE values
norm2_diff=zeros(1,max_iter); % PRMSE values
t=zeros(1,max_iter); % non-parallel time values
tp=zeros(1,max_iter); % parallel time estimations
u0=zeros(Ne*n,1); % auxiliary variable of the method
u0=reshape(u0,n,Ne);
uk=u0;

[edges_Q,edges_A]=create_edges(Network);
edges_A(:,1)=edges_A(:,1)+length(edges_Q);
edges=[edges_Q(:,2:3);edges_A(:,2:3)]; % indices of the edges

edges_matrix=zeros(K,K); %if (i,j) is the k-th edge, then the (i,j) entry is k if i<j and -k if j<i
for i=1:Ne
    edges_matrix(edges(i,1),edges(i,2))=i;
    edges_matrix(edges(i,2),edges(i,1))=-i;
end

%% ALGORITHM

%----------Warm Start----------%
if ~isempty(varargin)
    warm_start=varargin;
    if ~isempty(strfind(warm_start{1},'SF'))
        u0 = reshape(warm_start{4},Ne*n,1);
        outWarm=alg_SF(Network,x0,warm_start{2},u0);
        u0 = reshape(u0,n,Ne);
    elseif ~isempty(strfind(warm_start{1},'AG'))
        u0 = reshape(warm_start{4},n,Ne);
        x0 = reshape(x0,n,N);
        outWarm=alg_AG(Network,x0,u0,warm_start{2},warm_start{3});
    end
    
    fun_v(1:warm_start{2})=outWarm.function_values(2:end); % update function values without the F(x0,u0)
    norm2_diff(1:warm_start{2})=outWarm.norm2_diff(2:end);
    RMSE_v(1:warm_start{2})=outWarm.RMSE(2:end); % update RMSE values without RMSE(x0)
    t(1:warm_start{2})=outWarm.time; % update time entries
    tp(1:warm_start{2})=outWarm.estimated_parallel_time; % update estimated parallel time entries
    x1=reshape(outWarm.location_estimation(:,1:N),n,N); % new starting point
    [u1,~,~]=u_update(x1,u0,v_x_coef,v_const,Ne,n,K,[]); % new starting point for PG
    fun_v(warm_start{2})=F(x1,u1); % update function value at new starting point
    norm2_diff(warm_start{2})=sum(sum((x1(:,1:N)-Xreal(:,1:N)).^2));
    RMSE_v(warm_start{2})=RMSE(x1); % update RMSE value at new starting point
    xa_init=[x1 amatrix]; % starting point of the algorithm in case of a warm start
end
%--------Warm Start End--------%

for i=1:Ne % initializtion of the auxiliary variable for function value updates
    ind=edges(i,:);
    v=xa_init(:,ind(1))-xa_init(:,ind(2));
    if norm(v)>0
        u(:,i)=distances(i)*(v/norm(v));
    else
        u(:,i)=distances(i)*[1;0];
    end
end

if ~isempty(varargin) % set the iteration counter in case of warm start or not
    iter=warm_start{2};
else
    iter=0;
end
x=reshape(xa_init,n,K); % set the starting point in case of warm start or not

while iter<max_iter
    x_prev=x;
    u_prev=u;
    starta=tic;
    for i=1:K-m % xi update
        starti=tic;
        endi=zeros(K-m,1);
        s=zeros(n,1);
        for j=Network.node{1,i}.neighbors
            edge=abs(edges_matrix(i,j));
            s=s+sign(edges_matrix(i,j))*u_prev(:,edge)+x_prev(:,j);
        end
        x(:,i)=b(i)*x_prev(:,i)+(1/L)*s; % xi update
        
        for j=Network.node{1,i}.neighbors % auxiliary variable update
            edge=abs(edges_matrix(i,j));
            v=((L-1)/L)*u_prev(:,edge)+(1/L)*sign(edges_matrix(i,j))*(x_prev(:,i)-x_prev(:,j));
            if norm(v)>0
                u(:,edge)=distances(edge)*(v/norm(v));
            else
                u(:,edge)=distances(edge)*[1;0];
            end
        end
        endi(i)=toc(starti);
    end
    iter=iter+1;
    t(iter)=toc(starta);
    tp(iter)=max(endi);
    RMSE_v(iter)=RMSE(x(:,1:N));
    
    uk_prev=uk; % for function value updates
    [uk,~,~]=u_update(x(:,1:N),uk,v_x_coef,v_const,Ne,n,K,[]);
    fun_v(iter)=F(x(:,1:N),uk);
    norm2_diff(iter)=sum(sum((x(:,1:N)-Xreal(:,1:N)).^2));
    
    if ~mod(iter,50) || iter==1
        fprintf('PG:    Iter=%2d    RMSE=%10.10f    F(k)-F(k+1)=%10.10f\n', iter, RMSE_v(iter), F(x_prev(:,1:N),uk_prev)-F(x(:,1:N),uk))
    end
end

%% OUTPUT

if ~isempty(varargin) % with warm start
    if ~isempty(strfind(warm_start{1},'SF'))
        warm_start{end+1}={};
    end
    warm_start{end+1}=x1;
    out=output(Network,x(:,1:N),t,tp,iter,x0,u0,RMSE_v,fun_v,norm2_diff,0,warm_start);
else % no warm start
    out=output(Network,x(:,1:N),t,tp,iter,x0,u0,RMSE_v,fun_v,norm2_diff,0);
end

disp('PG is Done!')