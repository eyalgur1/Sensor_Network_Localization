function out=alg_SF(Network,x0,max_iter,varargin)
%% INFORMATION
% DESCRIPTION: implements the SF algorithm either as a warm start or not,
% as described in the paper https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7153574

% INPUTS:
% * Network - the true network structure
% * x0 - starting point
% * max_iter - maximum number of FISTA iterations
% * varargin - an optinal argument when used as a warm start:
%              * u0 - u starting point (necessary for function values updates)

%% INITIALIZATION

% setting matrices, objective function, etc.
[~,QDtilde,~,ADtilde,Btilde,~,~,a,Xreal,~,...
    IMnet,IMDnet,~,K,m,N,n,Ne,~,~,...
    ~,F,RMSE,~,~,v_x_coef,v_const]=general_init(Network);

nonanchor_indices=1:N;
Laplacian=sparse(IMDnet'*IMDnet); % Laplacian matrix of the entire network with distances
QTQtilde=sparse(IMnet'*IMnet); % Laplacian matrix of the entire network 

degrees=[]; % generates the sensor degree diagonal matrix
for i=m+1:K
    degrees=[degrees, QTQtilde(i,i)];
end
deltamax=max(degrees);

E={}; 
for i=1:K % sets the sets of neighbors
    E{end+1}=find(QTQtilde(i,:)==-1); %finding all the neigbors in one go
end

sizes_of_anchor=[]; 
for i=nonanchor_indices % finds the maximal number of anchor neighbors for the non-anchors
    E_i=E{i}; 
    A_i=setdiff(E_i,nonanchor_indices);
    sizes_of_anchor=[sizes_of_anchor, size(A_i,2)];
end

L=2*deltamax+max(sizes_of_anchor); % Lipschizt constant (according to the paper)

if ~isempty(varargin) % checks if SF is used as a warm start
    is_warm_start=true;
    u0=reshape(varargin{1},n,Ne);
    uk=reshape(varargin{1},n,Ne);
else
    is_warm_start=false;
    u0=zeros(Ne*n,1);
    u0=reshape(u0,n,Ne);
    uk=u0;
end

if size(x0,1)==2 % the input is a matrix of non anchors if used as warm start
    x0_matrix=x0;
    x0=reshape(x0,n*N,1);
end
if length(x0)>n*N % takes only the first n*(N-m) coordinates of x0
    x0=x0(1:n*N);
    x0_matrix=reshape(x0,n,N);
else
    x0_matrix=reshape(x0,n,N);
end
x0=[x0;a]; % the anchor locations are assumed to be exact
xk=x0;
xk_prev=xk;
fun_v=zeros(1,max_iter); % function values
RMSE_v=zeros(1,max_iter); % RMSE values
norm2_diff=zeros(1,max_iter); % PRMSE values
t=zeros(1,max_iter); % non-parallel time values
tp=zeros(1,max_iter); % parallel time estimations
iter=0;

%% Algorithm

while iter<max_iter
    starta=tic;
    w=xk+((iter-2)/(iter+1))*(xk-xk_prev); % Nesterov's optimal gradient method    
    grad=zeros(K*n,1);
    f_v=0;
    enda=toc(starta); % for parallel time
    
    endi=zeros(1,length(nonanchor_indices)); % time for each sensor i at each iteration    
    for i=nonanchor_indices % considering non-anchors only
        starti=tic; % for parallel time
        E_i=[E{i}];  % neighbor set of non-anchor i
        i_index=(i-1)*n+1:i*n;
        for j=E_i
            j_index=(j-1)*n+1:j*n;
            absLap=abs(Laplacian(i,j));
            sqrtabsLap=sqrt(absLap);
            normdiff=norm(w(i_index)-w(j_index));
            if normdiff>sqrtabsLap
                grad(i_index)=grad(i_index)+(w(i_index)-w(j_index))*(1-sqrtabsLap/normdiff); %and easy calculation of the gradient
                f_v=f_v+0.5*(normdiff-absLap);
            end
        end
        endi(i)=toc(starti); % for parallel time
    end
    
    % update the algorithm
    startb=tic; % for parallel time
    xk_prev=xk;
    xk=w-(1/L)*grad;    
    iter=iter+1;
    endb=toc(startb); % for parallel time
    
    t(iter)=toc(starta); % for non-parallel time
    tp(iter)=max(endi)+(enda+endb)/K; % for parallel time
    xk_matrix=reshape(xk,n,K); % for F evaluation
    xk_prev_matrix=reshape(xk_prev,n,K); % for F evaluation
    RMSE_v(iter)=RMSE(xk_matrix(:,1:N));
    norm2_diff(iter)=sum(sum((xk_matrix(:,1:N)-Xreal(:,1:N)).^2));
    
    if is_warm_start % updating function values and printing (warm start)
        fun_v(iter)=F(xk_matrix(:,1:N),uk); % no u update required
        if ~mod(iter,10) || iter==1
            fprintf('SF warm start:    Iter=%2d    RMSE=%10.10f    F(k)-F(k+1)=%10.10f\n', iter, RMSE_v(iter), F(xk_prev_matrix(:,1:N),u0)-F(xk_matrix(:,1:N),u0))
        end
    else % updating function values and printing (not a warm start)
        uk_prev=uk;
        [uk,~,~]=u_update(xk_matrix(:,1:N),uk,v_x_coef,v_const,Ne,n,K,[]);
        fun_v(iter)=F(xk_matrix(:,1:N),uk);
        if ~mod(iter,50) || iter==1
            fprintf('SF:    Iter=%2d    RMSE=%10.10f    F(k)-F(k+1)=%10.10f\n', iter, RMSE_v(iter), F(xk_prev_matrix(:,1:N),uk_prev)-F(xk_matrix(:,1:N),uk))
        end
    end
end

%% OUTPUT
if is_warm_start
    fun_v=fun_v-sum((sum(abs(QDtilde),2)/2).^2)-sum(sum(ADtilde,2).^2)-norm(kron(Btilde,eye(n))*a)^2;
end
out=output(Network,xk_matrix(:,1:N),t,tp,iter,x0_matrix,u0,RMSE_v,fun_v,norm2_diff,0);
out.is_warm_start=is_warm_start;
if ~is_warm_start
    disp('SF is done!')
end