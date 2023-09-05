function out=alg_AMRC(Network,x0,u0,STOP,q,varargin)
%% INFORMATION

% DESCRIPTION: implements AMRC (AMU with random clusters), with/without warm start.

% INPUTS:
% * Network - the true network structure
% * x0 - starting point
% * u0 - starting point
% * STOP - a cell array {eps,k}, {eps} or {k} that holds the following arguments:
%          * 0<eps<1 - used in the stopping_criteria function
%          * k>=1 - maximum number of iteration used in the stopping_criteria
%               function
% * q - number of required random clusters (q<N).
% * varargin - an optional argument for a warm start:
%              * empty - no warm start
%              * SF or AG - SF or AG warm start
%              * k - k>=1 is the maximum number of warm start iterations. If
%                    k is set, then k-l is the maximum number of AMCC iterations
%              * Lipschitz - a string variable that defines the Lipschitz
%                            constant of AG. See AG algorithm.

%% INITIALIZATION

% setting matrices, objective function, etc.
[~,QDtilde,Atilde,ADtilde,Btilde,~,Ptilde,~,Xreal,~,...
    ~,~,~,~,~,N,n,Ne,Neq,~,...
    amatrix,F,RMSE,~,~,v_x_coef,v_const]=general_init(Network);
edges=Network.Matrices.edges;

% setting the stopping criterion
[stopping,RMSE_v,fun_v,t,tp]=stopping_criteria(F,STOP);
norm2_diff=fun_v;

% setting the starting points
if length(x0)>n*N % takes only the first n*N coordinates of x0
    x0=x0(1:n*N);
end
x0=reshape(x0,n,N); % sets x0 in matrix form
u0=reshape(u0,n,Ne); % sets u0 in matrix form

% generating q clusters (each sensor is a cluster)
clusters_strct=MakeClusters(Network,{q,'r'});
clusters=clusters_strct.clusters;

% setting coefficients according to the clusters
for i=1:q
    S_i=clusters{i}; % cluster S_i
    clusters{i}.Ptilde=Ptilde(S_i.nodes,S_i.nodes);
    clusters{i}.invPtilde=inv(clusters{i}.Ptilde);
    invPtilde=clusters{i}.invPtilde;
    clusters{i}.x_const=amatrix(:,S_i.anchor_neigbhors)*(invPtilde*Atilde(S_i.anchor_edges,S_i.nodes)'*Btilde(S_i.anchor_edges,S_i.anchor_neigbhors))'; % linear coeeficient of x_i: invP*A'*B*a
    clusters{i}.x_u_coef=(invPtilde*[QDtilde(S_i.nonanchor_edges,S_i.nodes)',ADtilde(S_i.anchor_edges,S_i.nodes)'])'; % quadratic coefficient of x_i.u_i: invP*[QD' AD']
    clusters{i}.x_x_coef=-(invPtilde*Ptilde(S_i.nonanchor_neigbhors,S_i.nodes)')'; % quadratic coefficient of x_i.x_i: -invP
end

%% ALGORITHM

%----------Warm Start----------%
if ~isempty(varargin) % with warm start
    warm_start=varargin;
    if ~isempty(strfind(warm_start{1},'SF')) % Portugal warm start
        outWarm=alg_SF(Network,x0,warm_start{2},u0);
    elseif ~isempty(strfind(warm_start{1},'AG')) % FISTA warm start
        outWarm=alg_AG(Network,x0,u0,warm_start{2},warm_start{3});
    end
    
    fun_v(1:warm_start{2})=outWarm.function_values(2:end); % update function values without the F(x0,u0)
    norm2_diff(1:warm_start{2})=outWarm.norm2_diff(2:end);
    RMSE_v(1:warm_start{2})=outWarm.RMSE(2:end); % update RMSE values without RMSE(x0)
    t(1:warm_start{2})=outWarm.time; % update time entries
    tp(1:warm_start{2})=outWarm.estimated_parallel_time; % update estimated parallel time entries
    x1=reshape(outWarm.location_estimation(:,1:N),n,N); % new starting point for AM-E
    [u1,~,~]=u_update(x1,u0,v_x_coef,v_const,Ne,n,K,[]); % new starting point for AM-E
    fun_v(warm_start{2})=F(x1,u1); % update function value at new starting point
    norm2_diff(warm_start{2})=sum(sum((x1(:,1:N)-Xreal(:,1:N)).^2));
    RMSE_v(warm_start{2})=RMSE(x1); % update RMSE value at new starting point
    x=x1;
    u=u1;
    xk=x;
    uk=u;
    iter=warm_start{2};
else % no warm start
    warm_start=false;
    x=x0;
    u=u0;
    xk=x;
    uk=u;
    iter=0;
end
%--------Warm Start End--------%

first_iteration=true;
while all(stopping(x,u,xk,uk,iter))==0 || first_iteration
    first_iteration=false;
    starta=tic;
    x=xk;
    u=uk;
    
    % x update and parallel run time calculations
    endi=zeros(1,q); % setting run time of each cluster
    for i=1:q % update the x block
        starti=tic; % run time for each cluster
        S_i = clusters{i};
        rel_edges=[S_i.nonanchor_edges;S_i.anchor_edges+Neq]; % cluster edges
        if ~isempty(S_i.nonanchor_neigbhors) % if the cluster has non-anchor neighbors
            xk(:,S_i.nodes)=xk(:,S_i.nonanchor_neigbhors)*S_i.x_x_coef+uk(:,rel_edges)*S_i.x_u_coef+S_i.x_const; % update x of each cluster
        else % if the cluster has only anchor neighbors
            xk(:,S_i.nodes)=uk(:,rel_edges)*S_i.x_u_coef+S_i.x_const; % update x of each cluster
        end
        endi(i)=toc(starti); % run time of the cluster (sensor)
    end
    
    [uk,u_time,~]=u_update(xk,uk,v_x_coef,v_const,Ne,n,K,edges); % update u
    u_cluster_time=zeros(1,q); % u update run time for each cluster
    for i=1:q
        u_cluster_time(i)=sum(u_time(clusters{i}.nodes)); % u calculation time for each sensor in the clusters
    end
    u_max_time=max(u_cluster_time); % run time of cluster with greatest u calculation time
    
    iter=iter+1;
    
    t(iter)=toc(starta); % upper bound on the total non-parallelizable running time
    RMSE_v(iter)=RMSE(xk);
    fun_v(iter)=F(xk,uk);
    norm2_diff(iter)=sum(sum((xk(:,1:N)-Xreal(:,1:N)).^2));
    tp(iter)=sum(endi)+u_max_time; % parallel run time
    
    % print values every 50 iterations
    if ~mod(iter,50) || iter==1
        fprintf('AMRC with %2d clusters:    Iter=%2d    RMSE=%10.10f    F(k)-F(k+1)=%10.10f\n',q,iter,RMSE_v(iter),F(x,u)-F(xk,uk))
    end
    
end

%% OUTPUT

if ~isempty(varargin) % with warm start
    if ~isempty(strfind(warm_start{1},'SF'))
        warm_start{end+1}={};
    end
    warm_start{end+1}=x1;
    out=output(Network,xk,t,tp,iter,x0,u0,RMSE_v,fun_v,norm2_diff,0,clusters_strct,warm_start);
else % no warm start
    out=output(Network,xk,t,tp,iter,x0,u0,RMSE_v,fun_v,norm2_diff,0,clusters_strct);
end

disp('AMRC is Done!')