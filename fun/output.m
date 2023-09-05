function out = output(Network,xk,t,tp,iter,x0,u0,RMSE_v,fun_v,norm2_diff,cen,varargin)
%% INFORMATION
% Description: sets the output of the algorithms.

% INPUTS:
% * Network - the true network structure
% * xk - the output location estimation of the algorithm
% * t - array of times over the iterations
% * tp - array of estimated parallel time over iterationts
% * x0 - starting point of the algorithm or warm start algorithm
% * u0 - starting point of the algorithm
% * RMSE - RMSE function
% * RSMSE_v - array of RMSE values over the iterations
% * F - objective function
% * fun_v - array of function values values over the iterations
% * cen - if 1 the method is centralized and if 0 the method is distributable
% * varargin - an optional argument for clusters, random permutation and warm start:
%              * Clusters - struture of clusters used by AMU
%              * is_rand_perm - logical 1 if random permutation is used over the iterations in AMU
%              * warm_start - a cell for which:
%                             * warm_start{1} - 'FISTA' or 'Portugal'
%                             * warm_start{2} - number of warm start iterations
%                             * warm_start{3} - FISTA Lipschitz constant.
%                                               Set to {} if 'Portugal'.
%                             * warm_start{4} - warm start output.

%% INITIALIZATION

[~,QDtilde,~,ADtilde,Btilde,~,~,a,Xreal,~,...
    ~,~,~,~,~,N,n,~,~,~,...
    amatrix,F,RMSE,~,~,~,~]=general_init(Network);
xk=xk(1:n*N);
xk=reshape(xk,n,N);
X=[xk,amatrix];

%% OUTPUT

out=struct;
out.location_estimation=X; % output point
out.time=t; % this is an upper bound
out.total_time=sum(t); % this is an upper bound
if cen % only true for AMFC
    out.estimated_time=tp; % does not consider the updating of the previous iteration in the variables x and u 
    out.total_estimated_time=sum(tp);
else
    out.estimated_parallel_time=tp;
    out.total_estimated_parallel_time=sum(tp);
end
out.iterations=iter; % total number of iterations
out.RMSE=[RMSE(x0),RMSE_v];
out.norm2_diff=[sum(sum((x0(:,1:N)-Xreal(:,1:N)).^2)),norm2_diff];
out.function_values=[F(x0,u0),fun_v]+sum((sum(abs(QDtilde),2)/2).^2)+sum(sum(ADtilde,2).^2)+norm(kron(Btilde,eye(n))*a)^2; % adding the constant term sum(d_ij)^2+||Ba||^2 back to F;

%%%% Warm Start Section %%%%
if ~isempty(varargin) % warm start
    for i=1:length(varargin)
        var=class(varargin{i});
        switch var
            case 'struct'
                out.Clusters_info=varargin{i};
            case 'cell'
                % warm start
                if ~isempty(strfind(varargin{i}{1},'FISTA'))
                    out.warm_start='FISTA';
                    out.number_of_FISTA_iterations=varargin{i}{2};
                    out.FISTA_starting_point=x0;
                    out.FISTA_output_point=varargin{i}{4};
                    out.FISTA_Lipschitz_constant=varargin{i}{3};
                elseif ~isempty(strfind(varargin{i}{1},'ortugal'))
                    out.warm_start='Portugal';
                    out.Portugal_iterations=varargin{i}{2};
                    out.Portugal_starting_point=x0;
                    out.Portugal_output_point=varargin{i}{4};
                end
        end
    end
    %%%% End of Warm Start Section %%%%
else % no warm start
    out.starting_point=x0;
end

out.execution_date=date;