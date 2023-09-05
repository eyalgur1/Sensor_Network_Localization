function [stopping,RMSE_v,fun_v,t,tp]=stopping_criteria(F,STOP)
%% INFORMATION

% DESCRIPTION: returns a stopping criteirion according to STOP.

% INPUTS:
% * F - objective function
% * STOP - a scalar k or cell(1,2)={k,l}. If k>=1 than k is max_iter. If 
%          k<1 than k is epsilon.
%          Accroding to k,l the function returns one of the following 
%          three stopping criterias:
%          1. (abs(F(x,u)-F(xk,uk)) < epsilon) || (iter >= max_iter)
%          2. abs(F(x,u)-F(xk,uk)) < epsilon
%          3. iter >= max_iter

%% SETTING STOPPING CRITERION
if length(STOP)==2
    if STOP{1}<1
        epsilon=STOP{1};
        max_iter=STOP{2};
    else
        epsilon=STOP{2};
        max_iter=STOP{1};
    end
    stopping=@(x,u,xk,uk,iter)((abs(F(x,u)-F(xk,uk))<epsilon) || (iter>=max_iter));
    RMSE_v=zeros(1,max_iter); % preallocation of RMSE values array
    fun_v=zeros(1,max_iter); % preallocation of function values values array
    t=zeros(1,max_iter); % preallocation of time entries array
    tp=zeros(1,max_iter); % preallocation of estimated parallel time entries array
elseif length(STOP)==1
    if STOP<1
        epsilon=STOP;
        stopping=@(x,u,xk,uk,iter)(abs(F(x,u)-F(xk,uk))<epsilon);
    else
        max_iter=STOP;
        stopping=@(x,u,xk,uk,iter)(iter>=max_iter);
        RMSE_v=zeros(1,max_iter); % preallocation of RMSE values array
        fun_v=zeros(1,max_iter); % preallocation of function values values array
        t=zeros(1,max_iter); % preallocation of time entries array
        tp=zeros(1,max_iter); % preallocation of estimated parallel time entries array
    end
end

end