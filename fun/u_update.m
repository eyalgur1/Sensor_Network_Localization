function [uk,times,tot_time]=u_update(xk,uk,v_x_coef,v_const,Ne,n,K,edges)
%% INFORMATION

% DESCRIPTION: updates the u block in AM-T, AM-E and AM-G. Additionally,
%              gives a parellel computing time estimation.

% INPUTS:
% * xk - the x block after its the k-th iteration.
% * uk - the u block before its k-th iteration.
% * v_x_coef - [QD;AD]
% * v_const - zeros(size(QD,1),1);-BD*a]
% * Ne - total number of non-anchor edges
% * n - dimension
% * IMnet - incidence matrix of the entire network for AMFD, AMCC and AMU,
%           and [] otherwise.

%% INITIALIZATION

v = xk*v_x_coef+v_const; % [QD*xk; (AD*xk - BD*a)];
times=zeros(K,1); % size(IMnet,2) is the total number of sensors K
tot_time=0;

%% U UPDATE AND PARALLEL TIME COMPUTATION

for j=1:Ne
    timea=tic;
    norm_v_j=norm(v(:,j));
    if norm_v_j>0
        uk(:,j)=v(:,j)/norm_v_j;
    else
        uk(:,j)=zeros(n,1);
    end
    
    time=toc(timea);
    
    if ~isempty(edges)
        times(edges(j,2))=times(edges(j,2))+time;
        times(edges(j,3))=times(edges(j,3))+time;
    end
    tot_time=tot_time+time;
end

end