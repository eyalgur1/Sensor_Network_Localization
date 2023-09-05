% Copyright (C) 2016  Nicola Piovesan (npiovesan@cttc.es)
%                     Tomaso Erseghe (erseghe@dei.unipd.it)
%  
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function out = alg_ADMM(net,iter,c,c2,tauc,nomefile,st)
%--------2020--------%
[~,~,~,~,~,~,~,~,Xreal,~,...
    ~,~,~,K,~,N,n,Ne,~,~,...
    ~,F,~,~,~,v_x_coef,v_const]=general_init(net);
u0=zeros(Ne*n,1);
u0=reshape(u0,n,Ne);
uk=u0;
fun_v=zeros(1,iter);
norm2_diff=zeros(1,iter);
%--------------------%

%-----PARAMETERS-----%
SolverMaxIter = 3;
SolverTolFun = 1e-4;
dc = 1.01;
tc = 0.98;
lmax = 1e3;
%--------------------%

sv.simpar.iterations=iter;
sv.simpar.c=c;
sv.simpar.c2=c2;
sv.simpar.SolverMaxIter=SolverMaxIter;
sv.simpar.SolverTolFun=SolverTolFun;
%disp(['[NETinfo] N=' num2str(net.N) ', anchors= ' num2str(net.anchors)]);
%disp(sv.simpar)
%disp(nomefile)

% initial point
a_pos = find(net.is_anchor); % set known anchor nodes positions
if (~exist('st'))
    st = [];
end
if isempty(st)
    startx = (0.01*randn(net.N,2)*[1;1i]).';    
    startx(a_pos) = net.positions(net.is_anchor);
else
    %startx = st;
    %--------2020--------%
    if length(st)>n*N % takes only the first n*(N-m) coordinates of st
        st = st(1:n*N);
    end
    st=reshape(st,n,length(st)/n);
    startx=st(1,:)+st(2,:)*1i;
    startx(a_pos) = net.positions(net.is_anchor);
    %--------------------%
end

%--------2020--------%
startx_matrix=[real(startx);imag(startx)];
%--------------------%

% useful for local update function 
anv = zeros(1,net.K);
anv(a_pos) = net.positions(net.is_anchor);

% pre-allocate timing info
T = zeros(iter,net.K);

%Set fmincon/fminunc
options = optimset('GradObj','on','Display', 'off',...
    'TolFun',SolverTolFun,'MaxIter',SolverMaxIter,'Hessian','user-supplied');
                
for t=1:iter %Main loop (iterations)
%     fprintf(['IT: ' num2str(t) ' ~']);

    %--------2020--------%
    if ~mod(t,50) || t==1
        fprintf(['ADMM-H:    Iter=' num2str(t)]);
    end
    %--------------------%
    
    if t==1
        % Initialize local positions x
        for i = 1:net.K
            tic
            node = net.node{i}; % current node
            dt(i).Ni = node.N_neigh; % number of neigbors
            
            tmp = startx([i,node.neighbors]).';
            tmp = [real(tmp),imag(tmp)].';            
            dt(i).x = tmp(:); % vector with positions of current node and its neighbors
            dt(i).apos = anv(node.neighbors);
            dt(i).Ai = [ones(dt(i).Ni,1) -eye(dt(i).Ni); ones(dt(i).Ni,1) eye(dt(i).Ni)];
            dt(i).Ai = kron(dt(i).Ai,eye(2));
            
            if isinf(tauc)
                dt(i).lim = -inf; % start with non-convex
                dt(i).ci = c2; % cost value
            else
                dt(i).lim = 0; % start with convex
                dt(i).ci = c; % cost value
            end
            dt(i).id = [i node.neighbors]; % ID of each elements
            dt(i).r = node.distance.'; % ranging distances
            dt(i).is_anch = ismember(dt(i).id, a_pos); % is_anchor
            dt(i).neigh_pos = [0 node.neigh_pos]; % position of this node on neighbor list

            dt(i).m = dt(i).Ai*dt(i).x; % build messages
            dt(i).lambda = zeros(4*dt(i).Ni,1); % reset memory
            fin = startx; % starting guess

            T(t,i) = T(t,i)+toc;
        end
    else
        for i = 1:net.K
            tic
            % prepare constants (complex form)  
            y = (dt(i).Ai'*(dt(i).z-dt(i).lambda./dt(i).ci));
            y = (y(1:2:end)+1i*y(2:2:end))/2;
            y(1) = y(1)/dt(i).Ni;
            y([1==0, dt(i).apos~=0]) = dt(i).apos(dt(i).apos~=0);
            c_lo = ((dt(i).apos'==0)+2*dt(i).ci)*dt(i).Ni;
            r_lo = dt(i).r;
            l_lo = dt(i).lim;

            % update local positions            
            if net.node{i}.is_anchor==0
                xii = dt(i).x(1:2);
                xii = fminunc(@(x)FiMOD(x,r_lo,c_lo,y,l_lo),xii,options);
                dt(i).x(1:2) = xii;
                fin(i) = xii(1) + 1i*xii(2);
            end

            % update other positions            
            tmp = fin(i)-y(2:end);
            tmp2 = (1-r_lo./abs(tmp))/(1+2*dt(i).ci);
            tmp = y(2:end) + tmp.*tmp2.*(tmp2>l_lo);
            tmp(dt(i).apos~=0) = dt(i).apos(dt(i).apos~=0);            
            dt(i).x(3:2:end)=real(tmp);
            dt(i).x(4:2:end)=imag(tmp);

            % build messages
            dt(i).m = dt(i).Ai*dt(i).x+dt(i).lambda./dt(i).ci; 
            T(t,i) = T(t,i)+toc;
        end
        oldTi = Ti;
    end
    sv.x = fin;
    sv.RMSE(t) = sqrt((norm(net.positions-fin).^2)./length(fin));
    %--------2020--------%
     xk_matrix=[real(sv.x);-imag(sv.x)];
     uk=u_update(xk_matrix(:,1:N),uk,v_x_coef,v_const,Ne,n,K,[]);
     fun_v(t)=F(xk_matrix(:,1:N),uk);
     norm2_diff(t)=sum(sum((xk_matrix(:,1:N)-Xreal(:,1:N)).^2));
    %--------------------%
    
    %% Evaluate the mixing values
    for i=1:net.K
        tic
        dd=length(dt(i).m);
        z1=0.5.*(dt(i).m(1:dd/2));
        z2=0.5.*(dt(i).m(dd/2+1:end));
        for j=1:dt(i).Ni
            idn=dt(i).id(j+1);
            np=dt(i).neigh_pos(j+1);
            z1([2*j-1,2*j])=z1([2*j-1,2*j])-0.5*dt(idn).m([2*np-1,2*np]);
            z2([2*j-1,2*j])=z2([2*j-1,2*j])+0.5*dt(idn).m([2*(np+dt(idn).Ni)-1,2*(np+dt(idn).Ni)]);
        end
        dt(i).z = [z1;z2];
        % also evaluate the primal gap
        Ti(i) = max(abs(dt(i).Ai*dt(i).x-dt(i).z));
        T(t,i) = T(t,i)+toc;
    end
    
    
    %% Update memory
    if t>1
        c_vec = cat(1,dt.ci);
        for i=1:net.K
            tic
            % memory update
            tmp = dt(i).lambda+dt(i).ci.*(dt(i).Ai*dt(i).x-dt(i).z);
            tmp(tmp>lmax) = lmax;
            tmp(tmp<-lmax) = -lmax;
            dt(i).lambda = tmp;            
            if dt(i).lim~=0 
                % c update
                tmp = 1 + (dc-1)*(Ti(i)>tc*oldTi(i));
                dt(i).ci = max(c_vec(dt(i).id)) * tmp;
            else
                % transition
                if(t>5)&&((Ti(i)<tauc)||(max(c_vec(dt(i).id))>c))
                    dt(i).ci = c2;
                    dt(i).lim = -inf;
                end
            end
            T(t,i) = T(t,i)+toc;
        end
    end
    
    sv.perc(t) = 100*sum(cat(1,dt.lim)~=0)/net.K;
        
    
    %% Save and Plot informations
%     fprintf('%.2f sec - (~: %.3f, M: %.3f, P: %3.0f, G: %.4f)\n',...
%     sum(T(t,:)), mean(T(t,1:end-10)),max(T(t,1:end-10)),...
%     sv.perc(t),sv.RMSE(t));

    %--------2020--------%
    if ~mod(t,50) || t==1
        fprintf('    RMSE=%.10f\n' ,sv.RMSE(t));
    end
    %--------------------%

    assignin('base','T',T);
    if (mod(t,10)==0)
        save(nomefile);
    end
end

%--------2020--------%
out = output(net,[real(sv.x);-imag(sv.x)],sum(T,2),max(T,[],2),iter,startx_matrix(:,1:N),u0,sv.RMSE,fun_v,norm2_diff,0);
ADMM_parameters=struct;
ADMM_parameters.eps=c;
ADMM_parameters.zeta=c2;
ADMM_parameters.tau=tauc;
out.ADMM_parameters=ADMM_parameters;
out.sv=sv;
disp('ADMM is done!')
%--------------------%
end

%%
%%%%%%%%%%%%%%%%%%%% FiMOD function %%%%%%%%%%%%%%%%%%%%

% Copyright (C) 2016  Nicola Piovesan (npiovesan@cttc.es)
%                     Tomaso Erseghe (erseghe@dei.unipd.it)
%  
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%

function [f, g, H] = FiMOD(x,r,dv,y,lim)

xii = x(1)+1i*x(2);

q1 = xii-y(2:end); % e' la q
q3 = abs(q1); % e' l'abs
q2 = q3-r;
q1 = q1./q3;
q4 = q1.*(q2>lim).*sqrt(r./q3./dv); % ora e' q/|q|
q2 = q2.*(q2>lim); % e' quella con ^*

f = 0.5.*abs(xii-y(1)).^2+0.5*sum(q2.^2./dv);

g = xii-y(1) + sum((q2./dv).*q1);
g = [real(g) ; imag(g)];

re = real(q4);
im = imag(q4);
c = sum(re.*im);
H = eye(2)*(1+sum(q2./dv./q3)) + [sum(re.^2), c;c, sum(im.^2)];

f(isinf(f))=0;
g(isinf(g))=0;
H(isinf(H))=0;
f(isnan(f))=0;
g(isnan(g))=0;
H(isnan(H))=0;
end