% plot_fields={'ADMM','SF','PG','PG,SF','PG,AG','AMFC',...
%     'AMCC','AMCC,SF','AMCC,AG','AMFD','AMFD,SF','AMFD,AG',...
%     'AMGC','AMGC,SF','AMGC,AG'};

%% Plot PRMSE
figure(1)
nr=10;
N=980;
m=20;
R=0.061;
sigma=0.061*0.007;
avg_deg=11.0940;
max_iter=2000;

plot_fields={'AMAG','AMFC','AMCC_100AGl'};
num_methods=length(plot_fields);

for j=1:num_methods
    if contains(plot_fields{j},'AMAG')
        plotAMAG=sqrt(data.sum_norm2_diff.(plot_fields{j})/nr)/N;
        semilogy(plotAMAG(1:max_iter),'linewidth',1)
        hold on
    else
        semilogy(sqrt(data.sum_norm2_diff.(plot_fields{j})/nr)/N,'linewidth',1)
        hold on
    end
end
num_iter=length(data.sum_norm2_diff.(plot_fields{j}));
legend(strrep(plot_fields,'_',' '),'Interpreter','latex','FontSize',12)
xlim([0,num_iter])
xlabel('$\textrm{Iterations}$','Interpreter','latex','fontsize',14)
ylabel('$\frac{\textrm{PRMSE}}{N}$','Interpreter','latex','fontsize',14)
title({['$K=\ $',num2str(N+m),', $N=\ $',num2str(N),', $m=\ $',num2str(m),' random network in $\left[-\frac{1}{2},\frac{1}{2}\right]^2$'],...
    ['$R=\ $',num2str(R),', $L=\ $',num2str(nr),', $\sigma=\ $',num2str(sigma),...
    ', $\mathrm{Average\ degree}=\ $',num2str(avg_deg)]},'Interpreter','latex','fontsize',14)


%% Calculate function value at real location

R=data.original_network.net.GI.R;
PP=data.original_network.net.Matrices.X_real;
funv_real=zeros(nr,1);
for l=1:nr
    dd=data.realizations{l}.net.noised_distances;
    net_l=create_realization(m,PP,dd,0,R);
    [edges_nn,edges_na]=create_edges(net_l);
    edges_na(:,1)=edges_na(:,1)+length(edges_nn);
    edges=[edges_nn;edges_na]; % edges with real distances
    for e=1:size(edges,1) % calculate the real function value of the realization
        edge=edges(e,:);
        funv_real(l)=funv_real(l)+(norm(PP(:,edge(2))-PP(:,edge(3)))-edge(4))^2;
    end
end

%% Plot function values

figure(2)
nr=10;
N=980;
m=20;
R=0.061;
sigma=0.061*0.007;
avg_deg=11.0940;
max_iter=2000;

plot_fields={'AMAG','AMFC','AMCC_100AGl'};
num_methods=length(plot_fields);

for j=1:num_methods
    funv_avg=zeros(num_iter,1);
    for i=1:num_iter
        for l=1:nr
            funv_avg(i)=funv_avg(i)+data.realizations{l}.(plot_fields{j}).function_values(i)/nr;
        end
    end
    semilogy(funv_avg(1:max_iter),'linewidth',1)
    hold on
end

semilogy((sum(funv_real)/nr)*ones(num_iter,1),'color','black','linewidth',1)

num_iter=length(data.sum_norm2_diff.(plot_fields{j}));
legend([strrep(plot_fields,'_',' '),'Real location'],'Interpreter','latex','FontSize',12)
xlim([0,num_iter])
xlabel('$\textrm{Iterations}$','Interpreter','latex','fontsize',14)
ylabel('$\textrm{Average Function Value}$','Interpreter','latex','fontsize',14)
title({['$K=\ $',num2str(N+m),', $N=\ $',num2str(N),', $m=\ $',num2str(m),' random network in $\left[-\frac{1}{2},\frac{1}{2}\right]^2$'],...
    ['$R=\ $',num2str(R),', $L=\ $',num2str(nr),', $\sigma=\ $',num2str(sigma),...
    ', $\mathrm{Average\ degree}=\ $',num2str(avg_deg)]},'Interpreter','latex','fontsize',14)