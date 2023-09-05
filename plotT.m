nr=50;
max_iter=1001;

fields=fieldnames(data.realizations{end});
fields=fields(2:end); % methods in the data structure
num_methods=length(fields);

Xreal=data.original_network.net.Matrices.X_real;
[edges_nn,edges_na]=create_edges(data.original_network.net);
edges_na(:,1)=edges_na(:,1)+length(edges_nn);
edges=[edges_nn(:,1:3);edges_na(:,1:3)];

K=data.original_network.net.K;
m=data.original_network.net.anchors;
N=K-m;

funv_real_avg=0;
for l=1:nr
    nd=data.realizations{l}.net.noised_distances;
    for e=1:size(edges,1)
        edge=edges(e,:);
        funv_real_avg=funv_real_avg+((norm(Xreal(:,edge(2))-Xreal(:,edge(3)))-nd(edge(2),edge(3)))^2)/nr;
    end
end

ha=tight_subplot(2,1,[.01 .03],[.1 .1],[.1 .1]);

for j=1:num_methods
    if contains('EML',fields{j})
        fva=zeros(max_iter,1);
        for l=1:nr
            nd=data.realizations{l}.net.noised_distances;
            X_EML=data.realizations{l}.(fields{j}).location_estimation;
            for e=1:size(edges,1) % last function value of EML
                edge=edges(e,:);
                fva=fva+ones(max_iter,1)*(norm(X_EML(:,edge(2))-X_EML(:,edge(3)))-nd(edge(2),edge(3)))^2;
            end
        end
        snd=ones(max_iter,1)*(data.sum_norm2_diff.(fields{j}));
    else
        fva=zeros(max_iter,1);
        for l=1:nr
            fva=fva+((data.realizations{l}.(fields{j}).function_values)');
        end
        snd=data.sum_norm2_diff.(fields{j});
    end
    axes(ha(1))
    loglog(sqrt(snd/nr)/N)
    xticks([])
    ylabel('PRMSE/N')
    hold on
    axes(ha(2))
    loglog(fva/nr)
    xlim([1 max_iter])
    ylabel('Average Function Value')
    hold on   
end
xlabel('Iterations')
set(gcf, 'Position',  [500, 80, 750, 900])
hold off

%%
ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
function [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
% tight_subplot creates "subplot" axes with adjustable gaps and margins
%
% [ha, pos] = tight_subplot(Nh, Nw, gap, marg_h, marg_w)
%
%   in:  Nh      number of axes in hight (vertical direction)
%        Nw      number of axes in width (horizontaldirection)
%        gap     gaps between the axes in normalized units (0...1)
%                   or [gap_h gap_w] for different gaps in height and width
%        marg_h  margins in height in normalized units (0...1)
%                   or [lower upper] for different lower and upper margins
%        marg_w  margins in width in normalized units (0...1)
%                   or [left right] for different left and right margins
%
%  out:  ha     array of handles of the axes objects
%                   starting from upper left corner, going row-wise as in
%                   subplot
%        pos    positions of the axes objects
%
%  Example: ha = tight_subplot(3,2,[.01 .03],[.1 .01],[.01 .01])
%           for ii = 1:6; axes(ha(ii)); plot(randn(10,ii)); end
%           set(ha(1:4),'XTickLabel',''); set(ha,'YTickLabel','')
% Pekka Kumpulainen 21.5.2012   @tut.fi
% Tampere University of Technology / Automation Science and Engineering
if nargin<3; gap = .02; end
if nargin<4 || isempty(marg_h); marg_h = .05; end
if nargin<5; marg_w = .05; end
if numel(gap)==1;
    gap = [gap gap];
end
if numel(marg_w)==1;
    marg_w = [marg_w marg_w];
end
if numel(marg_h)==1;
    marg_h = [marg_h marg_h];
end
axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh;
axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
py = 1-marg_h(2)-axh;
% ha = zeros(Nh*Nw,1);
ii = 0;
for ih = 1:Nh
    px = marg_w(1);
    
    for ix = 1:Nw
        ii = ii+1;
        ha(ii) = axes('Units','normalized', ...
            'Position',[px py axw axh], ...
            'XTickLabel','', ...
            'YTickLabel','');
        px = px+axw+gap(2);
    end
    py = py-axh-gap(1);
end
if nargout > 1
    pos = get(ha,'Position');
end
ha = ha(:);
end