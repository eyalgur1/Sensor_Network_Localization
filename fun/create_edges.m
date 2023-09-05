function [edges_nn,edges_na]=create_edges(net)
%% INFORMATION

% DESCRIPTION: creates edgee matrix for a network.

% INPUTS:
% * net - network structure as obtained by the functions createNetwork or create_realization.

% OUTPUTS:
% * edges_nn - row size is the total number of non-anchor to non-anchor edges, column size is 4.
%              Column 1: row (edge) index starting from 1.
%              Column 2: first non-anchor index.
%              Column 3: second non-anchor index.
%              Column 4: measured length of edge.
% * edges_na - row size is the total number of non-anchor to anchor edges, column size is 4.
%              Column 1: row (edge) index starting from 1.
%              Column 2: non-anchor index.
%              Column 3: nchor index.
%              Column 4: measured length of edge.

%% IMPLEMENTATION

Qtilde=net.Matrices.QD_tilde; % sensor-edge incidence matrix of non-anchors to non-anchors with distances
Atilde=net.Matrices.AD_tilde;
Btilde=net.Matrices.BD_tilde;
ABtilde=[Atilde Btilde]; % sesnor-edge incidence matrix of non-anchors to anchors with distances

edges_nn=zeros(size(Qtilde,1),4); % non-anchor to non-anchor edges
edges_na=zeros(size(ABtilde,1),4); % non-anchor to anchor edges

% set edges_nn
for i=1:size(Qtilde,1) % total number of non-anchor to non-anchor edges
    edge=find(Qtilde(i,:)); % indices of the two non-anchor sensors
    edges_nn(i,1)=i; % edge number (1 through size(Qtilde,1))
    edges_nn(i,2)=edge(1); % index of first sensor
    edges_nn(i,3)=edge(2); % index of second sensor
    edges_nn(i,4)=sum(abs(Qtilde(i,:)))/2; % length of edge
end

% set edges_na
for i=1:size(ABtilde,1) % total number of non-anchor to anchor edges
    edge=find(ABtilde(i,:));
    edges_na(i,1)=i; % edge number (1 through size(Atilde,1))
    edges_na(i,2)=edge(1); % index of first sensor (always non-anchor)
    edges_na(i,3)=edge(2); % index of second sensor (always an anchor)
    edges_na(i,4)=sum(Atilde(i,:)); % length of edge
end