%% Remove unnecessary fields from the output structure to decrease file size

L=50;
m=10;
N=30;
radius_ticks=10;
radius=0.235:(0.55-0.235)/(radius_ticks-1):0.55;

strm=['m',num2str(m)];
strN=['net',num2str(N)];
remove_fields_from_output={'QD_tilde','AD_tilde','BD_tilde','IMD_net','dd_noise','noised_distances'};
for r=1:length(radius)
    R=radius(r);
    strR=['R',strrep(num2str(R),'.','')];
    net_structure=data.(strN).(strm).(strR).nets{1};
    for i=1:L
        dd_noise=sparse(data.(strN).(strm).(strR).nets{i}.Matrices.dd_noise);
        noised_distances=sparse(data.(strN).(strm).(strR).nets{i}.Matrices.noised_distances);        
        data.(strN).(strm).(strR).nets{i}=[];
        data.(strN).(strm).(strR).nets{i}.dd_noise=dd_noise;
        data.(strN).(strm).(strR).nets{i}.noised_distances=noised_distances;
    end
    data.(strN).(strm).(strR).net_structure=net_structure;
    data.(strN).(strm).(strR).net_structure.Matrices=rmfield(data.(strN).(strm).(strR).net_structure.Matrices,remove_fields_from_output);
end

%save(['output/cenExpAnc_',num2str(N),'_',num2str(m),datestr(now,'yyyy_mm_dd_HH_MM'),'.mat'],'data')