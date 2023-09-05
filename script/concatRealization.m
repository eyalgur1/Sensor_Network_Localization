%%
data1=data;
%%
sk=11;
realizations=data1.realizations;
for l=sk:sk+39
    realizations{end+1}=data.realizations{l};
end
data1.realizations=realizations;
%%
sum_norm2_diff=zeros(1001,1);
for i=1:1001
    for l=1:50
        sum_norm2_diff(i)=sum_norm2_diff(i)+data1.realizations{l}.ADMM.norm2_diff(i);
    end
end
%%
%EML only
sum_norm2_diff=0;
for l=1:50
    sum_norm2_diff=sum_norm2_diff+data1.realizations{l}.EML.norm2_diff;
end

%%
data1.sum_norm2_diff.EML=sum_norm2_diff;
%%
data=data1;