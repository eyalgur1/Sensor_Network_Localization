%%
sum_norm2_diff=zeros(1000,1);
for j=1:1000
    for i=1:10
        sum_norm2_diff(j)=sum_norm2_diff(j)+data.realizations{i}.SF.norm2_diff(j);
    end
end
plot(sqrt(sum_norm2_diff/10)/990)
hold off
figure()
funv_avg=zeros(1000,1);
for j=1:1000
    for i=1:10
        funv_avg(j)=funv_avg(j)+data.realizations{i}.SF.function_values(j)/10;
    end
end
plot(funv_avg)