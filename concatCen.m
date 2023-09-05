L=50;
m=10;
N=30;
radius=0.235:(0.55-0.235)/9:0.55;

strm=['m',num2str(m)];
strN=['net',num2str(N)];
for r=1:length(radius)
    R=radius(r);
    strR=['R',strrep(num2str(R),'.','')];
    data1.(strN).(strm).(strR).SF=data.(strN).(strm).(strR).SF;
end
    