Ts = 0.001;

t = 0:Ts:10;
L = length(t);
y = zeros(2,2,L);
y(1,1,:) = sin(t)+2;
y(1,2,:) = sin(t*5)*3+4;
y(2,1,:) = sin(t*20)+2;
y(2,2,:) = sin(t*10)+2;
%%
figure
treshold = 0.05;
[tp, tn, ip, in] = mat2inSpike.timeContrast(y,Ts,treshold);
plot(tp,ip,'.b');
hold on;
plot(tn,in,'.r');