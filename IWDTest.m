clc,clear,close all

t = 0:1e-4:30;

tA = sum(t <= 10);
tB = sum(t <= 20);

y = [chirp(t(1:tA),420,10,350) sin(550*2*pi*t(tA+1:tB-1)) chirp(t(tB:end)-t(tB),350,10,420)];
y1 = [chirp(t(1:tA),420,10,350) 0*sin(550*2*pi*t(tA+1:tB-1)) chirp(t(tB:end)-t(tB),350,10,420)];
y2 = [0*chirp(t(1:tA),420,10,350) sin(550*2*pi*t(tA+1:tB-1)) 0*chirp(t(tB:end)-t(tB),350,10,420)];

MoDAL.PlotTSWT(t,y,0,800)
hold on

m1 = (350-420)/(10);
PolyHigh1 = [m1*t(1:tA)+470 400*ones(1,length(t(tA+1:tB-1))) -m1*t(1:tA)+400]+50;
PolyLow1 = [m1*t(1:tA)+470 400*ones(1,length(t(tA+1:tB-1))) -m1*t(1:tA)+400]-100;

PolyHigh2 = [m1*t(1:tA)+470 400*ones(1,length(t(tA+1:tB-1))) -m1*t(1:tA)+400]+250;
PolyLow2 = [m1*t(1:tA)+470 400*ones(1,length(t(tA+1:tB-1))) -m1*t(1:tA)+400]+80;

plot(t,PolyHigh1,'r',t,PolyLow1,'r')
plot(t,PolyHigh2,'g',t,PolyLow2,'g')

Mode1 = MoDAL.IWDRange(t,y,PolyLow1,PolyHigh1);
Mode2 = MoDAL.IWDRange(t,y,PolyLow2,PolyHigh2);


figure
subplot(2,1,1)
plot(t,y1,'k',t,Mode1,'r--')
legend('Exact','IWD')

subplot(2,1,2)
plot(t,y2,'k',t,Mode2,'r--')
legend('Exact','IWD')

%%