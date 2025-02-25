clc,clear,close all

t = (0:1e-3:100)';

y(:,1) = 0.75*sin(2*t)+0.35*sin(3.5*t);
y(:,2) = 0.8*sin(2*t)-0.25*sin(3.5*t);
y(:,3) = 0.6*sin(2*t)+0.4*sin(3.5*t);

% figure
% plot(t,y)

% % lf = {[0 2.75],[0 2]};
lf = [0 2.75]/(2*pi);
uf = [2.75 4]/(2*pi);

charFreq{1} = [2;3.5]/(2*pi);
charFreq{2} = [2;3.5]/(2*pi);
charFreq{3} = [2;3.5]/(2*pi);

Modes = MoDAL.IFD(t',y,lf,uf);
Mirror.On = 0;
% Mirror.FMirror{1} = [2;3.5];
% Mirror.FMirror{2} = [2;3.5];
% Mirror.FMirror{3} = [2;3.5];
% Mirror.Chop1{1} = [0.2;0.2];
% Mirror.Chop1{2} = [0.2;0.2];
% Mirror.Chop1{3} = [0.2;0.2];
% Mirror.Chop2{1} = [0.2;0.2];
% Mirror.Chop2{2} = [0.2;0.2];
% Mirror.Chop2{3} = [0.2;0.2];
% Mirror.Ini{1} = [0;1];
% Mirror.Ini{2} = [0;1];
% Mirror.Ini{3} = [0;1];
% Mirror.Fin{1} = [0;0];
% Mirror.Fin{2} = [1;0];
% Mirror.Fin{3} = [0;1];
IMFs = MoDAL.WBEMD(t,y,0,5,charFreq,"Mirror",Mirror);
size(Modes)
%%
figure
plot(t,Modes(:,1,1).*Modes(:,2,1),'k',t,IMFs{1}(:,1).*IMFs{1}(:,2),'r')

r =((IMFs{1}(:,1).*IMFs{1}(:,2))./y(:,1).^2);
r(1) = 0;
sum(r)