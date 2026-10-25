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
Modes(:,:,1) = MoDAL.IFD_Mult(t',y(:,1),lf,uf);
Modes(:,:,2) = MoDAL.IFD_Mult(t',y(:,2),lf,uf);
Modes(:,:,3) = MoDAL.IFD_Mult(t',y(:,3),lf,uf);


charFreq{1} = [2;3.5]/(2*pi);
charFreq{2} = [2;3.5]/(2*pi);
charFreq{3} = [2;3.5]/(2*pi);


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