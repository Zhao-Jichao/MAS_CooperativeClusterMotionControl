% 数值仿真

clear
clc

% 领航者参数
ul(:,1) = [5 2]';
vl(:,1) = [0 0]';

% 跟随者参数
% PXf1(:,1)= 3; PYf1(:,1)=-2; VXf1(:,1)= 0; VYf1(:,1)= 0;
% PXf2(:,1)= 6; PYf2(:,1)= 4; VXf2(:,1)= 0; VYf2(:,1)= 0;
% PXf3(:,1)=-5; PYf3(:,1)=-5; VXf3(:,1)= 0; VYf3(:,1)= 0;
% PXf4(:,1)=-4; PYf4(:,1)= 3; VXf4(:,1)= 0; VYf4(:,1)= 0;

PX = [3  6 -5 -4]';
PY = [-2 4 -5  3]';

VX = [0  0  0  0]';
VY = [0  0  0  0]';

UX = [0  0  0  0]';
UY = [0  0  0  0]';

% 期望距离：d
% 通信距离：R
d = 0.5;
R = 1.0;

alpha = 10;

% 时间参数
tbegin = 0;
tfinal = 60;
dT = 0.001;
T(:,1) = 0; 

% 计算次数
times = (tfinal - tbegin) / dT;

for time = 1:1:times
    % 领航者轨迹
    ul(:,time+1) = [5 5]' -3 * cos( (T(:,time)+20)*pi/40 ) * [1 0]'...
                          -3 * sin( (T(:,time)+20)*pi/40 ) * [0 1]';
    vl(:,time+1) =         3 * sin( (T(:,time)+20)*pi/40 ) * [1 0]'...
                          -3 * cos( (T(:,time)+20)*pi/40 ) * [0 1]';
                      
    for i=1:4
        Sum1 = 0;
        Par2 = 0;
        Sum3 = 0;
        Sum4 = 0;
        Par5 = 0;
        Sum6 = 0;
        Sum7 = 0;
        Par8 = 0;
        
        % 第一行
        for j=1:4
            % 如果检测到两者之间符合通信关系
            if j~=i && abs(PX(i,time)-PX(j,time)) <= R
                % 计算梯度与势函数乘积求和
                Sum1 = Sum1 + nable_fun(PX(i,time),PX(j,time),d,R);% * Vij_Fun(PX(i,time),PX(j,time));
            end
        end
        if abs(PX(i,time)-ul(1,time)) >= R
            Par2 = nable_fun(PX(i,time),ul(1,time),d,R);% * Vij_Fun(PX(i,time),ul(1,time));
        end
        
        % 第二行
        for j=1:4
        if j~=i && abs(PX(i,time)-PX(j,time)) <= R
            for k=1:4
            if k~=i && abs(PX(i,time)-PX(k,time)) <= R
                Sum4 = Sum4 + (VX(i,time)-VX(k,time));
            end
            end
            if abs(PX(i,time)-ul(1,time)) >= R
                Par5 = VX(i,time) - vl(1,time);
            end
            Sum3 = Sum3 + sign(Sum4 + Par5);
        end
        end
        
        % 第三行
        for j=1:4
        if j~=i && abs(PX(i,time)-PX(j,time)) <= R
            for k=1:4
            if k~=i && abs(PX(i,time)-PX(k,time)) <= R
                Sum7 = Sum7 + (VX(j,time)-VX(k,time));
            end
            end
            if abs(PX(j,time)-ul(1,time)) >= R
                Par8 = VX(j,time) - vl(1,time);
            end
            Sum6 = Sum6 + sign(Sum7 + Par8);
        end
        end

        % 计算控制输入
        UX(i,time+1) = -Sum1 - Par2 - alpha*Sum3 + alpha*Sum6;
    end

    for i=1:4
        Sum1 = 0;
        Par2 = 0;
        Sum3 = 0;
        Sum4 = 0;
        Par5 = 0;
        Sum6 = 0;
        Sum7 = 0;
        Par8 = 0;
        
        % 第一行
        for j=1:4
            % 如果检测到两者之间符合通信关系
            if j~=i && abs(PY(i,time)-PY(j,time)) <= R
                % 计算梯度与势函数乘积求和
                Sum1 = Sum1 + nable_fun(PY(i,time),PY(j,time),d,R);% * Vij_Fun(PX(i,time),PX(j,time));
            end
        end
        if abs(PY(i,time)-ul(2,time)) >= R
            Par2 = nable_fun(PY(i,time),ul(2,time),d,R);% * Vij_Fun(PX(i,time),ul(1,time));
        end
        
        % 第二行
        for j=1:4
        if j~=i && abs(PY(i,time)-PY(j,time)) <= R
            for k=1:4
            if k~=i && abs(PY(i,time)-PY(k,time)) <= R
                Sum4 = Sum4 + (VY(i,time)-VY(k,time));
            end
            end
            if abs(PY(i,time)-ul(2,time)) >= R
                Par5 = VY(i,time) - vl(2,time);
            end
            Sum3 = Sum3 + sign(Sum4 + Par5);
        end
        end
        
        % 第三行
        for j=1:4
        if j~=i && abs(PY(i,time)-PY(j,time)) <= R
            for k=1:4
            if k~=i && abs(PY(i,time)-PY(k,time)) <= R
                Sum7 = Sum7 + (VY(j,time)-VY(k,time));
            end
            end
            if abs(PY(j,time)-ul(2,time)) >= R
                Par8 = VY(j,time) - vl(2,time);
            end
            Sum6 = Sum6 + sign(Sum7 + Par8);
        end
        end

        % 计算控制输入
        UY(i,time+1) = -Sum1 - Par2 - alpha*Sum3 + alpha*Sum6;
    end
    
        
    % 更新 X 参数
    VX(:,time+1) = VX(:,time) + dT * UX(:,time);
    PX(:,time+1) = PX(:,time) + dT * VX(:,time);
    
    % 更新 Y 参数
    VY(:,time+1) = VY(:,time) + dT * UY(:,time);
    PY(:,time+1) = PY(:,time) + dT * VY(:,time);
    
    % 记录时间
    T(:, time+1) = T(:, time) + dT;
end

% 绘制图像
figure(1)
plot3(T(:,1:times), ul(1,1:times), ul(2,1:times));
xlabel('T');
ylabel('X');
zlabel('Y');

% figure(2)
% subplot(2,1,1);
% plot( T(:,1:times),PX(1,1:times),...
%       T(:,1:times),PX(2,1:times),...
%       T(:,1:times),PX(3,1:times),...
%       T(:,1:times),PX(4,1:times),...
%       T(:,1:times),ul(1,1:times), 'linewidth',1.5); 
% legend('f1','f2','f3','f4', 'L');
% xlabel('T');
% ylabel('X');
% title('Position');
% grid on
% 
% subplot(2,1,2); 
% plot( T(:,1:times),VX(1,1:times),...
%       T(:,1:times),VX(2,1:times),...
%       T(:,1:times),VX(3,1:times),...
%       T(:,1:times),VX(4,1:times),'linewidth',1.5); 
% legend('f1','f2','f3','f4');
% xlabel('T');
% ylabel('X');
% title('Speed');
% grid on

figure(2)
subplot(2,1,1); 
plot3(T(:,1:times),PX(1,1:times),PY(1,1:times),...
      T(:,1:times),PX(2,1:times),PY(2,1:times),...
      T(:,1:times),PX(3,1:times),PY(3,1:times),...
      T(:,1:times),PX(4,1:times),PY(4,1:times),...
      T(:,1:times),ul(1,1:times),ul(2,1:times), 'linewidth',1.5); 
legend('f1','f2','f3','f4', 'L');
xlabel('T');
ylabel('X');
zlabel('Y');
title('Position');
grid on

subplot(2,1,2); 
plot3(T(:,1:times),VX(1,1:times),VY(1,1:times),...
      T(:,1:times),VX(2,1:times),VY(2,1:times),...
      T(:,1:times),VX(3,1:times),VY(3,1:times),...
      T(:,1:times),VX(4,1:times),VY(4,1:times),'linewidth',1.5); 
legend('f1','f2','f3','f4');
xlabel('T');
ylabel('X');
zlabel('Y');
title('Speed');
grid on




% 势函数
function Vij = Vij_Fun(xi, xj)
% 期望距离：d
% 通信距离：R
d = 0.5;
R = 1.0;
xij = abs(xi - xj);
Vij = (xij - d)^2 * (R - xij) / (xij + (R-xij)/240)...
    + (xij) * (xij - d)^2 / (R-xij+xij/240);
end

% 计算梯度函数
function nabla_PX = nable_fun(xi, xj, d, R)
    nabla_PX = (sign(xi - xj)*(d - abs(xi - xj))^2)/(R - (239*abs(xi - xj))/240) - (sign(xi - xj)*(d - abs(xi - xj))^2)/(R/240 + (239*abs(xi - xj))/240) - (239*sign(xi - xj)*(R - abs(xi - xj))*(d - abs(xi - xj))^2)/(240*(R/240 + (239*abs(xi - xj))/240)^2) - (2*abs(xi - xj)*sign(xi - xj)*(d - abs(xi - xj)))/(R - (239*abs(xi - xj))/240) + (239*abs(xi - xj)*sign(xi - xj)*(d - abs(xi - xj))^2)/(240*(R - (239*abs(xi - xj))/240)^2) - (2*sign(xi - xj)*(R - abs(xi - xj))*(d - abs(xi - xj)))/(R/240 + (239*abs(xi - xj))/240);
end













