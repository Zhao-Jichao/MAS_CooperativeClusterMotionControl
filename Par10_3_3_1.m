% 仿真实验
% 加上了算法 6 的实验程序
clear;
clc;

A = [0 1 0 0 0 1 1 1;
     1 0 1 0 0 1 0 1;
     0 1 0 1 1 0 0 1;
     0 0 1 0 1 0 1 1;
     0 0 1 1 0 1 1 0;
     1 1 0 0 1 0 1 0;
     1 0 0 1 1 1 0 0;
     1 1 1 1 0 0 0 0;];
 
D = [4 0 0 0 0 0 0 0;
     0 4 0 0 0 0 0 0;
     0 0 4 0 0 0 0 0;
     0 0 0 4 0 0 0 0;
     0 0 0 0 4 0 0 0;
     0 0 0 0 0 4 0 0;
     0 0 0 0 0 0 4 0;
     0 0 0 0 0 0 0 4;];
 
L = D - A;
 
Z0 = [10, 4, 8, 2, 9, 3, 5, 7]';
Zt(:,1) = Z0;
Ut(:,1) = Z0;
Et = zeros(8, 1);

epsilon = 0.001;

% 时间参数
tbegin = 0;
tfinal = 2;
dT = 0.001;
T(:,1) = 0;

% 计算次数
times = (tfinal - tbegin) / dT;

for time = 1:1:times
    % 记录时间刻度
    T(:,  time+1) = T(:, time) + dT;

    % 记录状态
    Ut(:, time) = epsilon * L * Zt(:, time);
    % 定义发生非合作行为，时间在 0.2s，此时 time=201
    if time >= 201
        % Et(4, time) = -Ut(4, time);                       % 毁坏型
        % Et(4, time) = -Ut(4, time) - 0.005;               % 失控型
        Et(4, time) = -Ut(4, time) - 0.03*(rand-0.5);     % 干扰型
        
        Ut(4, time) = Ut(4, time) + Et(4, time);
    end
    Zt(:, time+1) = Zt(:, time) - Ut(:, time);
    
end


plot(T(1,:),Zt(1,:),'-', T(1,:),Zt(2,:),':', T(1,:),Zt(3,:),'--', T(1,:),Zt(4,:),'-.',...
     T(1,:),Zt(5,:), T(1,:),Zt(6,:), T(1,:),Zt(7,:), T(1,:),Zt(8,:),...
     'linewidth',1.5)

grid on
legend('x_1', 'x_2', 'x_3', 'x_4', 'x_5', 'x_6', 'x_7', 'x_8');


