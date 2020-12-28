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

epsilon = 0.001;

% 时间参数
tbegin = 0;
tfinal = 10;
dT = 0.001;
T(:,1) = 0; 

% 计算次数
times = (tfinal - tbegin) / dT;


for time = 1:1:times
    U = epsilon * L * Zt(:, time);
    Zt(:, time+1) = Zt(:, time) - U;
    T(:,  time+1) = T(:, time) + dT;
    
    % 定义发生非合作行为，时间在 0.2s，此时 time=201
    if time == 201
        L(4,:) = 0;
    end
end


plot(T(1,:),Zt(1,:),'-', T(1,:),Zt(2,:),':', T(1,:),Zt(3,:),'--', T(1,:),Zt(4,:),'-.',...
     T(1,:),Zt(5,:), T(1,:),Zt(6,:), T(1,:),Zt(7,:), T(1,:),Zt(8,:),...
     'linewidth',1)
 
grid on
legend('x_1', 'x_2', 'x_3', 'x_4', 'x_5', 'x_6', 'x_7', 'x_8');






