clear
clc
% 此程序在之前的基础工作上，完全采用书中公式所写
% Author: Jichao-Zhao
% Date: 2021-03-18

% 定义实际值
p1(:,1) = [1; 3];
p2(:,1) = [1; 6];
p3(:,1) = [7; 1];
p4(:,1) = [4; 6];
p5(:,1) = [4; 9];
p6(:,1) = [6; 2];
p7(:,1) = [7; 4];
p8(:,1) = [9; 2];
p9(:,1) = [7; 6];
p10(:,1)= [9; 9];
p11(:,1)= [4; 1];
p12(:,1)= [6; 9];

v1(:,1) = [ 1.5;  1.5];
v2(:,1) = [  -2;  2.5];
v3(:,1) = [ 1.2;  1.5];
v4(:,1) = [ 1.5; -1.5];
v5(:,1) = [   0;    3];
v6(:,1) = [-2.4;  1.2];
v7(:,1) = [   3;    0];
v8(:,1) = [ 1.8;  1.2];
v9(:,1) = [ 0.8; -0.4];
v10(:,1)= [-0.6;   -3];
v11(:,1)= [ 2.3; -1.2];
v12(:,1)= [ 2.5;  2.5];

% 定义系统理论值
p1T(:,1) = [1; 3];
p2T(:,1) = [1; 6];
p3T(:,1) = [7; 1];
p4T(:,1) = [4; 6];
p5T(:,1) = [4; 9];
p6T(:,1) = [6; 2];
p7T(:,1) = [7; 4];
p8T(:,1) = [9; 2];
p9T(:,1) = [7; 6];
p10T(:,1)= [9; 9];
p11T(:,1)= [4; 1];
p12T(:,1)= [6; 9];

v1T(:,1) = [ 1.5;  1.5];
v2T(:,1) = [  -2;  2.5];
v3T(:,1) = [ 1.2;  1.5];
v4T(:,1) = [ 1.5; -1.5];
v5T(:,1) = [   0;    3];
v6T(:,1) = [-2.4;  1.2];
v7T(:,1) = [   3;    0];
v8T(:,1) = [ 1.8;  1.2];
v9T(:,1) = [ 0.8; -0.4];
v10T(:,1)= [-0.6;   -3];
v11T(:,1)= [ 2.3; -1.2];
v12T(:,1)= [ 2.5;  2.5];

% 定义系统测量值
p1M(:,1) = [1; 3];
p2M(:,1) = [1; 6];
p3M(:,1) = [7; 1];
p4M(:,1) = [4; 6];
p5M(:,1) = [4; 9];
p6M(:,1) = [6; 2];
p7M(:,1) = [7; 4];
p8M(:,1) = [9; 2];
p9M(:,1) = [7; 6];
p10M(:,1)= [9; 9];
p11M(:,1)= [4; 1];
p12M(:,1)= [6; 9];

v1M(:,1) = [ 1.5;  1.5];
v2M(:,1) = [  -2;  2.5];
v3M(:,1) = [ 1.2;  1.5];
v4M(:,1) = [ 1.5; -1.5];
v5M(:,1) = [   0;    3];
v6M(:,1) = [-2.4;  1.2];
v7M(:,1) = [   3;    0];
v8M(:,1) = [ 1.8;  1.2];
v9M(:,1) = [ 0.8; -0.4];
v10M(:,1)= [-0.6;   -3];
v11M(:,1)= [ 2.3; -1.2];
v12M(:,1)= [ 2.5;  2.5];

A =[0 1 0 0 0 0 0 0 0 0 1 0;
    1 0 0 1 1 0 0 0 0 0 0 0;
    0 0 0 0 0 1 1 1 0 0 0 0;
    0 1 0 0 1 0 1 0 1 0 1 0;
    0 1 0 1 0 0 0 0 0 0 0 1;
    0 0 1 0 0 0 1 0 0 0 1 0;
    0 0 1 1 0 1 0 1 1 0 0 0;
    0 0 1 0 0 0 1 0 0 1 0 0;
    0 0 0 1 0 0 1 0 0 1 0 1;
    0 0 0 0 0 0 0 1 1 0 0 1;
    1 0 0 1 0 1 0 0 0 0 0 0;
    0 0 0 0 1 0 0 0 1 0 0 1;];

t(1,1)= 0;

lambda = 3.3;
mu     = 2.6;
alpha  = 3.8;
beta   = 2.1;
kappa  = 2.3;

% time parameters
tBegin = 0;
tEnd   = 20;
dT     = 0.1;
times  = (tEnd-tBegin)/dT;
Anomaly = 80;
Detect  = 100;
Repair  = 120;

% global N
for N = 1:1:times
    t(:,N+1)  = t(:,N) + dT;
    
    u1 = input(p1(:,N), v1(:,N), [p2(:,N);p11(:,N)], [v2(:,N);v11(:,N)]);
    % 理论值
    v1T(:,N+1) = v1T(:,N) + dT * u1; 
    p1T(:,N+1) = p1T(:,N) + dT * v1T(:,N+1);
    % 实际值
    v1(:,N+1) = v1(:,N) + dT * u1; 
    p1(:,N+1) = p1(:,N) + dT * v1(:,N+1);
    % 测量值
    v1M(:,N+1) = v1(:,N+1); 
    p1M(:,N+1) = p1(:,N+1);
    
    % Choose Leader
    % u2 = input(p2(:,N), v2(:,N), [p1(:,N);p4(:,N);p5(:,N)], [v1(:,N);v4(:,N);v5(:,N)]);
    u2 = 0;
    % 理论值
    v2T(:,N+1) = v2T(:,N) + dT * u2; 
    p2T(:,N+1) = p2T(:,N) + dT * v2T(:,N+1);
    % 实际值
    v2(:,N+1) = v2(:,N) + dT * u2; 
    p2(:,N+1) = p2(:,N) + dT * v2(:,N+1);
    % 测量值
    v2M(:,N+1) = v2(:,N+1); 
    p2M(:,N+1) = p2(:,N+1);
    
    u3 = input(p3(:,N), v3(:,N), [p6(:,N);p7(:,N);p8(:,N)], [v6(:,N);v7(:,N);v8(:,N)]);
    % 理论值
    v3T(:,N+1) = v3T(:,N) + dT * u3; 
    p3T(:,N+1) = p3T(:,N) + dT * v3T(:,N+1);
    % 第 12s 开始修复
    if N>Repair
        u3 = input(p3(:,N), v3(:,N), [p6(:,N);p8(:,N)], [v6(:,N);v8(:,N)]);    
    end
    % 实际值
    v3(:,N+1) = v3(:,N) + dT * u3; 
    p3(:,N+1) = p3(:,N) + dT * v3(:,N+1);
    % 测量值
    v3M(:,N+1) = v3(:,N+1); 
    p3M(:,N+1) = p3(:,N+1);

    u4 = input(p4(:,N), v4(:,N), [p2(:,N);p5(:,N);p7(:,N);p9(:,N);p11(:,N)], [v2(:,N);v5(:,N);v7(:,N);v9(:,N);v11(:,N)]);
    % 理论值
    v4T(:,N+1) = v4T(:,N) + dT * u4; 
    p4T(:,N+1) = p4T(:,N) + dT * v4T(:,N+1);
    if N>Repair
        u4 = input(p4(:,N), v4(:,N), [p2(:,N);p5(:,N);p9(:,N);p11(:,N)], [v2(:,N);v5(:,N);v9(:,N);v11(:,N)]);
    end
    % 实际值
    v4(:,N+1) = v4(:,N) + dT * u4; 
    p4(:,N+1) = p4(:,N) + dT * v4(:,N+1);
    % 测量值
    v4M(:,N+1) = v4(:,N+1); 
    p4M(:,N+1) = p4(:,N+1);
    
    
    
    
    u5 = input(p5(:,N), v5(:,N), [p2(:,N);p4(:,N);p4(:,N)], [v2(:,N);v4(:,N);v12(:,N)]);
    % 理论值
    v5T(:,N+1) = v5T(:,N) + dT * u5; 
    p5T(:,N+1) = p5T(:,N) + dT * v5T(:,N+1);
    % 实际值
    v5(:,N+1) = v5(:,N) + dT * u5; 
    p5(:,N+1) = p5(:,N) + dT * v5(:,N+1);
    % 测量值
    v5M(:,N+1) = v5(:,N+1); 
    p5M(:,N+1) = p5(:,N+1);

    u6 = input(p6(:,N),v6(:,N), [p3(:,N);p7(:,N);p11(:,N)], [v3(:,N);v7(:,N);v11(:,N)]);
    % 理论值
    v6T(:,N+1) = v6T(:,N) + dT * u6; 
    p6T(:,N+1) = p6T(:,N) + dT * v6T(:,N+1);
    if N>Repair
        u6 = input(p6(:,N),v6(:,N), [p3(:,N);p11(:,N)], [v3(:,N);v11(:,N)]);
    end
    % 实际值
    v6(:,N+1) = v6(:,N) + dT * u6; 
    p6(:,N+1) = p6(:,N) + dT * v6(:,N+1);
    % 测量值
    v6M(:,N+1) = v6(:,N+1); 
    p6M(:,N+1) = p6(:,N+1);
    
    
    % 非合作节点
    u7 = input(p7(:,N), v7(:,N), [p3(:,N);p4(:,N);p6(:,N);p8(:,N);p9(:,N)], [v3(:,N);v4(:,N);v6(:,N);v8(:,N);v9(:,N)]);
    v7(:,N+1) = v7(:,N) + dT * u7; 
    % 异常发生
    if N>Anomaly
        u7 = 0;
        v7(:,N+1) = 0;
    end
    p7(:,N+1) = p7(:,N) + dT * v7(:,N+1);
    
    
    u8 = input(p8(:,N), v8(:,N), [p3(:,N);p7(:,N);p10(:,N)], [v3(:,N);v7(:,N);v10(:,N)]);
    % 理论值
    v8T(:,N+1) = v8T(:,N) + dT * u8; 
    p8T(:,N+1) = p8T(:,N) + dT * v8T(:,N+1);
    if N>Repair
        u8 = input(p8(:,N), v8(:,N), [p3(:,N);p10(:,N)], [v3(:,N);v10(:,N)]);
    end
    % 实际值
    v8(:,N+1) = v8(:,N) + dT * u8; 
    p8(:,N+1) = p8(:,N) + dT * v8(:,N+1);
    % 测量值
    v8M(:,N+1) = v8(:,N+1); 
    p8M(:,N+1) = p8(:,N+1);

    
    u9 = input(p9(:,N), v9(:,N), [p4(:,N);p7(:,N);p10(:,N);p12(:,N)], [v4(:,N);v7(:,N);v10(:,N);v12(:,N)]);
    % 理论值
    v9T(:,N+1) = v9T(:,N) + dT * u9; 
    p9T(:,N+1) = p9T(:,N) + dT * v9T(:,N+1);
    if N>Repair
        u9 = input(p9(:,N), v9(:,N), [p4(:,N);p10(:,N);p12(:,N)], [v4(:,N);v10(:,N);v12(:,N)]);
    end
    % 实际值
    v9(:,N+1) = v9(:,N) + dT * u9; 
    p9(:,N+1) = p9(:,N) + dT * v9(:,N+1);
    % 测量值
    v9M(:,N+1) = v9(:,N+1); 
    p9M(:,N+1) = p9(:,N+1);
    
    
    u10 = input(p10(:,N), v10(:,N), [p8(:,N);p9(:,N);p12(:,N)], [v8(:,N);v9(:,N);v12(:,N)]);
    % 理论值
    v10T(:,N+1) = v10T(:,N) + dT * u10; 
    p10T(:,N+1) = p10T(:,N) + dT * v10T(:,N+1);
    % 实际值
    v10(:,N+1) = v10(:,N) + dT * u10; 
    p10(:,N+1) = p10(:,N) + dT * v10(:,N+1);
    % 测量值
    v10M(:,N+1) = v10(:,N+1); 
    p10M(:,N+1) = p10(:,N+1);
    
    u11 = input(p11(:,N), v11(:,N), [p1(:,N);p4(:,N);p6(:,N)], [v1(:,N);v4(:,N);v6(:,N)]);
    % 理论值
    v11T(:,N+1) = v11T(:,N) + dT * u11; 
    p11T(:,N+1) = p11T(:,N) + dT * v11T(:,N+1);
    % 实际值
    v11(:,N+1) = v11(:,N) + dT * u11; 
    p11(:,N+1) = p11(:,N) + dT * v11(:,N+1);
    % 测量值
    v11M(:,N+1) = v11(:,N+1); 
    p11M(:,N+1) = p11(:,N+1);
    
    u12 = input(p12(:,N), v12(:,N), [p5(:,N);p9(:,N);p10(:,N)], [v5(:,N);v9(:,N);v10(:,N)]);
    % 理论值
    v12T(:,N+1) = v12T(:,N) + dT * u12; 
    p12T(:,N+1) = p12T(:,N) + dT * v12T(:,N+1);
    % 实际值
    v12(:,N+1) = v12(:,N) + dT * u12; 
    p12(:,N+1) = p12(:,N) + dT * v12(:,N+1);
    % 测量值
    v12M(:,N+1) = v12(:,N+1); 
    p12M(:,N+1) = p12(:,N+1);

end

huitu1 = 0; % 时刻点图
huitu2 = 0; % 轨迹图
huitu3 = 1; % 速度随时间变化图

if huitu1 == 1
    global T 
    T = 101;
    figure(1); hold on
    plot(p1(1,T),p1(2,T),'o');   hold on
    plot(p2(1,T),p2(2,T),'o');   hold on
    plot(p3(1,T),p3(2,T),'o');   hold on
    plot(p4(1,T),p4(2,T),'o');   hold on
    plot(p5(1,T),p5(2,T),'o');   hold on
    plot(p6(1,T),p6(2,T),'o');   hold on
    plot(p7(1,T),p7(2,T),'>');   hold on
    plot(p8(1,T),p8(2,T),'o');   hold on
    plot(p9(1,T),p9(2,T),'o');   hold on
    plot(p10(1,T),p10(2,T),'o'); hold on
    plot(p11(1,T),p11(2,T),'o'); hold on
    plot(p12(1,T),p12(2,T),'o'); hold on
    
    quiver(p1(1,T), p1(2,T), v1(1,T), v1(2,T));
    quiver(p2(1,T), p2(2,T), v2(1,T), v2(2,T));
    quiver(p3(1,T), p3(2,T), v3(1,T), v3(2,T));
    quiver(p4(1,T), p4(2,T), v4(1,T), v4(2,T));
    quiver(p5(1,T), p5(2,T), v5(1,T), v5(2,T));
    quiver(p6(1,T), p6(2,T), v6(1,T), v6(2,T));
    quiver(p7(1,T), p7(2,T), v7(1,T), v7(2,T));
    quiver(p8(1,T), p8(2,T), v8(1,T), v8(2,T));
    quiver(p9(1,T), p9(2,T), v9(1,T), v9(2,T));
    quiver(p10(1,T),p10(2,T),v10(1,T),v10(2,T));
    quiver(p11(1,T),p11(2,T),v11(1,T),v11(2,T));
    quiver(p12(1,T),p12(2,T),v12(1,T),v12(2,T));
    
%     lineFun(p1,  [p2; p11]);
%     lineFun(p2,  [p1; p4; p5]);
%     lineFun(p3,  [p6; p7; p8]);
%     lineFun(p4,  [p2; p5; p7; p9; p11]);
%     lineFun(p4,  [p2; p5; p7; p9; p11]);
%     lineFun(p6,  [p3; p7; p11]);
%     lineFun(p7,  [p3; p4; p6; p8;p9]);
%     lineFun(p8,  [p3; p7; p10]);
%     lineFun(p9,  [p4; p7; p10;p12]);
%     lineFun(p10, [p8; p9; p12]);
%     lineFun(p11, [p1; p4; p6]);
%     lineFun(p12, [p5; p9; p10]);
    
    lineFun(p1,  [p2; p11]);
    lineFun(p2,  [p1; p4; p5]);
    lineFun(p3,  [p6; p8]);
    lineFun(p4,  [p2; p5; p9; p11]);
    lineFun(p4,  [p2; p5; p9; p11]);
    lineFun(p6,  [p3;  p11]);
%     lineFun(p7,  [p3; p4; p6; p8;p9]);
    lineFun(p8,  [p3; p10]);
    lineFun(p9,  [p4; p10;p12]);
    lineFun(p10, [p8; p9; p12]);
    lineFun(p11, [p1; p4; p6]);
    lineFun(p12, [p5; p9; p10]);
    
    text(p1(1,T)-0.3,p1(2,T)-0.3,'1'); 
    text(p2(1,T)-0.3,p2(2,T)-0.3,'2');
    text(p3(1,T)-0.3,p3(2,T)-0.3,'3');
    text(p4(1,T)-0.3,p4(2,T)-0.3,'4');
    text(p5(1,T)-0.3,p5(2,T)-0.3,'5');
    text(p6(1,T)-0.3,p6(2,T)-0.3,'6');
    text(p7(1,T)-0.3,p7(2,T)-0.3,'7');
    text(p8(1,T)-0.3,p8(2,T)-0.3,'8');
    text(p9(1,T)-0.3,p9(2,T)-0.3,'9'); 
    text(p10(1,T)-0.3,p10(2,T)-0.3,'10');
    text(p11(1,T)-0.3,p11(2,T)-0.3,'11');
    text(p12(1,T)-0.3,p12(2,T)-0.3,'12');
end
% pause(0.5)
% end

if huitu2 == 1
    figure(2)
    plot(p1(1,:),p1(2,:),'color','r'); hold on
    plot(p2(1,:),p2(2,:)); hold on
    plot(p3(1,:),p3(2,:)); hold on
    plot(p4(1,:),p4(2,:)); hold on
    plot(p5(1,:),p5(2,:)); hold on
    plot(p6(1,:),p6(2,:)); hold on
    plot(p7(1,:),p7(2,:),'color','b','linewidth',1.5); hold on
    plot(p8(1,:),p8(2,:)); hold on
    plot(p9(1,:),p9(2,:)); hold on
    plot(p10(1,:),p10(2,:)); hold on
    plot(p11(1,:),p11(2,:)); hold on
    plot(p12(1,:),p12(2,:)); hold on
end

if huitu3 == 1
    plot(t,sqrt(v1(1,:).^2+v1(2,:).^2)); hold on
    plot(t,sqrt(v2(1,:).^2+v2(2,:).^2)); hold on
    plot(t,sqrt(v3(1,:).^2+v3(2,:).^2)); hold on
    plot(t,sqrt(v4(1,:).^2+v4(2,:).^2)); hold on
    plot(t,sqrt(v5(1,:).^2+v5(2,:).^2)); hold on
    plot(t,sqrt(v6(1,:).^2+v6(2,:).^2)); hold on
    plot(t,sqrt(v7(1,:).^2+v7(2,:).^2)); hold on
    plot(t,sqrt(v8(1,:).^2+v8(2,:).^2)); hold on
    plot(t,sqrt(v9(1,:).^2+v9(2,:).^2)); hold on
    plot(t,sqrt(v10(1,:).^2+v10(2,:).^2)); hold on
    plot(t,sqrt(v11(1,:).^2+v11(2,:).^2)); hold on
    plot(t,sqrt(v12(1,:).^2+v12(2,:).^2)); hold on
end

% 控制输入
function ui = input(pi, vi, pj,vj)
    lambda = 3.3;
    mu     = 2.6;
    alpha  = 3.8;
    beta   = 2.1;
    kappa  = 2.3;
    
    ui = 0;
    for k=1:1:(size(pj,1)/2)
        pk = pj( (2*k-1):2*k, :);
        vk = vj( (2*k-1):2*k, :);
        ui = ui + (lambda/norm(pi-pk,2)^alpha - mu/norm(pi-pk,2)^beta) * (pi-pk) - kappa * (vi-vk);
    end
end

% 绘制通信线段
function lineFun(pi, pj)
    global T
    pi = pi(:,T);                   % 取出对应时刻的坐标
    for k=1:1:(size(pj,1)/2)
        pk = pj( (2*k-1):2*k, :);
        pk = pk(:,T);
        line([pi(1,:),pk(1,:)],[pi(2,:),pk(2,:)]);
    end
end



