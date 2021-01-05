% 数值仿真

clear

% 领航者参数
ul(:,1) = [5 2]';

% 跟随者参数
PXf1(:,1)= 3; PYf1(:,1)=-2; VXf1(:,1)= 0; VYf1(:,1)= 0;
PXf2(:,1)= 6; PYf2(:,1)= 4; VXf2(:,1)= 0; VYf2(:,1)= 0;
PXf3(:,1)=-5; PYf3(:,1)=-5; VXf3(:,1)= 0; VYf3(:,1)= 0;
PXf4(:,1)=-4; PYf4(:,1)= 3; VXf4(:,1)= 0; VYf4(:,1)= 0;

PX = [PXf1 PXf2 PXf3 PXf4]';
PY = [PYf1 PYf2 PYf3 PYf4]';

VX = [VXf1 VXf2 VXf3 VXf4]';
VY = [VYf1 VYf2 VYf3 VYf4]';

% 系统关系
L = [ 3 -1 -1 -1;
     -1  3 -1 -1;
     -1 -1  3 -1;
     -1 -1 -1  3;];

% 时间参数
tbegin = 0;
tfinal = 60;
dT = 0.1;
T(:,1) = 0; 

% 计算次数
times = (tfinal - tbegin) / dT;

for time = 1:1:times
    % 记录时间
    T(:, time+1) = T(:, time) + dT;
    
    % 领航者轨迹
    ul(:,time) = [5 5]' -3 * cos( (T(:,time)+20)*pi/40 ) * [1 0]'...
                          -3 * sin( (T(:,time)+20)*pi/40 ) * [0 1]';
    
    % 跟随者轨迹
    alpha = 1;
    beta  = 0.01;    
    
    ufX = (alpha * ((-L) * PX(:,time) - (PX(:,time)-ul(1,time) - 0.1))...
        + beta * (-L) * VX(:,time));
    VX(:,time+1) = ufX;
    PX(:,time+1) = PX(:,time) + dT * VX(:,time+1);  % 更新
    
    ufY = (alpha * ((-L) * PY(:,time) - (PY(:,time)-ul(2,time)))...
        + beta * (-L) * VY(:,time));
    VY(:,time+1) = ufY;
    PY(:,time+1) = PY(:,time) + dT * VY(:,time+1);  % 更新
end

% 绘制图像
% figure(1)
% plot3(T(:,1:times), ul(1,1:times), ul(2,1:times));
% xlabel('T');
% ylabel('X');
% zlabel('Y');

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


u1 = - (Deltaxi*V(x1-x2) + Deltaxi*V(x1-x3) + Deltaxi*V(x1-x4))...
     - alpha * L(1,2) * sign( L(1,2)*V(v1-v2) + L(1,3)*V(v1-v3) + L(1,4)*V(v1-v4))...
     - alpha * L(1,3) * sign( L(1,3)*V(v1-v2) + L(1,3)*V(v1-v3) + L(1,4)*V(v1-v4))...
     - alpha * L(1,4) * sign( L(1,2)*V(v1-v2) + L(1,3)*V(v1-v3) + L(1,4)*V(v1-v4))...
     + alpha * L(1,2) * sign( L(2,1)*V(v2-v1) + L(2,4)*V(v2-v4))...
     + alpha * L(1,3) * sign( L(3,1)*V(v3-v1))...
     + alpha * L(1,4) * sign( L(4,1)*V(v4-v1) + L(4,2)*V(v4-v2));

u2 = - (Deltaxi*V(x2-x1) + Deltaxi*V(x2-x4))...
     - alpha * L(2,1) * sign( L(2,1)*V(v2-v1) + L(2,4)*V(v2-v4))...
     - alpha * L(2,4) * sign( L(2,1)*V(v2-v1) + L(2,4)*V(v2-v4))...
     + alpha * L(2,1) * sign( L(1,2)*V(v1-v2) + L(1,3)*V(v1-v3) + L(1,4)*V(v1-v4))...
     + alpha * L(2,4) * sign( L(4,1)*V(v4-v1) + L(4,2)*V(v4-v2));

u3 = - (Deltaxi*V(x3-x1))...
     - alpha * L(3,1) * sign( L(3,1)*V(v3-v1))...
     + alpha * L(3,1) * sign( L(1,2)*V(v1-v2) + L(1,3)*V(v1-v3) + L(1,4)*V(v1-v4));

u4 = - (Deltaxi*V(x4-x1) + Deltaxi*V(x4-x2))...
     - alpha * L(4,1) * sign( L(4,1)*V(v4-v1) + L(4,2)*V(v4-v2))...
     - alpha * L(4,2) * sign( L(4,1)*V(v4-v1) + L(4,2)*V(v4-v2))...
     + alpha * L(4,1) * sign( L(1,2)*V(v1-v2) + L(1,3)*V(v1-v3) + L(1,4)*V(v1-v4))...
     + alpha * L(4,2) * sign( L(2,1)*V(v2-v1) + L(2,4)*V(v2-v4));
 
xi=1;
xj=2;
xl=3;
t = xi-xj;
V(xi-xj)

Vij = V(xi-xj);
Vil = V(xi-xl);


% 时间参数
tbegin = 0;
tfinal = 2;
dT = 0.01;
T(:,1) = tbegin; 

% 计算次数
times = (tfinal - tbegin) / dT;

for time = 1:1:times
%     记录时间
    T(1, time+1) = T(1, time) + dT;
    T(2, time) = V(T(1, time));
end

plot(T(1,:),T(2,:))

function Vij = V(xij)
    % 参数
    d = 0.5;
    R = 1;
    c1 = 9.7;
    c2 = 9.7;
    psiMax = 50.3;
    % 运算程序
    Vij = ( (xij-d)^2 * (R-xij) )/( xij+d^2*(R-xij)/(c1+psiMax) )...
        + ( xij*(xij-d)^2 )/( (R-xij)+xij*(R-d)^2/(c2+psiMax) );
end



function ui = ui(Vij)
    % 参数
    L = [3 -1 -1 -1;
        -1  2  0 -1;
        -1  0 -1  0;
        -1 -1  0  2;];
    
    Ll= [1  1  0  0]'; 
    
    Deltaxi = 1;
    alpha = 0.5;
    
    for i=1:length(L)
        Sum1 = 0;
        Sum2 = 0;
        Sum3 = 0;
        Sum4 = 0;
        Sum5 = 0;
        for j=1:length(L)
            if L(i,j)==-1
                Sum1 = Sum1 + Deltaxi * V(xi-xj);
                Sum2 = Sum2 + L(i,j) * sign();
                Sum3 = Sum3 + L(i,j) * sign();
            end
        end
    end
    
    
    ui = -sum * Delta_xi * Vij - hi * Delta_xi * Vil ...

end





















