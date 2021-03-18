clear
clc
% 此程序在之前的基础工作上，完全采用书中公式所写
% Author: Jichao-Zhao
% Date: 2021-03-18

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

lambda = 3.3;
mu     = 2.6;
alpha  = 3.8;
beta   = 2.1;
kappa  = 2.3;

% time parameters
tBegin = 0;
tEnd   = 20;
dT     = 0.01;
times  = (tEnd-tBegin)/dT;

% global N
for N = 1:1:times
    
    u1 = input(p1(:,N), v1(:,N), [p2(:,N);p11(:,N)], [v2(:,N);v11(:,N)]);
    
    v1(:,N+1) = v1(:,N) + dT * u1; 
    p1(:,N+1) = p1(:,N) + dT * v1(:,N+1);
    
    u2 = input(p2(:,N), v2(:,N), [p1(:,N);p4(:,N);p5(:,N)], [v1(:,N);v4(:,N);v5(:,N)]);
    % Choose Leader
    u2 = 0;
    v2(:,N+1) = v2(:,N) + dT * u2; 
    p2(:,N+1) = p2(:,N) + dT * v2(:,N+1);
    
    u3 = input(p3(:,N), v3(:,N), [p6(:,N);p7(:,N);p8(:,N)], [v6(:,N);v7(:,N);v8(:,N)]);
    
    v3(:,N+1) = v3(:,N) + dT * u3; 
    p3(:,N+1) = p3(:,N) + dT * v3(:,N+1);

    u4 = input5(p4(:,N), v4(:,N), p2(:,N),v2(:,N), p5(:,N),v5(:,N), p7(:,N),v7(:,N), p9(:,N),v9(:,N), p11(:,N),v11(:,N));
    v4(:,N+1) = v4(:,N) + dT * u4; 
    p4(:,N+1) = p4(:,N) + dT * v4(:,N+1);
    
    u5 = input3(p5(:,N), v5(:,N), p2(:,N),v2(:,N), p4(:,N),v4(:,N), p12(:,N),v12(:,N));
    v5(:,N+1) = v5(:,N) + dT * u5; 
    p5(:,N+1) = p5(:,N) + dT * v5(:,N+1);
    
    u6 = input3(p6(:,N),v6(:,N), p3(:,N), v3(:,N), p7(:,N),v7(:,N), p11(:,N),v11(:,N));
    v6(:,N+1) = v6(:,N) + dT * u6; 
    p6(:,N+1) = p6(:,N) + dT * v6(:,N+1);
    
    u7 = input5(p7(:,N), v7(:,N), p3(:,N),v3(:,N), p4(:,N),v4(:,N), p6(:,N),v6(:,N), p8(:,N),v8(:,N), p9(:,N),v9(:,N));
    v7(:,N+1) = v7(:,N) + dT * u7; 
    p7(:,N+1) = p7(:,N) + dT * v7(:,N+1);
    
    u8 = input3(p8(:,N), v8(:,N), p3(:,N),v3(:,N), p7(:,N),v7(:,N), p10(:,N),v10(:,N));
    v8(:,N+1) = v8(:,N) + dT * u8; 
    p8(:,N+1) = p8(:,N) + dT * v8(:,N+1);
    
    u9 = input4(p9(:,N), v9(:,N), p4(:,N),v4(:,N), p7(:,N),v7(:,N), p10(:,N),v10(:,N), p12(:,N),v12(:,N));
    v9(:,N+1) = v9(:,N) + dT * u9; 
    p9(:,N+1) = p9(:,N) + dT * v9(:,N+1);
    
    u10 = input3(p10(:,N), v10(:,N), p8(:,N),v8(:,N), p9(:,N),v9(:,N), p12(:,N),v12(:,N));
    v10(:,N+1) = v10(:,N) + dT * u10; 
    p10(:,N+1) = p10(:,N) + dT * v10(:,N+1);
    
    u11 = input3(p11(:,N), v11(:,N), p1(:,N),v1(:,N), p4(:,N),v4(:,N), p6(:,N),v6(:,N));
    v11(:,N+1) = v11(:,N) + dT * u11; 
    p11(:,N+1) = p11(:,N) + dT * v11(:,N+1);

    u12 = input3(p12(:,N), v12(:,N), p5(:,N),v5(:,N), p9(:,N),v9(:,N), p10(:,N),v10(:,N));
    v12(:,N+1) = v12(:,N) + dT * u12; 
    p12(:,N+1) = p12(:,N) + dT * v12(:,N+1);

end

huitu1 = 1;
huitu2 = 0;

if huitu1 == 1
    global T 
    T = 1;
    hold on
    plot(p1(1,T),p1(2,T),'o');   quiver(p1(1,T), p1(2,T), v1(1,T), v1(2,T));  lineFun(p1,  [p2;p11]);  hold on
    plot(p2(1,T),p2(2,T),'o');   quiver(p2(1,T), p2(2,T), v2(1,T), v2(2,T));  lineFun(p2,  [p1;p4;p5]);      hold on
    plot(p3(1,T),p3(2,T),'o');   quiver(p3(1,T), p3(2,T), v3(1,T), v3(2,T));  lineFun(p3,  [p6;p7;p8]);      hold on
    plot(p4(1,T),p4(2,T),'o');   quiver(p4(1,T), p4(2,T), v4(1,T), v4(2,T));  lineFun(p4,  [p2;p5;p7; p9;p11]); hold on
    plot(p5(1,T),p5(2,T),'o');   quiver(p5(1,T), p5(2,T), v5(1,T), v5(2,T));  lineFun(p5,  [p2;p4;p12]);     hold on
    plot(p6(1,T),p6(2,T),'o');   quiver(p6(1,T), p6(2,T), v6(1,T), v6(2,T));  lineFun(p6,  [p3;p7;p11]); hold on
    plot(p7(1,T),p7(2,T),'>');   quiver(p7(1,T), p7(2,T), v7(1,T), v7(2,T));  lineFun(p7,  [p3;p4;p6; p8;p9]); hold on
    plot(p8(1,T),p8(2,T),'o');   quiver(p8(1,T), p8(2,T), v8(1,T), v8(2,T));  lineFun(p8,  [p3;p7;p10]); hold on
    plot(p9(1,T),p9(2,T),'o');   quiver(p9(1,T), p9(2,T), v9(1,T), v9(2,T));  lineFun(p9,  [p4;p7;p10;p12]); hold on
    plot(p10(1,T),p10(2,T),'o'); quiver(p10(1,T),p10(2,T),v10(1,T),v10(2,T)); lineFun(p10, [p8;p9;p12]); hold on
    plot(p11(1,T),p11(2,T),'o'); quiver(p11(1,T),p11(2,T),v11(1,T),v11(2,T)); lineFun(p11, [p1;p4;p6]); hold on
    plot(p12(1,T),p12(2,T),'o'); quiver(p12(1,T),p12(2,T),v12(1,T),v12(2,T)); lineFun(p12, [p5;p9;p10]); hold on
    text(p1(1,T)+0.3,p1(2,T)+0.3,'1');  hold on
    text(p2(1,T)+0.3,p2(2,T)+0.3,'2');  hold on
    text(p3(1,T)+0.3,p3(2,T)+0.3,'3');  hold on
    text(p4(1,T)+0.3,p4(2,T)+0.3,'4');  hold on
    text(p5(1,T)+0.3,p5(2,T)+0.3,'5');  hold on
    text(p6(1,T)+0.3,p6(2,T)+0.3,'6');  hold on
    text(p7(1,T)+0.3,p7(2,T)+0.3,'7');  hold on
    text(p8(1,T)+0.3,p8(2,T)+0.3,'8');  hold on
    text(p9(1,T)+0.3,p9(2,T)+0.3,'9');  hold on
    text(p10(1,T)+0.3,p10(2,T)+0.3,'10');  hold on
    text(p11(1,T)+0.3,p11(2,T)+0.3,'11');  hold on
    text(p12(1,T)+0.3,p12(2,T)+0.3,'12');  hold on
end
% pause(0.5)
% end

if huitu2 == 1
    plot(p1(1,:),p1(2,:),'*','color','r'); hold on
    plot(p2(1,:),p2(2,:),'o'); hold on
    plot(p3(1,:),p3(2,:),'o'); hold on
    plot(p4(1,:),p4(2,:),'o'); hold on
    plot(p5(1,:),p5(2,:),'o'); hold on
    plot(p6(1,:),p6(2,:),'o'); hold on
    plot(p7(1,:),p7(2,:),'>','color','b'); hold on
    plot(p8(1,:),p8(2,:),'o'); hold on
    plot(p9(1,:),p9(2,:),'o'); hold on
    plot(p10(1,:),p10(2,:),'o'); hold on
    plot(p11(1,:),p11(2,:),'o'); hold on
    plot(p12(1,:),p12(2,:),'o'); hold on
end

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

function ui = input2(pi, vi, pj,vj, pk,vk)
    lambda = 3.3;
    mu     = 2.6;
    alpha  = 3.8;
    beta   = 2.1;
    kappa  = 2.3;
    
    ui = (lambda/norm(pi-pj,2)^alpha - mu/norm(pi-pj,2)^beta) * (pi-pj) - kappa * (vi-vj)...
        +(lambda/norm(pi-pk,2)^alpha - mu/norm(pi-pk,2)^beta) * (pi-pk) - kappa * (vi-vk);
end

function ui = input3(pi, vi, pj,vj, pk,vk, pl,vl)
    lambda = 3.3;
    mu     = 2.6;
    alpha  = 3.8;
    beta   = 2.1;
    kappa  = 2.3;
    
    ui = (lambda/norm(pi-pj,2)^alpha - mu/norm(pi-pj,2)^beta) * (pi-pj) - kappa * (vi-vj)...
        +(lambda/norm(pi-pk,2)^alpha - mu/norm(pi-pk,2)^beta) * (pi-pk) - kappa * (vi-vk)...
        +(lambda/norm(pi-pl,2)^alpha - mu/norm(pi-pl,2)^beta) * (pi-pl) - kappa * (vi-vl);
end

function ui = input4(pi, vi, pj,vj, pk,vk, pl,vl, pm,vm)
    lambda = 3.3;
    mu     = 2.6;
    alpha  = 3.8;
    beta   = 2.1;
    kappa  = 2.3;
    
    ui = (lambda/norm(pi-pj,2)^alpha - mu/norm(pi-pj,2)^beta) * (pi-pj) - kappa * (vi-vj)...
        +(lambda/norm(pi-pk,2)^alpha - mu/norm(pi-pk,2)^beta) * (pi-pk) - kappa * (vi-vk)...
        +(lambda/norm(pi-pl,2)^alpha - mu/norm(pi-pl,2)^beta) * (pi-pl) - kappa * (vi-vl)...
        +(lambda/norm(pi-pm,2)^alpha - mu/norm(pi-pm,2)^beta) * (pi-pm) - kappa * (vi-vm);
end

function ui = input5(pi, vi, pj,vj, pk,vk, pl,vl, pm,vm, pn,vn)
    lambda = 3.3;
    mu     = 2.6;
    alpha  = 3.8;
    beta   = 2.1;
    kappa  = 2.3;
    
    ui = (lambda/norm(pi-pj,2)^alpha - mu/norm(pi-pj,2)^beta) * (pi-pj) - kappa * (vi-vj)...
        +(lambda/norm(pi-pk,2)^alpha - mu/norm(pi-pk,2)^beta) * (pi-pk) - kappa * (vi-vk)...
        +(lambda/norm(pi-pl,2)^alpha - mu/norm(pi-pl,2)^beta) * (pi-pl) - kappa * (vi-vl)...
        +(lambda/norm(pi-pm,2)^alpha - mu/norm(pi-pm,2)^beta) * (pi-pm) - kappa * (vi-vm)...
        +(lambda/norm(pi-pn,2)^alpha - mu/norm(pi-pn,2)^beta) * (pi-pn) - kappa * (vi-vn);
end


function lineFun(pi, pj)
    global T
    pi = pi(:,T);                   % 取出对应时刻的坐标
    for k=1:1:(size(pj,1)/2)
        pk = pj( (2*k-1):2*k, :);
        pk = pk(:,T);
        line([pi(1,:),pk(1,:)],[pi(2,:),pk(2,:)]);
    end
end



