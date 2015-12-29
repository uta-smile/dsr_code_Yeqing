clear all; close all; clc;
load mirt2D_data2.mat;
% load T;
load X;
T = res1.X;

% Main settings
main.similarity='dtv';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI
% main.mu=0.01;        % similarity measure parameter (e.g., alpha of RC)
main.TV = TVOP;
main.subdivide=3;       % use 3 hierarchical levels
main.okno=8;            % mesh window size, the smaller it is the more complex deformations are possible
main.lambda = 0.005;    % transformation regularization weight, 0 for none
main.single=0;          % show mesh transformation at every iteration

% Optimization settings
optim.maxsteps = 300;   % maximum number of iterations at each hierarchical level
optim.fundif = 1e-6;    % tolerance (stopping criterion)
optim.gamma = 1;       % initial optimization step size
optim.anneal=0.8;       % annealing rate on the optimization step

main2 = main;
% Main settings
main2.similarity='rc';   % similarity measure, e.g. SSD, CC, SAD, RC, CD2, MS, MI
main2.alpha=0.05;        % similarity measure parameter (e.g., alpha of RC)

trans.okno = 8;


K = 1:6;Iter=2; % it can be changed, 50 runs used in the paper
for k=K
    for i = 1:Iter
%         trans.X = T + 2*randn(size(T));
        trans.X = randpertube(T);
        im = refim;
        refim=mirt2D_transform(im, trans);
        
        refim_dis = add_distor(refim,k);

        
        t0=cputime();
        [res1, newim1]=mirt2D_register(refim_dis,im, main, optim);
        
        t1(k,i)=cputime()-t0;
        I_RMSE1(k,i) = RMSE(newim1,refim);
        T_RMSE1(k,i) = RMSE(res1.X,trans.X);
        
        t0=cputime();
        [res2, newim2]=mirt2D_register(refim_dis,im, main2, optim);
        t2(k,i)=cputime()-t0;
        I_RMSE2(k,i) = RMSE(newim2,refim);
        T_RMSE2(k,i) = RMSE(res2.X,trans.X);
        
    end
end

M_RMSE1 = mean(I_RMSE1,2); S_RMSE1 = std(I_RMSE1,1,2);
M_RMSE2 = mean(I_RMSE2,2); S_RMSE2 = std(I_RMSE2,1,2);

MT_RMSE1 = mean(T_RMSE1,2); ST_RMSE1 = std(T_RMSE1,1,2);
MT_RMSE2 = mean(T_RMSE2,2); ST_RMSE2 = std(T_RMSE2,1,2);


ls=3; ms=8; ts=20;
figure; hold on;box on;
errorbar(K,M_RMSE2, S_RMSE2, 'b*-.', 'linewidth',ls, 'markersize', ms);
errorbar(K,M_RMSE1, S_RMSE1, 'ro-', 'linewidth',ls, 'markersize', ms);
legend('RC','Proposed','Location','SouthEast');
xlabel('K');
ylabel('Intensity RMSE');
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts);
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts);
set(gca,'FontSize',16);
% ls=2; ms=8; ts=12;
figure; hold on;box on;
errorbar(K,MT_RMSE2, ST_RMSE2, 'b*-.', 'linewidth',ls, 'markersize', ms);
errorbar(K,MT_RMSE1, ST_RMSE1, 'ro-', 'linewidth',ls, 'markersize', ms);
legend('RC','Proposed','Location','SouthEast');
xlabel('K');
ylabel('Transformation RMSE');
textobj = findobj('type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(gca,'XLabel');
set(h_xlabel,'FontSize',ts);
h_xlabel = get(gca,'YLabel');
set(h_xlabel,'FontSize',ts);
set(gca,'FontSize',16);
% axis tight;


Time1 = mean(t1(:))
Time2 = mean(t2(:))
