tic
simulation_data = load('test.mat');
num_antenna = simulation_data.num_antenna;
lb = -2*pi*ones(1,num_antenna/2);
ub = -lb;
options = optimoptions('ga','PlotFcn',"gaplotbestf");
options = optimoptions(options,'MaxGenerations',100);
step_phi = simulation_data.step_phi;
range = -4.83*pi/180:step_phi:4.83*pi/180;
lambda = simulation_data.lambda;
k = 2*pi/lambda;
z = simulation_data.z;
d = simulation_data.d;
phase_error =   simulation_data.phase_error;
best_ang = simulation_data.bestang;
bestang_stg5 = zeros(9,num_antenna/2);
delete(gcp('nocreate'));
% parpool(9);

comp = zeros(1,64);
parfor steer = -4:4
    func1=@(x) obj_fun_steer(x,steer);
    [x,fval]=ga(func1,num_antenna/2,[],[],[],[],lb,ub,[],options); 
    bestang_stg5(steer+5,:) = x;
end
delete(gcp('nocreate'));

phi0 = zeros(length(varphi_),num_antenna);
for steer_i = -4:4
    xi_ = steer_i*pi/180;
    for i = 1:num_antenna
        phi0(:,i) = (i-1)*(k*d*sin(range));
    end
    varphi_ = range-xi_;
    comp_ga(1:2:num_antenna) = bestang_stg5(steer_i+5,:);
    comp_ga(2:2:num_antenna) = bestang_stg5(steer_i+5,:);
    d_phi = comp_ga+best_ang(steer_i+5,:);
    best_u = getfarfieldpattern(phi0,d_phi,phase_error,lambda,z,varphi_,d);
    psll_f = getPSLL(best_u);
    writematrix(comp_ga,"result_pso.xlsx","Sheet","compensated angle("+num2str(steer_i)+")");
    writematrix(psll_f,"result_pso.xlsx","Sheet","Final PSLL("+num2str(steer_i)+")");
    figure(3);
    plot(varphi_*180/pi, (best_u/maxau).^2);
    legend("psll="+num2str(psll_f)+"db");
    hold on
    saveas(gcf,"last_ga.fig");
end
save ga.mat
toc

