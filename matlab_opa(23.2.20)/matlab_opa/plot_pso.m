simulation_data = load('test.mat');
phi0 = simulation_data.phi0;
d_phi = simulation_data.d_phi;
dphi = simulation_data.dphi;
phase_error = simulation_data.phase_error;
lambda = simulation_data.lambda;
d = simulation_data.d;
z = simulation_data.z;
varphi = simulation_data.varphi;
dd_phi = load('pso.mat').x;
step = simulation_data.step;
au = getfarfieldpattern(phi0,d_phi+dphi+dd_phi,phase_error,lambda,z,varphi,d);
% [peakmainlobe, mainlobepeakpos] = max(au);
% 
% mainlobepos  = round(mainlobepeakpos-step/2):round(mainlobepeakpos+step/2);
% mainlobe = au(mainlobepos);

psll = getPSLL(au);
xlswrite("result.xlsx",dd_phi,"compensate_angle_pso");
xlswrite("result.xlsx",psll,"final psll after pso");
figure(21);
plot(varphi*180/pi,au.^2);
saveas(gcf,'pso_result1.fig');
% figure(22);
% plot((mainlobepos-length(varphi)/2)*step_phi*180/pi, mainlobe.^2/length(mainlobe));
% saveas(gcf,'pso_result2.fig');
