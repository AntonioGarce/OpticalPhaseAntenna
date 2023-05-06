simulation_data = load('test.mat');
phi0 = simulation_data.phi0;
d_phi = simulation_data.d_phi;
dphi = simulation_data.dphi;   
phase_error = simulation_data.phase_error;
lambda = simulation_data.lambda;
z = simulation_data.z;
varphi = simulation_data.varphi;
dd_phi = load('ga.mat').x;
step = simulation_data.step;
au = getfarfieldpattern(phi0,d_phi+dphi+dd_phi,phase_error,lambda,z,varphi);
[peakmainlobe, mainlobepeakpos] = max(au);

mainlobepos  = round(mainlobepeakpos-step/2):round(mainlobepeakpos+step/2);

mainlobe = au(mainlobepos);

psll = getPSLL(mainlobe);
figure(21);
plot(varphi*180/pi,au);
saveas(gcf,'ga_result1.fig');
figure(22);
plot((mainlobepos-length(varphi)/2)*step_phi*180/pi, mainlobe.^2/max(mainlobe.^2));
saveas(gcf,'ga_result2.fig');
