function y = obj_fun(dd_phi)
    simulation_data = load('test.mat');
    phi0 = simulation_data.phi0;
    d_phi = simulation_data.d_phi;
    phase_error = simulation_data.phase_error;
    lambda = simulation_data.lambda;
    z = simulation_data.z;
    varphi = simulation_data.varphi;
    step_phi = simulation_data.step_phi;
    d = simulation_data.d;
%     [m,~] = size(phi0);
    au = getfarfieldpattern(phi0,d_phi+dd_phi,phase_error,lambda,z,varphi,d);
           
%     a = diff(abs(au));%%%
%     
%     threshold_diff = max(a)*0.3;
%     
%     for k = 1:m-1
%         if(abs(a(k))<threshold_diff)
%             a(k)=0;
%         end
%     end
%     
%     b = zeros(1,m-1);
%     
%     step = 1/(10*step_phi);
%     
%     [~, mainlobepeakpos] = max(au);
%     
%     mainlobepos  = round(mainlobepeakpos-step/2):round(mainlobepeakpos+step/2);
%     mainlobe = au(mainlobepos);

% %     step = 1700;
% %     mainlobepeakpos = length(au)/2;
% %     mainlobepos  = round(mainlobepeakpos-step/2):round(mainlobepeakpos+step/2);
% %     mainlobe = au(mainlobepos);
% %     psll = getPSLL(mainlobe);
    psll = getPSLL(au);
    
    y = psll;
end

