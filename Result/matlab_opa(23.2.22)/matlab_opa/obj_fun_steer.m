%% updading date: 2/14/2023
function y = obj_fun_steer(dd_phi_stg,steer)
    simulation_data = load('test.mat');
    num_antenna = simulation_data.num_antenna;
    xi = steer*pi/180;

    bestang = simulation_data.bestang;
    phase_error = simulation_data.phase_error;
    lambda = simulation_data.lambda;
    k = 2*pi/lambda;
    z = simulation_data.z;
    d = simulation_data.d;
    step_phi = simulation_data.step_phi;
    range  = -4.83*pi/180:step_phi:4.83*pi/180;
    varphi = range - xi;
    phi0 = zeros(length(varphi),num_antenna);
    for i = 1:num_antenna
        phi0(:,i) = (i-1)*(k*d*sin(range));
    end
    d = simulation_data.d;
    dd_phi(1:2:num_antenna) = dd_phi_stg;
    dd_phi(2:2:num_antenna) = dd_phi_stg;
    au = getfarfieldpattern(phi0,bestang(steer+5,:)+dd_phi,phase_error,lambda,z,varphi,d);
    psll = getPSLL(au);
    y = psll;
end

