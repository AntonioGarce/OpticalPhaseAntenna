%% updading date: 2/14/2023
clear;
clc; close all;
num_antenna = 64;

z = 5;

lambda = 1.550e-6;
d = 9.5e-6;
% w = 0;
k = 2*pi/lambda;

phase_error = load('test.mat').phase_error;
% phase_error = 1.4*pi*rand(1,num_antenna)-pi*0.7;
% phase_error = 0;
step_phi = 0.0006; %(rad)
range = -4.83*pi/180:step_phi:4.83*pi/180;
reduced_range = -4.83*pi/180:step_phi:4.83*pi/180;
for i = 1:num_antenna
    ph0(:,i) = (i-1)*(k*d*sin(range));
end
bestu = zeros(9,length(reduced_range));
au = getfarfieldpattern(ph0,0,0,lambda,z,range,d);
orgpsll = getPSLL(au);
maxau = max(au);
writematrix(au,"result_origin.xlsx","Sheet","originfarfield");
writematrix(orgpsll,"result_origin.xlsx","Sheet","psll");


w=100;
bestang = zeros(9,num_antenna);
delete(gcp('nocreate'));
% parpool(9)
tic
parfor steer = -4:4

    minpsll = inf;
    xi = steer*pi/180;         
    varphi = reduced_range-xi;
    phi0 = zeros(length(varphi),num_antenna);
    for i = 1:num_antenna
        phi0(:,i) = (i-1)*(k*d*sin(varphi+xi));
    end
    
    [m,~] = size(phi0);
    comp = zeros(num_antenna,log2(num_antenna));
    comp_total = pi*ones(1,num_antenna);
   
    au_org = getfarfieldpattern(phi0,0,phase_error,lambda,z,varphi,d);
    psll = getPSLL(au_org);
    figure(1);
    plot(varphi*180/pi,(au_org/maxau).^2);
    saveas(gcf,"first("+num2str(steer)+").fig");
    firstpsll = psll;
    writematrix(firstpsll,"result(steer"+num2str(steer)+").xlsx","Sheet","first psll");
    resfilname = "result(steer"+num2str(steer)+").xlsx";
    writematrix(phase_error,resfilname,"Sheet","phase error");

    %% initialization for 1-dimension serach optimization
    epsilon = 0.0001;
    alpha = 2;
    %% optimization
    num_stage = log2(num_antenna)-1;

    for stg = 1:num_stage %stage
        numgroup = pow2(stg);
        numpergroup = num_antenna/numgroup;
        y = comp_total;
        comp_start = comp_total;
        delta = 0.05;
        
        while(delta>epsilon)
            au_xprev = getfarfieldpattern(phi0,comp_total,phase_error,lambda,z,varphi,d);
            psll = getPSLL(au_xprev);
            obj_val_xprev = psll - au_xprev(floor(length(au_xprev)/2))*w;

            au_prev = getfarfieldpattern(phi0,y,phase_error,lambda,z,varphi,d);
            psll = getPSLL(au_prev);
            
            obj_val_prev = psll - au_prev(floor(length(au_prev)/2))*w;

            for i = 1:numgroup %iterate subgroups
                if(mod(i,2)==1)
                    groupnum = (numgroup-i+1)/2;    %subgroup number
                else
                    groupnum = (numgroup+i)/2;      %subgroup number
                end
                subgroup = (numpergroup*(groupnum-1)+1):(numpergroup*groupnum);
                Di = zeros(1,num_antenna);
                Di(subgroup) = 1;
                y_p = y + delta*Di;
                
                au = getfarfieldpattern(phi0,y_p,phase_error,lambda,z,varphi,d);
                psll = getPSLL(au);
                obj_val_p = psll - au(floor(length(au)/2))*w;

                if obj_val_p < obj_val_prev
                    y = y_p;
                    obj_val_prev = obj_val_p;
                else 
                    y_n = y -delta*Di;
                    au = getfarfieldpattern(phi0,y_n,phase_error,lambda,z,varphi,d);
                    psll = getPSLL(au);
                    obj_val_n = psll - au(floor(length(au)/2))*w;
                    if obj_val_n < obj_val_prev
                        obj_val_prev = obj_val_n;
                        y = y_n;
                    end
                end
            end
            if obj_val_xprev>obj_val_prev
                comp_total_prev = comp_total;
                comp_total = y;
                y = comp_total + alpha*(comp_total - comp_total_prev);
            else
                delta = delta/2;
                y = comp_total;
            end
        end
        comp(:,stg) = y - comp_start;
    end
        
    for stg = 1:log2(num_antenna)-1
        sheet = ['stage_', num2str(stg)];
        comp_ang = comp(:,stg);
        comp_grp = comp_ang(1:2^stg);
        writematrix(comp_grp,resfilname,"Sheet",sheet);
    end
    
    writematrix(comp_total,resfilname,"Sheet","compensate_angle");
    
    au = getfarfieldpattern(phi0,comp_total,phase_error,lambda,z,varphi,d);

    comp_total_prev=comp_total;

    bestang(steer+5,:) = comp_total;

    prev_psll = psll;
    psll = getPSLL(au);
    
    writematrix(psll,resfilname,"Sheet","final psll");
    title("result with steer with "+num2str(steer));
    figure(2);
    plot(varphi*180/pi,(au/maxau).^2);
    saveas(gcf,"last(steer"+num2str(steer)+").fig"');
end

delete(gcp('nocreate'));
for steer_i = -4:4
    xi_ = steer_i*pi/180;
    varphi_ = range-xi_;
    for i = 1:num_antenna
        ph0(:,i) = (i-1)*(k*d*sin(varphi_+xi_));
    end
    best_ang = bestang(steer_i+5,:);
    best_u = getfarfieldpattern(ph0,best_ang,phase_error,lambda,z,varphi_,d);
    figure(3);
    plot(varphi_*180/pi, (best_u/maxau).^2);
    legend("psll="+num2str(getPSLL(best_u))+"db");
    hold on
    saveas(gcf,"last_newgrouping_origin.fig");
end

toc

