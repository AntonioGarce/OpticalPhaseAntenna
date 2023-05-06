clear;
clc; close all;
num_antenna = 64;

z = 1.8;

lambda = 1.550e-6;
d = 9.5e-6;

k = 2*pi/lambda;
d_phi = zeros(1,num_antenna);

phase_error = 2*pi*rand(1,length(d_phi));

resfilname = "result_iter.xlsx";
xi = 0*pi/180;
step_phi = 0.0001; %(rad)

varphi = (-xi-pi/360:step_phi:-xi+pi/360);
phi0 = zeros(length(varphi),num_antenna);
xlswrite(resfilname,phase_error,"phase error");
for i = 1:num_antenna
    phi0(:,i) = (i-1)*(k*d*sin(varphi+xi));
end

[m,~] = size(phi0);
step_dphi = pi*2e-3;

tic
comp = zeros(num_antenna,log2(num_antenna));
min_psll = inf;
min_objval = inf;
for iter=1:2
    if iter == 2
        min_objval = min_psll;
    end
    for i = 1:log2(num_antenna)-1 %stage
        numgroup = pow2(i);
        numpergroup = num_antenna/numgroup;
    
        best_dphi = zeros(1,num_antenna);
        minpsll=zeros(10,1);
        for group = 1:numgroup %iterate subgroups
            if(mod(group,2)==1)
                groupnum = (numgroup-group+1)/2;    %subgroup number
            else
                groupnum = (numgroup+group)/2;      %subgroup number
            end
            subgroup = (numpergroup*(groupnum-1)+1):(numpergroup*groupnum);

            dphi = zeros(1,num_antenna);
            for delta_phi = 0:step_dphi:2*pi
                preb_dphi = best_dphi;
                preb_dphi(subgroup) = 0;
                dphi(subgroup) = delta_phi;  
                au = getfarfieldpattern(phi0,d_phi+preb_dphi+dphi,phase_error,lambda,z,varphi,d);
                mainlobe = au;
                psll = getPSLL(mainlobe);
                if(i==1 && iter==1 && group ==1 && delta_phi==0)
                    figure(1);
                    plot(varphi*180/pi,mainlobe.^2);
                    save("original_iter.mat", 'mainlobe');
                    saveas(gcf,"first_iter.fig");
                    firstpsll = psll;
                    xlswrite(resfilname,firstpsll,"first psll");
                end
                if (iter==1)
                    obj_val = psll - au(floor(length(au)/2))*10;
                else
                    obj_val = psll;
                end
                if obj_val<min_objval
                    min_objval = obj_val
                    min_psll = psll;
                    minpsll(iter)=psll;
                    best_dphi(subgroup) = dphi(subgroup);
                    best_mainlobe = mainlobe;
                end
            end
            if (iter>=2 && minpsll(iter)-minpsll(iter-1)>-0.1)
                break;
            end
        end
        comp_ang = zeros(num_antenna,1);
        for group=1:numgroup
            comp_ang(group) = best_dphi(numpergroup*(group-1)+1);
        end
        comp(:,i) = comp(:,i) + comp_ang;

        d_phi = d_phi + best_dphi;
    end
%     plot(varphi*180/pi, best_mainlobe.^2);
    
    for stg = 1:log2(num_antenna)-1
        sheet = ['stage_', num2str(stg),'iter_',num2str(iter)];
        comp_ang = comp(:,stg);
        comp_grp = comp_ang(1:2^stg);
        xlswrite(resfilname,comp_grp,sheet);
    end
    d_phi = d_phi-pi*floor(d_phi/pi);
    xlswrite(resfilname,d_phi,"compensate_angle");
    xlswrite(resfilname,min_psll,"final psll"+"iter_"+num2str(iter));
    figure(2+iter*10);
    plot(varphi*180/pi, best_mainlobe.^2);
    saveas(gcf,"last_iter"+num2str(iter)+".fig");
end




save test.mat
toc

