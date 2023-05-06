clear;
clc; close all;
num_antenna = 64;

z = 3;

lambda = 1.550e-6;
d = 9.5e-6;
% w = 0;
k = 2*pi/lambda;


phase_error = 1.8*pi*rand(1,num_antenna)-pi*0.9;

range = -4.83*pi/180:0.0001:4.83*pi/180;
for i = 1:num_antenna
    ph0(:,i) = (i-1)*(k*d*sin(range));
end
au = getfarfieldpattern(ph0,0,0,lambda,z,range,d);
orgpsll = getPSLL(au);
maxau = max(au);
delete(gcp('nocreate'));
parpool(9)
tic
for w=0
    parfor steer = -4:4
%         bestu = zero
        d_phi = zeros(1,num_antenna);
        resfilname = "result(w"+num2str(w)+"steer"+num2str(steer)+").xlsx";
        xi = steer*pi/180;
        step_phi = 0.0006; %(rad)
%         
        varphi = (-xi-4.83*pi/180:step_phi:-xi+4.83*pi/180);
        phi0 = zeros(length(varphi),num_antenna);
        writematrix(phase_error,resfilname,"Sheet","phase error");
        for i = 1:num_antenna
            phi0(:,i) = (i-1)*(k*d*sin(varphi+xi));
        end
        
        [m,~] = size(phi0);
        step_dphi = pi*2e-3;
        

        comp = zeros(num_antenna,log2(num_antenna));
    
        for stg = 1:log2(num_antenna)-1 %stage
            numgroup = pow2(stg);
            numpergroup = num_antenna/numgroup;
        
            best_dphi = zeros(1,num_antenna);
            for group = 1:numgroup %iterate subgroups
                if(mod(group,2)==1)
                    groupnum = (numgroup-group+1)/2;    %subgroup number
                else
                    groupnum = (numgroup+group)/2;      %subgroup number
                end
                subgroup = (numpergroup*(groupnum-1)+1):(numpergroup*groupnum);
                min_psll = inf;
                min_objval = inf;
                dphi = zeros(1,num_antenna);
                for delta_phi = 0:step_dphi:2*pi
                    preb_dphi = best_dphi;
                    preb_dphi(subgroup) = 0;
                    dphi(subgroup) = delta_phi;  
                    au = getfarfieldpattern(phi0,d_phi+preb_dphi+dphi,phase_error,lambda,z,varphi,d);
                    mainlobe = au;
                    psll = getPSLL(mainlobe);
                    if(stg==1  && group ==1 && delta_phi==0)
                        figure(1);
                        plot(varphi*180/pi,mainlobe.^2);
                        saveas(gcf,"first("+num2str(steer)+").fig");
                        firstpsll = psll;
                        writematrix(firstpsll,resfilname,"Sheet","first psll");
                    end
                    obj_val = psll - au(floor(length(au)/2))*w;
                    if obj_val<min_objval
                        min_objval = obj_val;
                        min_psll = psll;
                        best_dphi(subgroup) = dphi(subgroup);
                        best_mainlobe = mainlobe;
                    end
                end
            end
            comp_ang = zeros(num_antenna,1);
            for group=1:numgroup
                comp_ang(group) = best_dphi(numpergroup*(group-1)+1);
            end
            comp(:,stg) = comp(:,stg) + comp_ang;
    
            d_phi = d_phi + best_dphi;
        end

        
        for stg = 1:log2(num_antenna)-1
            sheet = ['stage_', num2str(stg)];
            comp_ang = comp(:,stg);
            comp_grp = comp_ang(1:2^stg);
            writematrix(comp_grp,resfilname,"Sheet",sheet);
        end
        d_phi = d_phi-pi*floor(d_phi/pi);
        writematrix(d_phi,resfilname,"Sheet","compensate_angle");
        writematrix(min_psll,resfilname,"Sheet","final psll");

        figure(2+w);
        plot(varphi*180/pi, (best_mainlobe/maxau).^2);
        hold on
        saveas(gcf,"last(w"+num2str(w)+").fig"');
    end
end
delete(gcp('nocreate'));
save test.mat
toc

