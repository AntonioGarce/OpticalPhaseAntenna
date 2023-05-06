%% updading date: 2/27/2023(steering angle -9 to 9 + distribution range + individual figs + phase comp)
clear;
clc; close all;
num_antenna = 64;

z = 5;

ymin = -pi;
ymax = pi;
ycm = 0;

lambda = 1.550e-6;
d = 9.5e-6;
% w = 0;
k = 2*pi/lambda;

% phase_error = load("test_0.mat").phase_error;
phase_error = 1.4*pi*rand(1,num_antenna)-pi*0.7;
% phase_error = 0;
step_phi = 0.0006; %(rad)
range = -4.83*pi/180:step_phi:4.83*pi/180;
reduced_range = -4.83*pi/180:step_phi:4.83*pi/180;
ph0 = zeros(length(range),num_antenna);
for i = 1:num_antenna
    ph0(:,i) = (i-1)*(k*d*sin(range));
end
bestu = zeros(9,length(reduced_range));
au = getfarfieldpattern(ph0,0,0,lambda,z,range,d);
orgpsll = getPSLL(au);
maxau = max(au);
writematrix(au,"result_origin.xlsx","Sheet","originfarfield");
writematrix(orgpsll,"result_origin.xlsx","Sheet","psll");

tic
w=100;
bestang = zeros(9,num_antenna);

for steer = 0
    minpsll = inf;
    xi = steer*pi/180;         
    varphi = reduced_range-xi;
    phi0 = zeros(length(varphi),num_antenna);
    for i = 1:num_antenna
        phi0(:,i) = (i-1)*(k*d*sin(varphi+xi));
    end
    
    [m,~] = size(phi0);
    comp = zeros(num_antenna,log2(num_antenna));
    best_phase = zeros(num_antenna,log2(num_antenna));
    comp_total = ycm*ones(1,num_antenna);
   
    au_org = getfarfieldpattern(phi0,0,phase_error,lambda,z,varphi,d);
    psll = getPSLL(au_org);
    figure(1);
    plot(varphi*180/pi,(au_org/maxau).^2);
    saveas(gcf,"first("+num2str(steer)+").fig");
    firstpsll = psll;
    writematrix(firstpsll,"result(steer"+num2str(steer)+").xlsx","Sheet","first psll");
    w1=1;
    w0 = 100;
    wstep = 20;

    prev_psll = 100;
    itenum = 10;
    w0num = 0;
    w0flag = 0;
    for iter =1:10
        w = w0-wstep*iter; 
        if w<0 
            w=0;
        end
        resfilname = "result(steer"+num2str(steer)+"iter"+num2str(iter)+").xlsx";
        writematrix(phase_error,resfilname,"Sheet","phase error");
        %% initialization for optimization
        epsilon = 0.0001;
        alpha = 2;
        %% optimization
        num_stage = log2(num_antenna)-1;

        for stg = 1:num_stage %stage
            numgroup = pow2(stg);
            numpergroup = num_antenna/numgroup;
%                 comp_total_prev = comp_total;
            y = comp_total;
            delta = 0.05;
            comp_start = comp_total;
            while(delta>epsilon)
                au_xprev = getfarfieldpattern(phi0,comp_total,phase_error,lambda,z,varphi,d);
                psll = getPSLL(au_xprev);
                obj_val_xprev = psll - au_xprev(floor(length(au_xprev)/2))*w*w1;

                au_prev = getfarfieldpattern(phi0,y,phase_error,lambda,z,varphi,d);
                psll = getPSLL(au_prev);
                
                obj_val_prev = psll - au_prev(floor(length(au_prev)/2))*w*w1;

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
                    obj_val_p = psll - au(floor(length(au)/2))*w*w1;

                    if obj_val_p < obj_val_prev
                        y = y_p;
                        obj_val_prev = obj_val_p;
                    else 
                        y_n = y -delta*Di;
                        au = getfarfieldpattern(phi0,y_n,phase_error,lambda,z,varphi,d);
                        psll = getPSLL(au);
                        obj_val_n = psll - au(floor(length(au)/2))*w*w1;
                        if obj_val_n < obj_val_prev
                            obj_val_prev = obj_val_n;
                            y = y_n;
                        end
                    end
                end
                if (obj_val_xprev>obj_val_prev) && (sum(y<ymin)==0) && (sum(y>ymax)==0)
                    comp_total_prev = comp_total;
                    comp_total = y;
                    y = comp_total + alpha*(comp_total - comp_total_prev);
                else
                    delta = delta/2;
                    y = comp_total;
                end
            end
            comp(:,stg) = comp(:,stg) + comp_total' - comp_start';
            best_phase(:,stg) = comp_total;
        end
        
        writematrix(comp_total,resfilname,"Sheet","compensate_angle");
        
        au = getfarfieldpattern(phi0,comp_total,phase_error,lambda,z,varphi,d);
        if  psll>=prev_psll || iter>=8 
            w0flag = 1;
        end

        comp_total_prev=comp_total;
        if w0flag ==1
            w0num = w0num+1;
        end
        if minpsll>psll
            minpsll = psll;
            bestang(steer+5,:) = comp_total;
        end
        prev_psll = psll;
        psll = getPSLL(au);
        
        writematrix(psll,resfilname,"Sheet","final psll");
%         figure(iter+10);
%         plot(varphi*180/pi,(au/maxau).^2);
        if w0num ==3 
            break;
        end
    end
    best_comp = 0;
    for stg = 1:log2(num_antenna)-1
        sheet = ['stage_', num2str(stg)];
        comp_ang = comp(:,stg);
        best_comp = best_comp+comp_ang;
        grp = 1:2^(6-stg):64;
        comp_grp = comp_ang(grp);
        writematrix(comp_grp,resfilname,"Sheet",sheet);
        figure(20+stg+iter*100);
        subplot(2,1,1);
        plot(comp_grp,'o-');
        title("phase correction at stage "+num2str(stg)+" of iteration "+num2str(iter));
        au = getfarfieldpattern(phi0,best_comp',phase_error,lambda,z,varphi,d);
        psll = getPSLL(au);
        subplot(2,1,2);
        plot(reduced_range,(au/maxau).^2);
        legend("psll = "+num2str(psll));
        title("farfiled pattern at stage "+num2str(stg)+" of iteration "+num2str(iter));
    end
end

toc

