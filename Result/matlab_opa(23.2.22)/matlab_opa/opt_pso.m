tic
d_phi = load('test.mat').d_phi;
lb = -2*pi/3*ones(1,length(d_phi));
ub = -lb;
options = optimoptions('particleswarm','PlotFcn',"pswplotbestf",'SwarmSize',50);
options = optimoptions(options,'MaxIterations',500);
[x,fval]=particleswarm(@obj_fun,length(d_phi),lb,ub,options);

save pso.mat
toc

