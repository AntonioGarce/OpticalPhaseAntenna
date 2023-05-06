tic
d_phi = load('test.mat').d_phi;
lb = -2*pi/3*ones(1,length(d_phi));
ub = -lb;
options = optimoptions('ga','PlotFcn',"gaplotbestf");
options = optimoptions(options,'MaxGenerations',100);
[x,fval]=ga(@obj_fun,length(d_phi),[],[],[],[],lb,ub,[],options); 
save ga.mat
toc

