% clear;

U 		= casadi.SX.sym('U',2,1);
d 		= casadi.SX.sym('d',2,1);
theta 	= casadi.SX.sym('theta',3,1);


obj = 3*norm(U - [theta(1);theta(2)])^2 + d'*U;
g 	= norm(U)^2 - theta(3);

NLP = struct('f',obj,'x',U,'g',g,'p',[theta;d]);
lbg = -inf(1,1);
ubg = zeros(1,1);
lbx = -inf(2,1);
ubx = inf(2,1);
x0 = zeros(2,1);

opts = struct;
opts.ipopt.max_iter = 100;
opts.ipopt.print_level = 0;
opts.print_time = 0;
opts.ipopt.acceptable_tol = 1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-8;
opts.ipopt.mu_target= 10^-8;

solver = casadi.nlpsol('solver','ipopt',NLP,opts);



sigma = 0.5;
mu = 0;

rng = 5000;
u = NaN(rng,6);
param = [0.5,0.5,1;0.64,0.64,1;0.7071,0.7071,1];


for i = 1:rng
	d = sigma.*randn(2,1) + mu;
	for j = 1:3
		sol = solver('x0',x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p',[param(j,:)';d]);
		u(i,(j-1)*2 + 1:j*2) = full(sol.x);
		
		disp(i);
	
	end
	

end
sol = solver('x0',x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p',[param(1,:)';0;0]);
uopt = full(sol.x);
sol = solver('x0',x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p',[param(2,:)';0;0]);
uopt = [uopt ;full(sol.x)];
sol = solver('x0',x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p',[param(3,:)';0;0]);
uopt = [uopt ;full(sol.x)];

n = 16;
n2 = 20;
figure(3); clf(3);
scatter(  u(:,1),      u(:,2),n,'b.'); grid on; hold on;
plot(mean(u(:,1)),mean(u(:,2)),'kx','Markersize',n2,'LineWidth',2);
sol = solver('x0',x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p',[param(1,:)';0;0]);
plot(full(sol.x(1)),full(sol.x(2)),'k+','Markersize',n2,'LineWidth',1.5);

scatter(  u(:,3),      u(:,4),n,'r.');
plot(mean(u(:,3)),mean(u(:,4)),'kx','Markersize',n2,'LineWidth',2)
sol = solver('x0',x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p',[param(2,:)';0;0]);
plot(full(sol.x(1)),full(sol.x(2)),'k+','Markersize',n2,'LineWidth',1.5);

scatter(  u(:,5),      u(:,6),n,'g.');
plot(mean(u(:,5)),mean(u(:,6)),'kx','Markersize',n2,'LineWidth',2)
sol = solver('x0',x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p',[param(3,:)';0;0]);
plot(full(sol.x(1)),full(sol.x(2)),'k+','Markersize',n2,'LineWidth',1.5);
u1 = 0:0.01:1;
plot(u1,sqrt(param(1,3) - u1.^2),'k');

for i = 1:3
    n1 = (i-1)*2 + 1; n2 = n1 + 1
    A = cov(u(:,n1),u(:,n2))*10;
    t = 0:0.01:2*pi;
    u1 = cos(t);
    u2 = sin(t);
    utmp = [u1;u2];
    y = A*utmp + [mean(u(:,n1));mean(u(:,n2))];
    figure(3);
    plot(y(1,:),y(2,:),'k')
    axis equal
end

