m = 3
sol = solver('x0',x0,'lbx',lbx,'ubx',ubx,'lbg',lbg,'ubg',ubg,'p',[param(2,:)';0;0]);

figure(3);clf(3)
e = u(:,m:m+1)-full(sol.x)';
ne = mean(e)
tmp = [0 0; 0 0];
for i = 1:5000
	tmp = tmp + (e(i,:)' - ne')*(e(i,:)' - ne')';
	i	
end
tmp = (tmp);

t = 0:0.01:2*pi;
u1 = cos(t);
u2 = sin(t);
utmp = [u1;u2];
y = tmp*utmp + [mean(u(:,m));mean(u(:,m+1))];
figure(3);
plot(y(1,:),y(2,:)); hold on;
scatter(u(:,m),u(:,m+1),4,'.')
axis equal