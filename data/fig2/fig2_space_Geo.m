function [u, ds, ns, Mat] = fig2_space_Geo(dt,nt)
    geo_space = load('data/fig2/fig2_x.mat');
    x = geo_space.x;
    ds = geo_space.ds;
    ns = geo_space.ns;
    Mat = geo_space.Mat;
    u = zeros(size(x,1), size(x,2), nt);
	lx = max(x(:,1));
	ly = max(x(:,2));
	w = lx/15;
	A = 1;
	for ti = 1:nt
		f  = [lx/4,ly/2];
		r1 = [+1, 0];
		g1 = exp(-(x(:,1)-f(1)).^2/(2*w^2)-(x(:,2)-f(2)).^2/(2*w^2))*A.*r1;

		f  = [lx/4+lx/2,ly/2];
		r2 = [-1, 0];
		g2 = exp(-(x(:,1)-f(1)).^2/(2*w^2)-(x(:,2)-f(2)).^2/(2*w^2))*A.*r2;

		f  = [lx/2,ly/4];
		r3 = [0, +1];
		g3 = exp(-(x(:,1)-f(1)).^2/(2*w^2)-(x(:,2)-f(2)).^2/(2*w^2))*A.*r3;

		f  = [lx/2,ly/4+ly/2];
		r4 = [0, -1];
		g4 = exp(-(x(:,1)-f(1)).^2/(2*w^2)-(x(:,2)-f(2)).^2/(2*w^2))*A.*r4;

		u(:,[1,2],ti)= ds(1)*(g1+g2+g3+g4)/6;
	end
end