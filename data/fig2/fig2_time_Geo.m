function [u, ds, ns, Mat] = fig2_time_Geo(dt,nt)
    f   = 1/10;
    w   = 2*pi*f;
    geo_space = load('data/fig2/fig2_x.mat');
    x = geo_space.x;
    ds = geo_space.ds;
    ns = geo_space.ns;
    Mat = geo_space.Mat;
    u = zeros(size(x,1), size(x,2), nt);
	for ti = 1:nt
		t = (ti-1)*dt;
		u(:,1,ti) = ds(1)*0.0156*sin(t*w)*10;
	end
end