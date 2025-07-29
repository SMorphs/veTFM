function Mat = baseMat()
    Mat.model = 'rk4zener';
    Mat.elast = 'hookean';
    Mat.nu = 0.4;
    Mat.E =   500*(1+Mat.nu);
    Mat.c = [ 500,  1*500; 4000 4000*20];
	Mat.visco = 'dev';
    Mat.ve_elems = 2;
end