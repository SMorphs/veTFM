clc; close all; clear
addpath('src');

%% FIG. 2 - EXACT
addpath("data/fig2/");
h = 2e-6;
dts = logspace(-1,log10(20),21);
time    = 40*1;
for dti = 1:length(dts)
    dt = dts(dti);
    nt = round(time/dt);
    [u, ds, ns, Mat] = fig2_space_Geo(dt,nt);
    t      = viscoI_s(u,ds,ns,dt,h,32,Mat);
    u_star = viscoF_s(t,ds,ns,dt,h,32,Mat);
end

for dti = 1:length(dts)
    dt = dts(dti);
    nt = round(time/dt);
    [u, ds, ns, Mat] = fig2_time_Geo(dt,nt);
    t      = viscoI_s(u,ds,ns,dt,h,32,Mat);
    u_star = viscoF_s(t,ds,ns,dt,h,32,Mat);
end

%% FIG. 3
addpath("data/fig3/");
h = 1.0775e-04;

% CM
cm_vetfm = load('data/fig3/cm_vetfm.mat');
veTFMMat = cm_vetfm.Mat;
eTFM0Mat = veTFMMat; eTFM0Mat.E = veTFMMat.E+veTFMMat.c(1,1)*(1+veTFMMat.nu)+veTFMMat.c(2,1)*(1+veTFMMat.nu);
eTFMIMat = veTFMMat;
t_cm_veTFM_Num = viscoI_s(cm_vetfm.u,cm_vetfm.ds,cm_vetfm.ns,cm_vetfm.dt,h,32,veTFMMat); % veTFM - Numerical inversion
t_cm_veTFM_Ana = viscoI_t(cm_vetfm.u,cm_vetfm.ds,cm_vetfm.ns,cm_vetfm.dt,h,veTFMMat); % veTFM - Analytical inversion
t_cm_eTFM0 = elastI_t(cm_vetfm.u,cm_vetfm.ds,cm_vetfm.ns,h,eTFM0Mat); % eTFM0
t_cm_eTFMI = elastI_t(cm_vetfm.u,cm_vetfm.ds,cm_vetfm.ns,h,eTFMIMat); % eTFMI

% MCF10a
mcf10a_vetfm = load('data/fig3/mcf10a_vetfm.mat');
veTFMMat = mcf10a_vetfm.Mat;
eTFM0Mat = veTFMMat; eTFM0Mat.E = veTFMMat.E+veTFMMat.c(1,1)*(1+veTFMMat.nu)+veTFMMat.c(2,1)*(1+veTFMMat.nu);
eTFMIMat = veTFMMat;
t_mcf10a_veTFM_Num = viscoI_s(mcf10a_vetfm.u,mcf10a_vetfm.ds,mcf10a_vetfm.ns,mcf10a_vetfm.dt,h,32,veTFMMat); % veTFM - Numerical inversion
t_mcf10a_veTFM_Ana = viscoI_t(mcf10a_vetfm.u,mcf10a_vetfm.ds,mcf10a_vetfm.ns,mcf10a_vetfm.dt,h,veTFMMat); % veTFM - Analytical inversion
t_mcf10a_eTFM0 = elastI_t(mcf10a_vetfm.u,mcf10a_vetfm.ds,mcf10a_vetfm.ns,h,eTFM0Mat); % eTFM0
t_mcf10a_eTFMI = elastI_t(mcf10a_vetfm.u,mcf10a_vetfm.ds,mcf10a_vetfm.ns,h,eTFMIMat); % eTFMI

% HDF
hdf_vetfm = load('data/fig3/hdf_vetfm.mat');
veTFMMat = hdf_vetfm.Mat;
eTFM0Mat = veTFMMat; eTFM0Mat.E = veTFMMat.E+veTFMMat.c(1,1)*(1+veTFMMat.nu)+veTFMMat.c(2,1)*(1+veTFMMat.nu);
eTFMIMat = veTFMMat;
t_hdf_veTFM_Num = viscoI_s(hdf_vetfm.u,hdf_vetfm.ds,hdf_vetfm.ns,hdf_vetfm.dt,h,32,veTFMMat); % veTFM - Numerical inversion
t_hdf_veTFM_Ana = viscoI_t(hdf_vetfm.u,hdf_vetfm.ds,hdf_vetfm.ns,hdf_vetfm.dt,h,veTFMMat); % veTFM - Analytical inversion
t_hdf_eTFM0 = elastI_t(hdf_vetfm.u,hdf_vetfm.ds,hdf_vetfm.ns,h,eTFM0Mat); % eTFM0
t_hdf_eTFMI = elastI_t(hdf_vetfm.u,hdf_vetfm.ds,hdf_vetfm.ns,h,eTFMIMat); % eTFMI

%% FIG. 4 - CLOSE ENOUGH
addpath("data/fig3/"); addpath("data/fig4/");
h = 1.0775e-04;
alphas1 = linspace(0.05,0.95,10);
alphas2 = linspace(0.05,0.95,10);
taus1 = logspace(log10(0.5),log10(30),10);
taus2 = logspace(log10(0.5),log10(30),10);
Mat0 = baseMat();
Einf0 = Mat0.E/(Mat0.nu+1);
E1   = Mat0.c(1,1);
eta1 = Mat0.c(1,2);
E2   = Mat0.c(2,1);
eta2 = Mat0.c(2,2);
tau1_0 = eta1/E1;
tau2_0 = eta2/E2;
a1_0   = E1/(Einf0+E1+E2);
a2_0   = E2/(Einf0+E1+E2);
alphas1 = alphas1([1,4]);
alphas2 = alphas2([1,4]);
common_tFree = 5*max([taus1,taus2]);
common_tN = 6*max([taus1,taus2]);
for ai = 1:length(alphas1)
    for aj = 1:length(alphas2)
	    Mat = Mat0;
	    a1 = alphas1(ai);
	    a2 = alphas2(aj);
        t1 = tau1_0;
	    t2 = tau2_0;
	    Ei = a1*Einf0./(1-a1-a2);
	    Ej = a2*Einf0./(1-a1-a2);
	    Mat.c(1,1) = Ei;
	    Mat.c(2,1) = Ej;
	    Mat.c(1,2) = t1*Ei;
	    Mat.c(2,2) = t2*Ej;
        MatH = Mat;
        MatH.E = Mat.E+Mat.c(1,1)*(1+Mat.nu)+Mat.c(2,1)*(1+Mat.nu);
        % CM field
        u_vetfm = load('data/fig3/cm_vetfm.mat');
        % MCF10a field
        % u_vetfm = load('data/fig3/mcf10a_vetfm.mat');
        % HDF field
        % u_vetfm = load('data/fig3/hdf_vetfm.mat');
        time = (0:(size(u_vetfm.u,3)-1))*u_vetfm.dt;
        endidx = find(time>=common_tN,1);
        time = time(1:endidx);
        tAnalysisIdx = time > common_tFree;
        t_ij_veTFM = viscoI_s(u_vetfm.u,u_vetfm.ds,u_vetfm.ns,u_vetfm.dt,h,32,Mat);
        t_ij_eTFM0 = elastI_t(u_vetfm.u,u_vetfm.ds,u_vetfm.ns,h,MatH);
        t_ij_eTFMI = elastI_t(u_vetfm.u,u_vetfm.ds,u_vetfm.ns,h,Mat);
    end
end

