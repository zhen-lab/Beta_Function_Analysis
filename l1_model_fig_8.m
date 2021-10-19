function l1_model_fig_8()
    %wild type model
    ein_link = 1;
    iin_link = 1;
    [t,x] = l1_model(ein_link,iin_link);
    subplot(3,1,1);    
    imagesc(t,1:6,x(:,5:5:30)');
    xlim([1400 1800]);
    set(gca,'xticklabels',...
        cellfun(@(x) num2str(str2num(x)-1400),get(gca,'Xticklabels'),...
            'UniformOutput',false));
    colorbar;caxis([-6 6]);
    title('Model Wildtype (No Bias)');
   
    %model without inhibition
    ein_link = 1;
    iin_link = 0;
    [t,x] = l1_model(ein_link,iin_link);
    subplot(3,1,2);    
    imagesc(t,1:6,x(:,5:5:30)');
    xlim([1400 1800]);
    set(gca,'xticklabels',...
        cellfun(@(x) num2str(str2num(x)-1400),get(gca,'Xticklabels'),...
            'UniformOutput',false));
    colorbar;caxis([-6 6]);
    title('Model DD Ablation (Ventral Bias)');
    
    %model without extrasynaptic excitation
    ein_link = 0;
    iin_link = 1;
    [t,x] = l1_model(ein_link,iin_link);
    subplot(3,1,3);    
    imagesc(t,1:6,x(:,5:5:30)');
    xlim([1400 1800]);
    set(gca,'xticklabels',...
        cellfun(@(x) num2str(str2num(x)-1400),get(gca,'Xticklabels'),...
            'UniformOutput',false));
    colorbar;caxis([-6 6]);
    title('Model Extrasynaptic Input Ablation (Dorsal Bias)');
end
function [t,x] = l1_model(ein_link,iin_link)
​
    t=0:0.001:1800;   % time stamp
    
    initial_x    = [-30;0.3;-1;1;-1; ...
                    0;0;0;0;0; ...
                    0;0;0;0;0; ...
                    0;0;0;0;0; ...
                    0;0;0;0;0; ...
                    0;0;0;0;0; ...
                    ];
   
    [t,x]=ode45( @(t,x) model_eqns(t,x,ein_link,iin_link), t, initial_x );
    
end
​
%PNAS
function dxdt=model_eqns(t,x,ein_link,iin_link)
​
    %membrane capacitance
    Cm = 3;
    
    %conductances
    gl = 100; %leak
    gca = 400; %calcium
    gk = 500; %potassium 
    %g = 100;
    
    %reversal potentials
    El = -60; %leak
    Eca = 60; %calcium
    Ek = -70; %potassium
    Eplus = -10; %excitatory
    Estar = 33.5; %extrasynaptic
    Eminus = -67; %inhibitory
    
    
    vnn = -33.5;   %synaptic conductivity parameter
    vavb = -40; %avb input membrane voltge
    
    %timescales
    tau_n = 30;
    tau_u = 85;
    tau_b = 10;
    
    %proprioceptive coupling
    c = 5;
    
    %segment 1
    vd1 = x(1);
    nd1 = x(2);    
    md1 = x(3);
    mv1 = x(4);
    kap1 = x(5);
    
    d_vd_dt1 = (1/Cm)*( -gl*(vd1 - El) ...
                - gca*minf(vd1)*(vd1 - Eca) ...
                - gk*nd1*(vd1 - Ek) ...
                - syncon(vavb,100,-40,30)*(vd1 - Eplus )...
                   ); %dorsal  MN
    d_nd_dt1 = (1/tau_n)*(-nd1 + ninf(vd1));
    
    %dorsal muscle
    d_md_dt1 = (1/tau_u)*( -gl*(md1 - El ) - syncon(vd1,1000,-30,20)*(md1 - Eplus) ); 
    %ventral muscle
    d_mv_dt1 = (1/tau_u)*( -gl*(mv1 - El ) - iin_link*syncon(vd1,1000,vnn,20)*(mv1 - Eminus) ...
                           - ein_link*syncon(vavb,2000,-30,20)*(mv1 - Estar) );
    %curvature
    d_kap_dt1 = (1/tau_b)*( -kap1 + sigma_def(mus(md1)) - sigma_def(mus(mv1)) );
    
    %segment 2
    vd2 = x(6);
    nd2 = x(7);    
    md2 = x(8);
    mv2 = x(9);
    kap2 = x(10);
    
    d_vd_dt2 = (1/Cm)*( -gl*(vd2 - El) ...
                - gca*minf(vd2)*(vd2 - Eca) ...
                - gk*nd2*(vd2 - Ek) ...
                - syncon(vavb,100,-40,30)*(vd2 - Eplus )...
                + c*kap1  ); %dorsal  MN
    d_nd_dt2 = (1/tau_n)*(-nd2 + ninf(vd2));
    d_md_dt2 = (1/tau_u)*( -gl*(md2 - El ) - syncon(vd2,1000,-30,20)*(md2 - Eplus) ); %dorsal muscle
    d_mv_dt2 = (1/tau_u)*( -gl*(mv2 - El ) - iin_link*syncon(vd2,1000,vnn,20)*(mv2 - Eminus) ...
                           - ein_link*syncon(vavb,1000,-30,20)*(mv2 - Estar) ); %ventral muscle
    d_kap_dt2 = (1/tau_b)*( -kap2 + sigma_def(mus(md2)) - sigma_def(mus(mv2)) ); %curvature
   
    %segment 3
    vd3 = x(11);
    nd3 = x(12);    
    md3 = x(13);
    mv3 = x(14);
    kap3 = x(15);
    
    d_vd_dt3 = (1/Cm)*( -gl*(vd3 - El) ...
                - gca*minf(vd3)*(vd3 - Eca) ...
                - gk*nd3*(vd3 - Ek) ...
                - syncon(vavb,100,-40,30)*(vd3 - Eplus )...
                + c*kap2  ); %dorsal  MN
    d_nd_dt3 = (1/tau_n)*(-nd3 + ninf(vd3));
    d_md_dt3 = (1/tau_u)*( -gl*(md3 - El ) - syncon(vd3,1000,-30,20)*(md3 - Eplus) ); %dorsal muscle
    d_mv_dt3 = (1/tau_u)*( -gl*(mv3 - El ) - iin_link*syncon(vd3,1000,vnn,20)*(mv3 - Eminus) ...
                           - ein_link*syncon(vavb,1000,-30,20)*(mv3 - Estar) ); %ventral muscle
    d_kap_dt3 = (1/tau_b)*( -kap3 + sigma_def(mus(md3)) - sigma_def(mus(mv3)) ); %curvature
    
    %segment 4
    vd4 = x(16);
    nd4 = x(17);    
    md4 = x(18);
    mv4 = x(19);
    kap4 = x(20);
    
    d_vd_dt4 = (1/Cm)*( -gl*(vd4 - El) ...
                - gca*minf(vd4)*(vd4 - Eca) ...
                - gk*nd4*(vd4 - Ek) ...
                - syncon(vavb,100,-40,30)*(vd4 - Eplus )...
                + c*kap3  ); %dorsal  MN
    d_nd_dt4 = (1/tau_n)*(-nd4 + ninf(vd4));
    d_md_dt4 = (1/tau_u)*( -gl*(md4 - El ) - syncon(vd4,1000,-30,20)*(md4 - Eplus) ); %dorsal muscle
    d_mv_dt4 = (1/tau_u)*( -gl*(mv4 - El ) - iin_link*syncon(vd4,1000,vnn,20)*(mv4 - Eminus) ...
                           - ein_link*syncon(vavb,1000,-30,20)*(mv4 - Estar) ); %ventral muscle
    d_kap_dt4 = (1/tau_b)*( -kap4 + sigma_def(mus(md4)) - sigma_def(mus(mv4))); %curvature
    
    %segmeng 5
    vd5 = x(21);
    nd5 = x(22);    
    md5 = x(23);
    mv5 = x(24);
    kap5 = x(25);
    
    d_vd_dt5 = (1/Cm)*( -gl*(vd5 - El) ...
                - gca*minf(vd5)*(vd5 - Eca) ...
                - gk*nd5*(vd5 - Ek) ...
                - syncon(vavb,100,-40,30)*(vd5 - Eplus )...
                + c*kap4  ); %dorsal  MN
    d_nd_dt5 = (1/tau_n)*(-nd5 + ninf(vd5));
    d_md_dt5 = (1/tau_u)*( -gl*(md5 - El ) - syncon(vd5,1000,-30,20)*(md5 - Eplus) ); %dorsal muscle
    d_mv_dt5 = (1/tau_u)*( -gl*(mv5 - El ) - iin_link*syncon(vd5,1000,vnn,20)*(mv5 - Eminus) ...
                           - ein_link*syncon(vavb,1000,-30,20)*(mv5 - Estar) ); %ventral muscle
    d_kap_dt5 = (1/tau_b)*( -kap5 + sigma_def(mus(md5)) - sigma_def(mus(mv5)) ); %curvature
    
    
    %segment 6
    vd6 = x(26);
    nd6 = x(27);    
    md6 = x(28);
    mv6 = x(29);
    kap6 = x(30);
    
    d_vd_dt6 = (1/Cm)*( -gl*(vd6 - El) ...
                - gca*minf(vd6)*(vd6 - Eca) ...
                - gk*nd6*(vd6 - Ek) ...
                - syncon(vavb,100,-40,30)*(vd6 - Eplus )...
                + c*kap5  ); %dorsal  MN
    d_nd_dt6 = (1/tau_n)*(-nd6 + ninf(vd6));
    d_md_dt6 = (1/tau_u)*( -gl*(md6 - El ) - syncon(vd6,1000,-30,20)*(md6 - Eplus) ); %dorsal muscle
    d_mv_dt6 = (1/tau_u)*( -gl*(mv6 - El ) - iin_link*syncon(vd6,1000,vnn,20)*(mv6 - Eminus) ...
                           - ein_link*syncon(vavb,1000,-30,20)*(mv6 - Estar) ); %ventral muscle
d_kap_dt6 = (1/tau_b)*( -kap6 + sigma_def(mus(md6)) - sigma_def(mus(mv6)) ); %curvature
    
    dxdt=[d_vd_dt1; d_nd_dt1; d_md_dt1; d_mv_dt1;d_kap_dt1;...
          d_vd_dt2; d_nd_dt2; d_md_dt2; d_mv_dt2;d_kap_dt2;...
          d_vd_dt3; d_nd_dt3; d_md_dt3; d_mv_dt3;d_kap_dt3;...
          d_vd_dt4; d_nd_dt4; d_md_dt4; d_mv_dt4;d_kap_dt4;...
          d_vd_dt5; d_nd_dt5; d_md_dt5; d_mv_dt5;d_kap_dt5;...
          d_vd_dt6; d_nd_dt6; d_md_dt6; d_mv_dt6;d_kap_dt6        
        ];
end
​
​
%activation functions
function sig = sigma_def(x)
    cs = 0.01;
    a0 = 0;
    sig = (tanh(cs*(x-a0))+1)*1000;
end
​
function mi = minf(v)
    theta_m = 10.25;
    vm = -29;
    mi = 1./(1 + exp((vm-v)/(theta_m)));
end
​
function ni = ninf(v)
    theta_n = 20;
    vn = -55;
    ni = 1./(1 + exp((vn-v)/(theta_n)));
end
​
function ni = syncon(v,gbar,vn,theta_n)
    K=4.3944;
    ni = gbar./(1 + exp(K*((vn-v)/(theta_n))));
end
​
function mu = mus(v)
    theta_mus = 10;
    vmus = -45;
    mu = 1./(1 + exp((vmus-v)/(theta_mus)));
end
