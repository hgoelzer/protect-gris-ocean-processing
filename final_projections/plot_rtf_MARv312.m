load runoff_tf_MARv312.mat

t=[1950:2100];

% retreat time series plot
lspace = 0.06;
rspace = 0.025;
bspace = 0.055;
tspace = 0.03;
hspace = 0.03;
vspace = 0.06;
pw = (1-lspace-rspace-3*hspace)/4;
ph = (1-bspace-tspace-vspace)/2;
lwthick = 1;
lwthin = 0.5;
transp = 0.3;
fs = 8;

% colours
cols = [0.000,0.447,0.741;...
        0.850,0.325,0.098;...
        0.929,0.694,0.125;...
        0.494,0.184,0.556;...
        0.466,0.674,0.188;...
        0.301,0.745,0.933;...
        0.635,0.078,0.184;...
        0.000,0.000,0.000;...
        0.500,0.500,0.500;...
        1,0,0;...
        0,1,0;...
        0,0,1;...
        1,0,1;...
        1,1,0;...
        0,1,1];

figure();
a(1) = axes('position',[lspace+0*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(2) = axes('position',[lspace+1*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(3) = axes('position',[lspace+2*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(4) = axes('position',[lspace+3*(hspace+pw),bspace+1*(vspace+ph),pw,ph]);
a(5) = axes('position',[lspace+0*(hspace+pw),bspace+0*(vspace+ph),pw,ph]);
a(6) = axes('position',[lspace+1*(hspace+pw),bspace+0*(vspace+ph),pw,ph]);
a(7) = axes('position',[lspace+2*(hspace+pw),bspace+0*(vspace+ph),pw,ph]);

for l=1:7,
    
    axes(a(l)); hold on;
    plot(t,0*t,'--','linewidth',lwthin,'color',0.5*[1,1,1]);
    % ACCESS RCP8.5
    plot(t,retreat.ACCESS.RCP85.Qmed(l,:),'color',cols(1,:),'linewidth',lwthick);
    % CESM2 ssp585
    plot(t,retreat.CESM2Leo.ssp585.Qmed(l,:),'color',cols(2,:),'linewidth',lwthick);
    % CNRM-CM6-1 ssp585
    plot(t,retreat.CNRMCM6.ssp585.Qmed(l,:),'color',cols(3,:),'linewidth',lwthick);
    %% CNRM-ESM2-1 ssp585
    plot(t,retreat.CNRMESM2.ssp585.Qmed(l,:),'color',cols(4,:),'linewidth',lwthick);
    % MPIESM12HR ssp585
    plot(t,retreat.MPIESM12HR.ssp585.Qmed(l,:),'color',cols(5,:),'linewidth',lwthick);
    % MPIESM12HR ssp245
    plot(t,retreat.MPIESM12HR.ssp245.Qmed(l,:),'color',cols(6,:),'linewidth',lwthick);
    % MPIESM12HR ssp126
    plot(t,retreat.MPIESM12HR.ssp126.Qmed(l,:),'color',cols(7,:),'linewidth',lwthick);
    %% UKESM1-0-LL ssp585
    plot(t,retreat.UKESM1Robin.ssp585.Qmed(l,:),'color',cols(8,:),'linewidth',lwthick);

    xlim([2015 2100]); ylim([-1000 5]);
end

% Check runoff vs tf
