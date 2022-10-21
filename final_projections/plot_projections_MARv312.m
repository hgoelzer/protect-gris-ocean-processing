%% plotting

load projected_retreat_MARv312.mat 

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
    plot(t,retreat.ACCESS.RCP85.high(l,:),'color',cols(1,:),'linewidth',lwthick);
    plot(t,retreat.ACCESS.RCP85.med(l,:),'color',cols(1,:),'linewidth',lwthick);
    plot(t,retreat.ACCESS.RCP85.low(l,:),'color',cols(1,:),'linewidth',lwthick); 
    % CESM2 ssp585
    plot(t,retreat.CESM2.ssp585.high(l,:),'color',cols(2,:),'linewidth',lwthick);
    plot(t,retreat.CESM2.ssp585.med(l,:),'color',cols(2,:),'linewidth',lwthick);
    plot(t,retreat.CESM2.ssp585.low(l,:),'color',cols(2,:),'linewidth',lwthick); 
    % CNRM-CM6-1 ssp585
    plot(t,retreat.CNRMCM6.ssp585.high(l,:),'color',cols(3,:),'linewidth',lwthick);
    plot(t,retreat.CNRMCM6.ssp585.med(l,:),'color',cols(3,:),'linewidth',lwthick);
    plot(t,retreat.CNRMCM6.ssp585.low(l,:),'color',cols(3,:),'linewidth',lwthick); 
    %% CNRM-ESM2-1 ssp585
    plot(t,retreat.CNRMESM2.ssp585.high(l,:),'color',cols(4,:),'linewidth',lwthick);
    plot(t,retreat.CNRMESM2.ssp585.med(l,:),'color',cols(4,:),'linewidth',lwthick);
    plot(t,retreat.CNRMESM2.ssp585.low(l,:),'color',cols(4,:),'linewidth',lwthick); 
    % MPIESM12HR ssp585
    plot(t,retreat.MPIESM12HR.ssp585.high(l,:),'color',cols(5,:),'linewidth',lwthick);
    plot(t,retreat.MPIESM12HR.ssp585.med(l,:),'color',cols(5,:),'linewidth',lwthick);
    plot(t,retreat.MPIESM12HR.ssp585.low(l,:),'color',cols(5,:),'linewidth',lwthick); 
    % MPIESM12HR ssp245
    plot(t,retreat.MPIESM12HR.ssp245.high(l,:),'color',cols(6,:),'linewidth',lwthick);
    plot(t,retreat.MPIESM12HR.ssp245.med(l,:),'color',cols(6,:),'linewidth',lwthick);
    plot(t,retreat.MPIESM12HR.ssp245.low(l,:),'color',cols(6,:),'linewidth',lwthick); 
    % MPIESM12HR ssp126
    plot(t,retreat.MPIESM12HR.ssp126.high(l,:),'color',cols(7,:),'linewidth',lwthick);
    plot(t,retreat.MPIESM12HR.ssp126.med(l,:),'color',cols(7,:),'linewidth',lwthick);
    plot(t,retreat.MPIESM12HR.ssp126.low(l,:),'color',cols(7,:),'linewidth',lwthick); 
    %% UKESM1-0-LL ssp585
    plot(t,retreat.UKESM1.ssp585.high(l,:),'color',cols(8,:),'linewidth',lwthick);
    plot(t,retreat.UKESM1.ssp585.med(l,:),'color',cols(8,:),'linewidth',lwthick);
    plot(t,retreat.UKESM1.ssp585.low(l,:),'color',cols(8,:),'linewidth',lwthick); 

    xlim([2015 2100]); ylim([-45 5]);
    set(gca,'fontsize',fs,'box','on');
    if l==1 | l==5, ylabel('proj. retreat $\Delta L$ (km)','fontsize',fs); end
    if l<5, set(gca,'xticklabel',[]); end
    if l~=1 & l~=5, set(gca,'yticklabel',[]); end
    text(2020,-23,sectorname{l},'fontsize',fs,'color','k');
    
end

tx = 0.035;
ty = 0.51;
dx = -0.02;
annotation('textbox','position',[tx+0*(pw+hspace)+dx,ty+1*(ph+vspace),0,0],'string','a','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+1*(pw+hspace),ty+1*(ph+vspace),0,0],'string','b','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+2*(pw+hspace),ty+1*(ph+vspace),0,0],'string','c','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+3*(pw+hspace),ty+1*(ph+vspace),0,0],'string','d','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+0*(pw+hspace)+dx,ty+0*(ph+vspace),0,0],'string','e','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+1*(pw+hspace),ty+0*(ph+vspace),0,0],'string','f','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+2*(pw+hspace),ty+0*(ph+vspace),0,0],'string','g','fontsize',fs,'fontweight','bold','edgecolor','none');

% manual legend
tx = 0.78;
ty = 0.15;
dy = 0.04;
fs0 = 4;
annotation('textbox','position',[tx,ty+1*dy,0,0],'string','ACCESS-RCP8.5','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(1,:));
annotation('textbox','position',[tx,ty+2*dy,0,0],'string','CESM2-ssp585','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(2,:));
annotation('textbox','position',[tx,ty+3*dy,0,0],'string','CNRM-CM6-1-ssp585','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(3,:));
annotation('textbox','position',[tx,ty+4*dy,0,0],'string','CNRM-ESM2-1-ssp585','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(4,:));
annotation('textbox','position',[tx,ty+5*dy,0,0],'string','MPIESM12HR-ssp585','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(5,:));
annotation('textbox','position',[tx,ty+6*dy,0,0],'string','MPIESM12HR-ssp245','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(6,:));
annotation('textbox','position',[tx,ty+7*dy,0,0],'string','MPIESM12HR-ssp126','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(7,:));
annotation('textbox','position',[tx,ty+8*dy,0,0],'string','UKESM1-0-LL-ssp585','fontsize',fs0,'interpreter','latex','edgecolor','none','color',cols(8,:));

saveplot(17,8,300,'retreat_projections_MARv312.png');
close all;

% TF timeseries plot
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
    plot(t,retreat.ACCESS.RCP85.TF(l,:),'color',cols(1,:),'linewidth',lwthick);
    plot(t,retreat.CESM2.ssp585.TF(l,:),'color',cols(2,:),'linewidth',lwthick);
    plot(t,retreat.CNRMCM6.ssp585.TF(l,:),'color',cols(3,:),'linewidth',lwthick);
    plot(t,retreat.CNRMESM2.ssp585.TF(l,:),'color',cols(4,:),'linewidth',lwthick);
    plot(t,retreat.MPIESM12HR.ssp585.TF(l,:),'color',cols(5,:),'linewidth',lwthick);
    plot(t,retreat.MPIESM12HR.ssp245.TF(l,:),'color',cols(6,:),'linewidth',lwthick);
    plot(t,retreat.MPIESM12HR.ssp126.TF(l,:),'color',cols(7,:),'linewidth',lwthick);
    plot(t,retreat.UKESM1.ssp585.TF(l,:),'color',cols(8,:),'linewidth',lwthick);

    xlim([2015 2100]); ylim([1 11]);
    set(gca,'fontsize',fs,'box','on');
    if l==1 | l==5, ylabel('proj. TF ($^{\circ}$C)','fontsize',fs); end
    if l<5, set(gca,'xticklabel',[]); end
    if l~=1 & l~=5, set(gca,'yticklabel',[]); end
    text(2020,9.8,sectorname{l},'fontsize',fs,'color','k');
    
end

tx = 0.035;
ty = 0.51;
dx = -0.02;
annotation('textbox','position',[tx+0*(pw+hspace)+dx,ty+1*(ph+vspace),0,0],'string','a','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+1*(pw+hspace),ty+1*(ph+vspace),0,0],'string','b','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+2*(pw+hspace),ty+1*(ph+vspace),0,0],'string','c','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+3*(pw+hspace),ty+1*(ph+vspace),0,0],'string','d','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+0*(pw+hspace)+dx,ty+0*(ph+vspace),0,0],'string','e','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+1*(pw+hspace),ty+0*(ph+vspace),0,0],'string','f','fontsize',fs,'fontweight','bold','edgecolor','none');
annotation('textbox','position',[tx+2*(pw+hspace),ty+0*(ph+vspace),0,0],'string','g','fontsize',fs,'fontweight','bold','edgecolor','none');

% manual legend
tx = 0.78;
ty = 0.5;
dy = 0.03;
fs = 5;
annotation('textbox','position',[tx,ty-1*dy,0,0],'string','ACCESS-RCP8.5','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(1,:));
annotation('textbox','position',[tx,ty-2*dy,0,0],'string','CESM2-ssp585','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(2,:));
annotation('textbox','position',[tx,ty-3*dy,0,0],'string','CNRM-CM6-1-ssp585','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(3,:));
annotation('textbox','position',[tx,ty-4*dy,0,0],'string','CNRM-ESM2-1-ssp585','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(4,:));
annotation('textbox','position',[tx,ty-5*dy,0,0],'string','MPIESM12HR-ssp585','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(5,:));
annotation('textbox','position',[tx,ty-6*dy,0,0],'string','MPIESM12HR-ssp245','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(6,:));
annotation('textbox','position',[tx,ty-7*dy,0,0],'string','MPIESM12HR-ssp126','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(7,:));
annotation('textbox','position',[tx,ty-8*dy,0,0],'string','UKESM1-0-LL-ssp585','fontsize',fs,'interpreter','latex','edgecolor','none','color',cols(8,:));

saveplot(17,8,300,'TF_projections_MARv312.png');
close all;

