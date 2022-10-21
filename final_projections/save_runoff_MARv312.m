% get runoff time series per sector

close all

load ../glaciers/glaciers_MARv312.mat

%
for l=1:7,
    qsum = zeros(1,151);
    for g=ids(l).inds
        qsum = qsum + glaciers(g).ACCESS.RCP85.QJJA;
    end
    QJJAsect.ACCESS.RCP85.QJJA(l,:) = qsum; 
    qsum = zeros(1,151);
    for g=ids(l).inds
        qsum = qsum + glaciers(g).CESM2.ssp585.QJJA;
    end
    QJJAsect.CESM2.ssp585.QJJA(l,:) = qsum;
    qsum = zeros(1,151);
    for g=ids(l).inds
        qsum = qsum + glaciers(g).CNRMCM6.ssp585.QJJA;
    end
    QJJAsect.CNRMCM6.ssp585.QJJA(l,:) = qsum;
    qsum = zeros(1,151);
    for g=ids(l).inds
        qsum = qsum + glaciers(g).CNRMESM2.ssp585.QJJA;
    end
    QJJAsect.CNRMESM2.ssp585.QJJA(l,:) = qsum;
    qsum = zeros(1,151);
    for g=ids(l).inds
        qsum = qsum + glaciers(g).MPIESM12HR.ssp585.QJJA;
    end
    QJJAsect.MPIESM12HR.ssp585.QJJA(l,:) = qsum;
    qsum = zeros(1,151);
    for g=ids(l).inds
        qsum = qsum + glaciers(g).MPIESM12HR.ssp245.QJJA;
    end
    QJJAsect.MPIESM12HR.ssp245.QJJA(l,:) = qsum;
    qsum = zeros(1,151);
    for g=ids(l).inds
        qsum = qsum + glaciers(g).MPIESM12HR.ssp126.QJJA;
    end
    QJJAsect.MPIESM12HR.ssp126.QJJA(l,:) = qsum;
    qsum = zeros(1,151);
    for g=ids(l).inds
        % pad for missing 1950-1959 with 1960-1969
        qsum = qsum + [glaciers(g).UKESM1.ssp585.QJJA(1:10), glaciers(g).UKESM1.ssp585.QJJA];
    end
    QJJAsect.UKESM1.ssp585.QJJA(l,:) = qsum;
end

% Plotting

% Q timeseries plot
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
    plot(t,QJJAsect.ACCESS.RCP85.QJJA(l,:),'color',cols(1,:),'linewidth',lwthick);
    plot(t,QJJAsect.CESM2.ssp585.QJJA(l,:),'color',cols(2,:),'linewidth',lwthick);
    plot(t,QJJAsect.CNRMCM6.ssp585.QJJA(l,:),'color',cols(3,:),'linewidth',lwthick);
    plot(t,QJJAsect.CNRMESM2.ssp585.QJJA(l,:),'color',cols(4,:),'linewidth',lwthick);
    plot(t,QJJAsect.MPIESM12HR.ssp585.QJJA(l,:),'color',cols(5,:),'linewidth',lwthick);
    plot(t,QJJAsect.MPIESM12HR.ssp245.QJJA(l,:),'color',cols(6,:),'linewidth',lwthick);
    plot(t,QJJAsect.MPIESM12HR.ssp126.QJJA(l,:),'color',cols(7,:),'linewidth',lwthick);
    plot(t,QJJAsect.UKESM1.ssp585.QJJA(l,:),'color',cols(8,:),'linewidth',lwthick);

    xlim([1950 2100]); ylim([0 50000]);
    set(gca,'fontsize',fs,'box','on');
    if l==1 | l==5, ylabel('proj. QJJA (m3 s-1)','fontsize',fs); end
    if l<5, set(gca,'xticklabel',[]); end
    if l~=1 & l~=5, set(gca,'yticklabel',[]); end
    text(2020,48000,sectorname{l},'fontsize',fs,'color','k');
    
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

saveplot(17,8,300,'QJJA_projections_MARv312.png');
%close all;

