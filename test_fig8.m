
clc;clear;close all;
% addpath(genpath('~/chenyk/matlibcyk'));

%% load data
eq=zeros(2000,960);
[n1,n2]=size(eq);
ii=3;%reasonable
ii=1;
if ~ismember(ii,[14,16,17,27,47,52])
    rsf_read(eq,strcat('../downloads_read/seg-',num2str(ii),'.rsf'));
end
d1=eq;
save(strcat('mat_raw/eq-',num2str(ii),'.mat'),'d1','-v7.3');

load(strcat('mat_raw/eq-',num2str(ii),'.mat'));

load(strcat('mat_bpsomffk/eq-',num2str(ii),'.mat'));
figure;yc_imagesc([eq,d1,eq-d1]);

inds=20:20:n2; %works fine
inds=100:10:n2;
inds=inds(1:30);
traces=d1(:,inds);
traces0=eq(:,inds);

dn=traces;
d0=traces0;
nsta=30;nlta=300;
O=ones(size(inds));O0=ones(size(inds));
[ O,R ] = das_picker_stalta(dn,nsta, nlta);
[ O0,R0 ] = das_picker_stalta(d0,nsta, nlta);
% figure;subplot(1,2,1);yc_wigb_p(dn,O);
% subplot(1,2,2);yc_wigb_p(d0,O0);

times0=[O0-1]*0.0005;
times=[O-1]*0.0005;
%name
for ii=1:30
    stname{ii}=strcat('Channel:',num2str(inds(ii)));
end

%% begin plotting
figure('units','normalized','Position',[0.0 0.0 0.45, 1],'color','w');
nr=15;%number of stations in the first column
x0=0.1;y0=0.1;dy=0;dx=0;
%length: 0.5x0.5, 0.5x0.25, 0.25x0.5
%% axis XY
dh=(1-0.2)/nr;dw=0.37;
dh1=dh/2;%axis height
t=[0:1999]*0.0005;
for ir=nr:-1:1
%     fid=fopen(strcat('/Users/chenyk/chenyk/texnet/texnet2020galz/process/stations_binary/',stname{ir},'.bin'),'r');
%     wav=fread(fid,[12000,1],'double');
    wav0=traces0(:,ir);
    wav=traces(:,ir);
%     wav=yc_scale(wav(:),1);
    a1=axes('Parent',gcf,'Position',[x0,y0+dy+dh*(nr-ir),dw,dh1]);
    plot(t,wav0,'k','linewidth',2); hold on; axis off;
    plot(t,wav,'r','linewidth',2);
    plot([times0(ir),times0(ir)],[min(wav),max(wav)],'g','linewidth',2);
    plot([times(ir),times(ir)],[min(wav),max(wav)],'b','linewidth',2);
    ylim([min(wav),max(wav)]);
%     if ir~=7
%     ylim([-1,1]);
%     else
%     ylim([-1,1]);
%     end
    
%     fid=fopen(strcat('/Users/chenyk/chenyk/texnet/texnet2020galz/process/stations_binary/',stname{ir+15},'.bin'),'r');
%     wav=fread(fid,[12000,1],'double');
    wav0=traces0(:,ir+15);
    wav=traces(:,ir+15);
%     wav=yc_scale(wav(:),1);
    a1=axes('Parent',gcf,'Position',[x0+0.5,y0+dy+dh*(nr-ir),dw,dh1]);
    plot(t,wav0,'k','linewidth',2);hold on; axis off; 
    plot(t,wav,'r','linewidth',2);
    plot([times0(ir+15),times0(ir+15)],[min(wav),max(wav)],'g','linewidth',2);
    plot([times(ir+15),times(ir+15)],[min(wav),max(wav)],'b','linewidth',2);
    ylim([min(wav),max(wav)]);
    
end
legend('Raw waveform','Denoised waveform','Picked arrival from raw data','Picked arrival from denoised data','Position',[x0+0.15,y0-0.1,0.6,0.1],'NumColumns',4);
legend('boxoff');
% 
%% add station name
for ir=nr:-1:1
a1=axes('Parent',gcf,'Position',[0.02,y0+dh*(nr-ir)+dh/4,dw,dh1]);
text(-0.035,0,stname{ir},'color','k','Fontsize',10,'fontweight','bold');axis off;

a1=axes('Parent',gcf,'Position',[0.02+0.5,y0+dh*(nr-ir)+dh/4,dw,dh1]);
text(-0.035,0,stname{ir+15},'color','k','Fontsize',10,'fontweight','bold');axis off;
end
%
%
%% add source info
dw2=(1-x0)/5.0;
a1=axes('Parent',gcf,'Position',[0,0.93,1,dh1]);
text(0.5,0,'Earthquake detection of FORGE\_78-32\_iDASv3-P11\_UTC190423150554.sgy','color','k','Fontsize',14,'fontweight','bold','HorizontalAlignment','center');axis off;
print(gcf,'-depsc','-r300','fig8.eps');



% 
% s=d1(81:160,750:850);%s=d1(:,750:800);
% nn=d1(1:80,750:850);%nn=d1(:,750:800);
% nn2=eq(1:80,750:850);%nn2=eq(:,750:800);
% figure;yc_imagesc([nn2,s,nn2-s]);
% 
% yc_snr(s,nn) %denoised
% yc_snr(s,nn2)%raw
%-3.9598->4.0248



