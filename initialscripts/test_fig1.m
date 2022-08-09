% Script to plot Figure 1
% Aug, 2020
% Yangkang Chen

clc;clear;close all;
addpath(genpath('../subroutines/'));

nt=1000;dt=0.004;t=[0:nt-1]*dt-1.8;
d=zeros(nt,1);
d(500,1)=1;
d1=das_bandpass(d,dt,0,50,6,6,0,0);

%% trapezoid BP
f1=0;f2=0;f3=45;f4=55;
nt=length(d);
k = nextpow2(nt);
nf = 4*(2^k);
if ~(f1==0 && f2==0)
    i1 = floor(nf*f1*dt)+1;
    i2 = floor(nf*f2*dt)+1;
    up =  (1:1:(i2-i1))/(i2-i1);
end
i3 = floor(nf*f3*dt)+1;
i4 = floor(nf*f4*dt)+1;
down = (i4-i3:-1:1)/(i4-i3);
if f1==0 && f2==0
    aux = [ones(1,i3), down, zeros(1,nf/2+1-i4) ];
else
    aux = [zeros(1,i1), up, ones(1,i3-i2), down, zeros(1,nf/2+1-i4) ];
end
aux2 = fliplr(aux(1,2:nf/2));
c = 0; 
F = ([aux,aux2]');
Phase = (pi/180.)*[0.,-c*ones(1,nf/2-1),0.,c*ones(1,nf/2-1)];
Transfer = F.*exp(-i*Phase');
D = fft(d,nf,1);
Do = Transfer.*D;
d2 = ifft(Do,nf,1);
d2 = real(d2(1:nt,:));


%% Frequency spectrum
[fd, f] = das_fft1(d, nt, dt);
[fd1, f] = das_fft1(d1,nt, dt);
[fd2, f] = das_fft1(d2,nt, dt);

%% visualize
figure('units','normalized','Position',[0.2 0.4 0.4, 1],'color','w');
subplot(6,2,1);plot(t,d,'linewidth',2);xlim([0,0.4]);ylim([-0.2,1]);ylabel('Amplitude','Fontsize',15,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-0.05,1.6,'(a)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
subplot(6,2,3);plot(t,d1,'linewidth',2);xlim([0,0.4]);ylim([-0.2,1]);ylabel('Amplitude','Fontsize',15,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-0.05,1.6,'(c)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
subplot(6,2,5);plot(t,d2,'linewidth',2);xlim([0,0.4]);ylim([-0.2,1]);ylabel('Amplitude','Fontsize',15,'fontweight','bold');set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-0.05,1.6,'(e)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
xlabel('Time (s)','Fontsize',15,'fontweight','bold');
annotation(gcf,'textarrow',[0.333 0.322],...
    [0.611 0.582],'String',{'fluctuation'},'linewidth',2,'color','r','Fontsize',10,'fontweight','bold');

subplot(6,2,2);plot(f,abs(fd),'linewidth',2);ylim([0,1]);set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-20,1.6,'(b)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
subplot(6,2,4);plot(f,abs(fd1),'linewidth',2);ylim([0,1]);set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-20,1.6,'(d)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
subplot(6,2,6);plot(f,abs(fd2),'linewidth',2);ylim([0,1]);set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');text(-20,1.6,'(f)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
xlabel('Frequency (Hz)','Fontsize',15,'fontweight','bold');
%% comparison
nphis=[2,4,6,8,10,12,14,16,18,20];
D=zeros(nt,length(nphis));
les=[];
for i=1:length(nphis)
    D(:,i)=das_bandpass(d,dt,0,50,6,nphis(i),0,0);
end
FD= das_fft1(D, nt, dt);

subplot(2,1,2);
plot(f,abs(FD(:,1)),'linewidth',2);hold on;
les{1}=strcat('N=',num2str(nphis(1)));
for i=2:length(nphis)
    plot(f,abs(FD(:,i)),'linewidth',2);
    les{end,end+1}=strcat('N=',num2str(nphis(i)));
end
legend(les);
ylim([0,1.01]);
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');
xlabel('Frequency (Hz)','Fontsize',15,'Fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',15,'Fontweight','bold');%ytickformat('%.1e');
ylabel('Amplitude','Fontsize',15,'fontweight','bold');
text(-10,1.1,'(g)','color','k','Fontsize',15,'fontweight','bold','HorizontalAlignment','center');
print(gcf,'-depsc','-r300','fig1.eps');


