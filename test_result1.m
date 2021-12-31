clc;clear;close all;
%please download seistr package from https://github.com/chenyk1990/seistr
addpath(genpath('seistr/'));
addpath(genpath('subroutines/'));
% dir('../downloads_read/seg*.rsf')
% !mkdir -p bpfk
% !mkdir -p fk
% !mkdir -p mf
% !mkdir -p somf
% !mkdir -p bpfk
% !mkdir -p bpmffk
% !mkdir -p bpmf
% !mkdir -p bpsomffk
eq=zeros(2000,960);
% for ii=1:60
[n1,n2]=size(eq);
for ii=3
    if ~ismember(ii,[14,16,17,27,47,52])
        strcat('mat_raw/eq-',num2str(ii),'.mat')
        load(strcat('mat_raw/eq-',num2str(ii),'.mat'));
    end
    d1=d1;
    eq=d1;
    d1=das_bandpass(d1,0.0005,0,200,6,6,0,0);%
    d_bp=d1;
    figure(1);das_imagesc([eq,d1,eq-d1]);
    
    %     d11=das_mf(d1,20,1,2);%MF
    %     figure(2);das_imagesc([eq,d11,eq-d11]);
    
    %% SOMF
        [pp]=str_dip2d(d1,2,10,2,0.01, 1, 0.000001,[50,50,1],1);%figure;das_imagesc(pp);colormap(jet);
    ns=8;
    order=2;
    eps=0.01;
    type_mf=0;ifsmooth=0;
        [~,d1]=das_pwsmooth_lop_mf(pp,[],n1,n2,ns,order,eps,n1*n2,n1*n2,type_mf,ifsmooth,d1,[]);%SOMF
    	d1=reshape(d1,n1,n2);
%     load(strcat('mat_bpsomf/eq-',num2str(ii),'.mat'));
    %     save(strcat('mat_bpsomf/eq-',num2str(ii),'.mat'),'d1','-v7.3');
    figure(3);das_imagesc([eq,d1,eq-d1]);
    d_bpsomf=d1;
    
    d1=d1-das_fk_dip(d1,0.02);%
    d_bpsomffk=d1;
    %     load(strcat('mat_bpsomffk/eq-',num2str(ii),'.mat'));
    %     save(strcat('mat_bpsomffk/eq-',num2str(ii),'.mat'),'d1','-v7.3');
    
    figure(4);das_imagesc([eq,d1,eq-d1]);
    %     print(gcf,'-djpeg','-r300',strcat('bpsomffk/eq-',num2str(ii),'.jpg'));
end

%ii=3: FORGE_78-32_iDASv3-P11_UTC190423213209.sgy, 1484, 3.394402, 0.910045

%ii=10 is good

t=[0:n1]*0.0005;
x=1:n2;
ngap=50;

eq2=[eq,zeros(n1,ngap),zeros(size(eq))];
d_bp2=[d_bp,zeros(n1,ngap),eq-d_bp];
d_bpsomf2=[d_bpsomf,zeros(n1,ngap),eq-d_bpsomf];
d_bpsomffk2=[d_bpsomffk,zeros(n1,ngap),eq-d_bpsomffk];
x=1:ngap+n2*2;

figure('units','normalized','Position',[0.1 0.1 0.8, 0.9],'color','w');
subplot(2,2,1);das_imagesc(eq2,100,2,x,t);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Channel','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
text(n2/2,-0.05,'Raw data','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');


subplot(2,2,2);das_imagesc(d_bp2,100,2,x,t);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Channel','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
text(n2/2,-0.05,'BP','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.05,'Removed noise','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.679 0.700],...
    [0.772 0.866],'Color','r','TextColor','r','HorizontalAlignment','center',...
    'String',{'High-amplitude erratic noise'},...
    'LineWidth',2,...
    'FontSize',20,'fontweight','bold');
text(-200,-0.1,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

subplot(2,2,3);das_imagesc(d_bpsomf2,100,2,x,t);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Channel','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
text(n2/2,-0.05,'BP+SOMF','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.05,'Removed noise','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.169 0.190],...
    [0.332 0.426],'Color','r','TextColor','r','HorizontalAlignment','center',...
    'String',{'Horizontal noise'},...
    'LineWidth',2,...
    'FontSize',20,'fontweight','bold');
text(-200,-0.1,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

subplot(2,2,4);das_imagesc(d_bpsomffk2,100,2,x,t);
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Channel','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
text(n2/2,-0.05,'BP+SOMF+FK','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.05,'Removed noise','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

%zomming area
% annotation(gcf,'rectangle',...
%     [0.13 0.85 0.335 0.065]);
print(gcf,'-depsc','-r300','result1.eps');

inds1=1:400;
figure('units','normalized','Position',[0.1 0.1 0.8, 0.9],'color','w');
subplot(2,2,1);das_imagesc(eq2(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Channel','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
text(n2/2,-0.01,'Raw data','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(a)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');

subplot(2,2,2);das_imagesc(d_bp2(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Channel','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
text(n2/2,-0.01,'BP','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.01,'Removed noise','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(b)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.679 0.700],...
    [0.772 0.866],'Color','r','TextColor','r','HorizontalAlignment','center',...
    'String',{'High-amplitude erratic noise'},...
    'LineWidth',2,...
    'FontSize',20,'fontweight','bold');

subplot(2,2,3);das_imagesc(d_bpsomf2(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Channel','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
text(n2/2,-0.01,'BP+SOMF','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.01,'Removed noise','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(c)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.179 0.200],...
    [0.272 0.366],'Color','r','TextColor','r','HorizontalAlignment','center',...
    'String',{'Horizontal noise'},...
    'LineWidth',2,...
    'FontSize',20,'fontweight','bold');

subplot(2,2,4);das_imagesc(d_bpsomffk2(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',20,'fontweight','bold');
xlabel('Channel','Fontsize',20,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',20,'Fontweight','bold');
text(n2/2,-0.01,'BP+SOMF+FK','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.01,'Removed noise','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(d)','color','k','Fontsize',20,'fontweight','bold','HorizontalAlignment','center');
print(gcf,'-depsc','-r300','result2.eps');

%% combined figure
figure('units','normalized','Position',[0.0 0.0 0.5, 1],'color','w');
subplot(4,2,1);das_imagesc(eq2,100,2,x,t);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.05,'Raw data','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(a)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',[0.13 0.894 0.334 0.030],'linewidth',2,'color','g');

subplot(4,2,2);das_imagesc(d_bp2,100,2,x,t);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.05,'BP','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.05,'Removed noise','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.679 0.700],...
    [0.85 0.90],'Color','r','TextColor','r','HorizontalAlignment','center',...
    'String',{'High-amplitude erratic noise'},...
    'LineWidth',2,...
    'FontSize',10,'fontweight','bold');
text(-200,-0.1,'(b)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',[0.57 0.894 0.334 0.030],'linewidth',2,'color','g');

subplot(4,2,3);das_imagesc(d_bpsomf2,100,2,x,t);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.05,'BP+SOMF','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.05,'Removed noise','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.190 0.192],...
    [0.65 0.694],'Color','r','TextColor','r','HorizontalAlignment','center',...
    'String',{'Horizontal noise'},...
    'LineWidth',2,...
    'FontSize',10,'fontweight','bold');
text(-200,-0.1,'(c)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',[0.13 0.675 0.334 0.030],'linewidth',2,'color','g');

subplot(4,2,4);das_imagesc(d_bpsomffk2,100,2,x,t);
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.05,'BP+SOMF+FK','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.05,'Removed noise','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.1,'(d)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'rectangle',[0.57 0.675 0.334 0.030],'linewidth',2,'color','g');
%zomming area
% annotation(gcf,'rectangle',...
%     [0.13 0.85 0.335 0.065]);


inds1=1:400;
subplot(4,2,5);das_imagesc(eq2(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.01,'Raw data (zoomed)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(e)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');

subplot(4,2,6);das_imagesc(d_bp2(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.01,'BP (zoomed)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.01,'Removed noise','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(f)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.679 0.700],...
    [0.42 0.47],'Color','r','TextColor','r','HorizontalAlignment','center',...
    'String',{'High-amplitude erratic noise'},...
    'LineWidth',2,...
    'FontSize',10,'fontweight','bold');

subplot(4,2,7);das_imagesc(d_bpsomf2(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.01,'BP+SOMF (zoomed)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.01,'Removed noise','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(g)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
annotation(gcf,'textarrow',[0.189 0.210],...
    [0.202 0.246],'Color','r','TextColor','r','HorizontalAlignment','center',...
    'String',{'Horizontal noise'},...
    'LineWidth',2,...
    'FontSize',10,'fontweight','bold');

subplot(4,2,8);das_imagesc(d_bpsomffk2(inds1,:),100,2,x,t(inds1));
ylabel('Time (s)','Fontsize',10,'fontweight','bold');
xlabel('Channel','Fontsize',10,'fontweight','bold');
set(gca,'Linewidth',2,'Fontsize',10,'Fontweight','bold');
text(n2/2,-0.01,'BP+SOMF+FK (zoomed)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(n2+ngap+n2/2,-0.01,'Removed noise','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
text(-200,-0.02,'(h)','color','k','Fontsize',10,'fontweight','bold','HorizontalAlignment','center');
print(gcf,'-depsc','-r300','fig1.eps');



