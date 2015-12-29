% clear all; 
close all;

load lena;
T = -50:50;

i = 0;
for t=T
    i = i+1;
    
    I = zeros(m,n);
    
    if t>=0
    I(:,t+1:end) = I_ref(:,1:end-t);
    else
        t = -t;
        I(:,1:end-t) = I_ref(:,t+1:end);
    end
    
    f_ssd(i) = reg_similarity(I_source2,I,'ssd');
    f_rc(i) = reg_similarity(I_source2,I,'rc');
    f_sad(i) = reg_similarity(I_source2,I,'sad');
    f_cc(i) = reg_similarity(I_source2,I,'cc');
    f_cd2(i) = reg_similarity(I_source2,I,'cd2');
    f_ms(i) = reg_similarity(I_source2,I,'ms');
    f_mi(i) = reg_similarity(I_source2,I,'mi');
    f_dtv(i) = reg_similarity(I_source2,I,'dtv');
  
end

ls = 2.5;
figure;
h = subplot(2,4,1);plot(T,f_ssd,'r','linewidth',ls); title('SSD', 'fontsize', 14,'Color', 'b');
set(h,'Tag','1');
h =subplot(2,4,2);plot(T,f_rc,'r','linewidth',ls); title('RC', 'fontsize', 14,'Color', 'b');
set(h,'Tag','2');
h =subplot(2,4,3);plot(T,f_sad,'r','linewidth',ls); title('SAD', 'fontsize', 14,'Color', 'b');
set(h,'Tag','3');
h =subplot(2,4,4);plot(T,f_cc,'r','linewidth',ls); title('CC', 'fontsize', 14,'Color', 'b');
set(h,'Tag','4');
h =subplot(2,4,5);plot(T,f_cd2,'r','linewidth',ls); title('CD2', 'fontsize', 14,'Color', 'b');
set(h,'Tag','5');
h =subplot(2,4,6);plot(T,f_ms,'r','linewidth',ls); title('MS', 'fontsize', 14,'Color', 'b');
set(h,'Tag','6');
h =subplot(2,4,7);plot(T,f_mi,'r','linewidth',ls); title('MI', 'fontsize', 14,'Color', 'b');
set(h,'Tag','7');
h =subplot(2,4,8);plot(T,f_dtv,'r','linewidth',ls); title('Proposed', 'fontsize', 14,'Color', 'b');
set(h,'Tag','8');

i = 0;
for t=T
    i = i+1;
    
    I = zeros(m,n);
    
    if t>=0
    I(:,t+1:end) = I_ref(:,1:end-t);
    else
        t = -t;
        I(:,1:end-t) = I_ref(:,t+1:end);
    end
    
    f_ssd(i) = reg_similarity(I_source,I,'ssd');
    f_rc(i) = reg_similarity(I_source,I,'rc');
    f_sad(i) = reg_similarity(I_source,I,'sad');
    f_cc(i) = reg_similarity(I_source,I,'cc');
    f_cd2(i) = reg_similarity(I_source,I,'cd2');
    f_ms(i) = reg_similarity(I_source,I,'ms');
    f_mi(i) = reg_similarity(I_source,I,'mi');
    f_dtv(i) = reg_similarity(I_source,I,'dtv');
    
end

ls = 1.0; ts = 12;
% figure;
h = findobj('Tag','1'); % get handle to object tagged as 'left'
set(h,'NextPlot','add'); % set 'NextPlot' property to 'add'
plot(h,T,f_ssd,'linewidth',ls); 
xlabel(h,'Translation (pixels)');
ylabel(h,'Function value');
textobj = findobj(h,'type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(h,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(h,'YLabel');
set(h_xlabel,'FontSize',ts); 

h = findobj('Tag','2'); % get handle to object tagged as 'left'
set(h,'NextPlot','add'); % set 'NextPlot' property to 'add'
plot(h,T,f_rc,'linewidth',ls); 
xlabel(h,'Translation (pixels)');
ylabel(h,'Function value');
textobj = findobj(h,'type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(h,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(h,'YLabel');
set(h_xlabel,'FontSize',ts); 

h = findobj('Tag','3'); % get handle to object tagged as 'left'
set(h,'NextPlot','add'); % set 'NextPlot' property to 'add'
plot(h,T,f_sad,'linewidth',ls); 
xlabel(h,'Translation (pixels)');
ylabel(h,'Function value');
textobj = findobj(h,'type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(h,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(h,'YLabel');
set(h_xlabel,'FontSize',ts); 

h = findobj('Tag','4'); % get handle to object tagged as 'left'
set(h,'NextPlot','add'); % set 'NextPlot' property to 'add'
plot(h,T,f_cc,'linewidth',ls); title('CC');
xlabel(h,'Translation (pixels)');
ylabel(h,'Function value');
textobj = findobj(h,'type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(h,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(h,'YLabel');
set(h_xlabel,'FontSize',ts); 

h = findobj('Tag','5'); % get handle to object tagged as 'left'
set(h,'NextPlot','add'); % set 'NextPlot' property to 'add'
plot(h,T,f_cd2,'linewidth',ls); title('CD2');
xlabel(h,'Translation (pixels)');
ylabel(h,'Function value');
textobj = findobj(h,'type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(h,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(h,'YLabel');
set(h_xlabel,'FontSize',ts); 

h = findobj('Tag','6'); % get handle to object tagged as 'left'
set(h,'NextPlot','add'); % set 'NextPlot' property to 'add'
plot(h,T,f_ms,'linewidth',ls); title('MS');
xlabel(h,'Translation (pixels)');
ylabel(h,'Function value');
textobj = findobj(h,'type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(h,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(h,'YLabel');
set(h_xlabel,'FontSize',ts); 

h = findobj('Tag','7'); % get handle to object tagged as 'left'
set(h,'NextPlot','add'); % set 'NextPlot' property to 'add'
plot(h,T,f_mi,'linewidth',ls); title('MI');
xlabel(h,'Translation (pixels)');
ylabel(h,'Function value');
textobj = findobj(h,'type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(h,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(h,'YLabel');
set(h_xlabel,'FontSize',ts); 

h = findobj('Tag','8'); % get handle to object tagged as 'left'
set(h,'NextPlot','add'); % set 'NextPlot' property to 'add'
plot(h,T,f_dtv,'linewidth',ls); title('Proposed', 'fontsize', 14);

% legend('without intensity distortions','with intensity distortions');
xlabel(h,'Translation (pixels)');
ylabel(h,'Function value');

textobj = findobj(h,'type', 'text');
set(textobj, 'fontsize', ts);
h_xlabel = get(h,'XLabel');
set(h_xlabel,'FontSize',ts); 
h_xlabel = get(h,'YLabel');
set(h_xlabel,'FontSize',ts); 

