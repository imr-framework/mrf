function display_epg(om_store,seq,annot,ax)
%DISPLAY_EPG(om_store, seq, annot, ax)
% Displays extended phase graph (EPG) given configiration state values and
% pulse sequence.
%
% om_store: cell array of config. state history, genereated by EPG_custom.m
% seq: struct containing sequence information (see EPG_custom.m)
% annot: 1 or 0, whether to display k-state populations
% ax: axis to display on (optional); if none, it defaults to a new figure

% Created by Sairam Geethanath
% Modified by Gehua Tong, Nov 14 2018

kmax = size(om_store{length(om_store)},2);
if kmax ~= 1
    kstates = -kmax+1:kmax-1; 
else 
    kstates = [-1,1];
end
grad_cnt =1; % indexing gradient
rf_cnt=1; % indexing rf pulse
timing = seq.time; 
uniqtimes = unique(timing);
grad = seq.grad;

if nargin < 4
    ax = gca;
end
%% For t==0 case - must be true for all sequences
% Plot horizontal axis (k=0)
plot(ax,[0 seq.time(end)],[0 0],'-k','LineWidth',1.5)
hold on
% Plot first RF pulse
t=0.*ones(1,length(kstates));
plot(ax,t,kstates,'-','Color',[0.5 0.5 0.5],'LineWidth',3);
flip = seq.rf(:,rf_cnt).'; % why .'???
text(ax,t(1),kmax-1,['(',num2str(flip(1)), '^{\circ}',',',num2str(flip(2)),'^{\circ}',')'] ,'FontSize',10);
rf_cnt = rf_cnt + 1;
% set axis for entire graph
%axis(ax,[0 timing(end) -kmax+1 kmax-1]);
%% For t > 0, plot on!
for seq_read = 2:length(seq.events) % for all events after the first pulse
% Get data
%           om_current = omega{seq_read};
            om_past  = om_store{seq_read -1};
            
            % Fp - states - 
%             Fpc = squeeze(om_current(1,:));
            Fpp = squeeze(om_past(1,:)); % all the +k states and k=0
            
            % Fm - states - 
%             Fmc = squeeze(om_current(2,:));
            Fmp = squeeze(om_past(2,:)); % all the -k states and k=0
               
            % Zk - states - 
%             Zc = squeeze(om_current(3,:));
            Zp = squeeze(om_past(3,:)); % all the Z states, k>=0
            
    switch (seq.events{seq_read})
        % --- Event type : RF pulse ---  
        case 'rf' %exchanges populations among three states; depict only 2
            % Draw vertical line spanning all k-states
            t = seq.time(seq_read).*ones(1,length(kstates));
            plot(ax,t,kstates,'Color',[0.5 0.5 0.5],'LineWidth',3); hold on; 
            % Get and label RF angles
            flip = seq.rf(:,rf_cnt).';
            text(ax,t(1),max(kmax-1,1),['(',num2str(flip(1)), '^{\circ}', ... 
                ',',num2str(flip(2)),'^{\circ}',')'] ,'FontSize',10);
            % Increase RF count
            rf_cnt = rf_cnt + 1;

       % --- Event type: Gradient ---
        case 'grad' % important: at a given time, grad always happens before rf
             % Increase gradient count
             grad_cnt = grad_cnt + 1;
             %(+) Fp state plot
             Fpp_kstates= find(abs(Fpp)> 5*eps) -1; % find nonzero states
             for k=1:length(Fpp_kstates) % for each +k state
                % vertical locations - [last_k, last_k + grad(in units of delk)]  
                  Fp_plot = [Fpp_kstates(k) Fpp_kstates(k)+grad(grad_cnt-1)];
                  t =  [uniqtimes(grad_cnt-1)  uniqtimes(grad_cnt)];
                  plot(ax,t,Fp_plot,'k-');hold on;
                  %-----------------------------------------------------------------
                  %Anotation of config. state value (a complex number for each line)
                  if(annot==1) 
                     intensity = round(Fpp(Fpp_kstates(k)+1)*100)/100;
                     text(ax,t(1),Fp_plot(1)+0.5,num2str(intensity),...
                              'Color',[0.01 0.58 0.53],'FontSize',9);
                  end
                  %-----------------------------------------------------------------
             end
             %(-) Fm state plot
             Fmp_kstates= -1*(find(abs(Fmp)> 5*eps) -1);
             for k=1:length(Fmp_kstates)
                 Fp_plot = [Fmp_kstates(k) Fmp_kstates(k)+grad(grad_cnt-1)];
                 t =  [uniqtimes(grad_cnt-1)  uniqtimes(grad_cnt)];
                 plot(ax,t,Fp_plot,'k-');hold on;
                 
                 % Echos
                 Fmp_echo = find(Fp_plot==0,1);
                 if ~isempty(Fmp_echo)
                     plot(ax,t(Fmp_echo), 0, '--ro','LineWidth',2,...
                                    'MarkerEdgeColor','k',...
                                    'MarkerFaceColor','g',...
                                    'MarkerSize',10);
                 end
                 if(annot==1)
                     intensity = round(Fmp(-Fmp_kstates(k)+1)*100)/100;
                     text(ax,t(1),Fp_plot(1)-0.5,num2str(intensity),...
                               'Color',[0.02 0.02 0.67],'FontSize',9);
                 end    
             end
             %Zp state plot
             Zp_kstates= (find(abs(Zp)> 5*eps) -1); 
             for k=1:length(Zp_kstates)
                 Fp_plot = [Zp_kstates(k) Zp_kstates(k)];
                 t =  [uniqtimes(grad_cnt-1)  uniqtimes(grad_cnt)];
                 plot(ax,t,Fp_plot,'--k');hold on;

                 if(annot==1) 
                     intensity = round(Zp(Zp_kstates(k)+1)*100)/100; 
                     text(ax,t(1),Fp_plot(1),num2str(intensity),...
                           'Color',[1 0.47 0.42],'FontSize',9);
                 end
             end
    end                
end

title(ax,seq.name,'fontsize',12);
xlabel(ax,'Time (ms)','fontsize',12);ylabel(ax,'k states','fontsize',12);
grid(ax,'ON');

%% Plot gradient 
baseline = -kmax-1;
M = length(uniqtimes);
for m = 2:M
    if grad(m-1)>0
       col = 'g'; % positive gradient in green
    else
       col = 'r'; % negative gradient in red
    end
    area(ax,[uniqtimes(m-1),uniqtimes(m)],...
        [baseline+grad(m-1),baseline+grad(m-1)],'FaceColor',col,'BASEVALUE',baseline)
    
end
% kmax here is the # columns in the last omega matrix
% so if kmax is 1, it means no nonzero gradients were used.  
if kmax == 1
    axis(ax,[0 timing(end) -1 1]);
else 
    axis(ax,[0 timing(end)  -kmax-1-max(abs(grad)) kmax-1]);
end

xticks(ax,uniqtimes)
