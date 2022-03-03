function dtagaudit_plotRES(AXc,RES,XLIMS)
% Plot audit labels within defined (invisible) axis
%
% 

global foc_labels nonfoc_labels

axes(AXc), plot(0,0,'k*-') ;

if ~isempty(RES.cue)
   hold on
   kk = find(sum(RES.cue')'>XLIMS(1) & RES.cue(:,1)<=XLIMS(2)) ;
   if ~isempty(kk)
      for k=kk'
          cue = RES.cue(k,:);
          label = RES.stype{k};
          label_type = strtok(label);
          if any(contains(lower(label_type),nonfoc_labels))
              col = [238 045 046]/255 ; % Red
              pos = 0.1;
          elseif any(contains(lower(label_type),foc_labels))
              col = [000 139 071]/255 ; % Green      
              pos = 0.8;
          else
              col = [10 10 10]/255; % Black
              pos = 0.45;
          end
          plot([cue(1) sum(cue)],pos*ones(2,1),'k*-','color',col) ;
          text(max([XLIMS(1) cue(1)+0.1]),pos,label,'FontSize',8,'color',col,'VerticalAlignment','bottom') ;
      end
   end
   hold off
end

set(AXc,'XLim',XLIMS,'YLim',[0 1]) ;
bc = get(gcf,'Color') ;
set(AXc,'Box','off','XTick',[],'YTick',[],'XColor',bc,'YColor',bc,'Color',bc) ;
return