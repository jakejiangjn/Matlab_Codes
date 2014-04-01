function btitle(haxes)
% Make the title under the figure(at the bottom))
if nargin == 0
	haxes = gca;  % Default status: get current axis handle
end
%%
% Re-arrange axis's position
%set(gca,'position',[0.13 0.11 0.1 0.815]) % Default Value
set(haxes,'position',[0.13 0.15 0.775 0.815]); 
%% 
%Re-arrange title's position
htitle=get(gca,'title');
AXIS=axis;
pos=[AXIS(1) AXIS(3)]+[AXIS(2)-AXIS(1) AXIS(4)-AXIS(3)].*[0.5 -0.1];
set(htitle,'position',pos)%Position = [0.498975 1.01597 1.00005]
set(htitle,'VerticalAlignment','top');