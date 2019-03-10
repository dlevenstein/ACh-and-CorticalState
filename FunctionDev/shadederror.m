function [h] = shadederror(x,y,z,col,edg_col,lw,mean_col)
%x=range of values over which to plot the shaded region (column vector)
%y=means (row vector)
%z=errors (row vector)
%col=color to fill shadederror
%edg_col=edge color of shadederror
%lw=linewidth for xy
%mean_col=line color of xy
%10/29/2013

%could specify some defaults



tempx=fliplr(x);
newx=cat(2,x,tempx);

temp1=y-z;
temp2=y+z;
newy=cat(1,temp1,flipud(temp2));

h=fill(newx,newy,col);
set(h,'EdgeColor',edg_col)
set(h,'FaceAlpha',0.5)
hold on

hold on
plot(x,y,'color', mean_col,'linewidth',lw)
xlim([min(x) max(x)])

end

