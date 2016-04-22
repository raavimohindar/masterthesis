close all;
clear all;

x = -4:0.01:4;

alpha = [0.5 1 2 4];

for i = 1:length(alpha)

	plot(x,tanh(x*alpha(i)))
	hold on;

end
x = -1:0.01:1;
y = -1:0.01:1;

x1 = -.2:0.0001:0.2;
y1 = -1:4.9988e-04:1;
hold on;
	
	d = plot(2,tanh(2*0.5),'-*');
	f = plot(0.75,tanh(0.75*1),'-o');
	e = plot(0.5,tanh(0.5*2),'-+');
	c = plot(0.25,tanh(0.25*4),'-x');
	g = plot(x,y,'--');
	h = plot(x1,y1,'linewidth',1,'linestyle','-');
legend([d,f,e,c,g,h],'aaaaaa','bbbbbb','cccccc','dddddd','eeeeee','ffffff','location','northeastoutside');

axis([-4 4 -1 1]);
grid on
set(gca,'ytick',[-1 -0.5 0 0.5 1]);
set(gca,'xtick',[-4 -3 -2 -1 0 1 2 3 4]);

xlabel('x');
ylabel('y');
title('t');

set(gcf,'paperposition',[0 0 4.5 3.0]);
set(gca,'fontname','times','fontsize',10);

print -deps non_linear_dev_psfrag.eps

