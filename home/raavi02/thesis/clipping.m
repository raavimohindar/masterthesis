function [y]=clipping(x,alpha1,alpha2)
ind1=find(abs(x)>alpha1);
ind2=find(abs(x)<alpha2);
y=x;
y(ind1)=sign(x(ind1));
y(ind2)=0;
