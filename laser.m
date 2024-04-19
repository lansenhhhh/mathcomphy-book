function [y,ff]=laser(t,w,fm)
%tao=Constants.tao;
%ff=exp(-2*log(2)*t.^2/tao^2);
ncycle=Constants.ncycle;
ff=(sin((pi*t)./(ncycle*2*pi/w))).^2;
ff=fm*ff;
y=ff.*sin(w*t);