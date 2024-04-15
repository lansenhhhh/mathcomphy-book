 %Modeling %%Transfer Function 
 mt = 1;
 mp = 0.2; b = 0.1; I = 0.03;
 g = 9.8;
 l = 0.5;
 q = (mt+mp)*(I+mp*l^2)-(mp*l)^2;
 s = tf('s');
 P_cart = (((I+mp*l^2)/q)*s^2 - (mp*g*l/q))/(s^4 + (b*(I + mp*l^2))*s^3/q - ((mt + mp)*mp*g*l)*s^2/q - b*mp*g*l*s/q);
 P_pend = (mp*l*s/q)/(s^3 + (b*(I + mp*l^2))*s^2/q - ((mt + mp)*mp*g*l)*s/q - b*mp*g*l/q);
  sys_tf = [P_cart ; P_pend] 
  inputs = {'u'}; outputs = {'x'; 'phi'};
  set(sys_tf,'InputName',inputs)
  set(sys_tf,'OutputName',outputs)