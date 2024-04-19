classdef Constants
   properties (Constant)
      w=0.056;
      fm=0.05337605126836238;%0.09245;
      ncycle=20;
      tao=4*41.34;
      n=256*2;
      m=256*2;
      x10=20;
      y10=0.2;
      z10=0;
      x20=-20;
      y20=0.2;
      z20=0;
      E0=2.9;
      I10=0.9;
      I20=2;
      Up=Constants.fm^2/4/Constants.w.^2;
      T0=2*pi/Constants.w;
      t0 = -1*Constants.T0; 
      tr=2*Constants.T0;
      V0=-2/sqrt(Constants.x10.^2+Constants.y10.^2+Constants.z10.^2)-...
          2/sqrt(Constants.x20.^2+Constants.y20.^2+Constants.z20.^2)+...
          1/sqrt((Constants.x10-Constants.x20).^2+(Constants.y10-Constants.y20).^2+(Constants.z10-Constants.z20).^2);
      v10=0.2;
      v20=0.2;
      y0 = [Constants.x10 Constants.y10 Constants.z10 ...
          0  0    Constants.v10  ...
          Constants.x20 Constants.y20 Constants.z20 ...
          0   0   0 ]';
   end

%     methods (Static)
%       end
%       end
end
