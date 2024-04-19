clear;
T=Constants.T0;
omg=Constants.w;
fm=Constants.fm;
T0=Constants.T0;
n=Constants.n;
m=Constants.m;
x10=Constants.x10;
y10=Constants.y10;
z10=Constants.z10;
x20=Constants.x20;
y20=Constants.y20;
z20=Constants.z20;
V0=Constants.V0;
ncycle=Constants.ncycle;
v10=Constants.v10;
Up=Constants.Up;
Ip=0.9;
y0 =Constants.y0;
E2_data=[];
E1_data=[];
E_data=[];
E2_data1=[];
E1_data1=[];
E_data1=[];
E2_data2=[];
E1_data2=[];
E_data2=[];
t0_start=8.7;
t0_end=10.7;
t0_allpoints=(t0_end-t0_start)*n;
t0_step=1./n;
tr=4;
t_allpoints=tr*m;
t_step=1/m;
options = odeset('AbsTol',1e-6,'RelTol', 1e-6);

for i=1:t0_allpoints
    t0=(i*t0_step+t0_start)*T0;
    yend=y0;
    return_num=0;
    for j=1:t_allpoints
        ts=t0+(j-1)*t_step*T0;
        te=t0+(j)*t_step*T0;
        tspan = linspace(ts,te,2);
        ytemp1=yend(3);
        ytemp2=yend(9);
        [t,y] = ode45(@(t,y)he3d(t,y,omg,fm), tspan, yend,options);
        yend=y(end,:);
        E1 = 1/2 * (yend(4)^2 + yend(5)^2 + yend(6)^2) - 2/sqrt(yend(1)^2 + yend(2)^2 + yend(3)^2);
        E2 = 1/2 * (yend(10)^2 + yend(11)^2 + yend(12)^2) - 2/sqrt(yend(7)^2 + yend(8)^2 + yend(9)^2);
        E=E1+E2+ 1/sqrt((yend(1)-yend(7))^2 + (yend(2)-yend(8))^2 + (yend(3)-yend(9))^2);
        IsReturn1=ytemp1*(yend(3))<0;
        IsReturn2=ytemp2*(yend(9))<0;
        if( IsReturn1||IsReturn2)
            %             E_data=[E_data,E1(j)];
            %             t_data=[t_data,t(j)-t0];
            if(return_num==0)
                if(IsReturn1)
                    E1_data=[E1_data;[t(end),E1]];
                end
                if(IsReturn2)
                    E2_data=[E2_data;[t(end),E2]];
                end
                E_data=[E_data;[t(end),E]];
            elseif(return_num==1)
                if(IsReturn1)
                    E1_data1=[E1_data1;[t(end),E1]];
                end
                if(IsReturn2)
                    E2_data1=[E2_data1;[t(end),E2]];
                end
                E_data1=[E_data1;[t(end),E]];
            elseif(return_num==2)
                if(IsReturn1)
                    E1_data2=[E1_data2;[t(end),E1]];
                end
                if(IsReturn2)
                    E2_data2=[E2_data2;[t(end),E2]];
                end
                E_data2=[E_data2;[t(end),E]];
            end
            %             if(return_num==0)
            %                 plot(t(end)/T,(E1+Ip)/omg,'.','Color', [1,0,0])
            %                 hold on;
            %             elseif(return_num==1)
            %                 plot(t(end)/T,(E1+Ip)/omg,'.','Color', [0,1,0])
            %                 hold on;
            %             else
            %                 plot(t(end)/T,(E1+Ip)/omg,'.','Color', [0,0,1])
            %                 hold on;
            %             end
            return_num=return_num+1;
        end

    end

%     

end

%% plot eneygy
% % 
% hold on
% plot(E_data(:,1)/T,(E_data(:,2)+Ip)/omg,'r.')
% hold on;
% plot(E1_data(:,1)/T-0.25,(E1_data(:,2)+Ip)/omg,'b.')
% hold on
% plot(E2_data(:,1)/T,(E2_data(:,2)+Ip)/omg,'g.')
% % %return-2
% hold on;
% plot(E_data1(:,1)/T,(E_data1(:,2)+Ip)/omg,'r.')
% hold on;
% plot(E1_data1(:,1)/T-0.25,(E1_data1(:,2)+Ip)/omg,'b.')
% hold on
% plot(E2_data1(:,1)/T,(E2_data1(:,2)+Ip)/omg,'g.')
%return>2

% hold on
% plot(E_data2(:,1)/T,(E_data2(:,2)+Ip)/omg,'r.')
% hold on;
% plot(E1_data2(:,1)/T,(E1_data2(:,2)+Ip)/omg,'b.')
% hold on
% plot(E2_data2(:,1)/T,(E2_data2(:,2)+Ip)/omg,'g.')

hold on;
plot(E_data(:,1)/T,(E_data(:,2)+Ip)/omg,'r.')
hold on;


%% electron1
figure
plot(t/T0,y(:,3),'b-',MarkerSize=5)
xlabel('times(o.c.)')
ylabel('z(t)')
title('return cycle=2')
legend('z(t)' ,'laser')
%% eletron2
figure
plot(t/T0,y(:,9),'b-',MarkerSize=5)
ff=exp(-2*log(2)*t.^2/Constants.tao^2);
xlabel('times(o.c.)')
ylabel('y(t)')
title('return cycle=2')