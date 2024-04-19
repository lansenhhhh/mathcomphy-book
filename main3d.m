E0=0.05337605126836238;
omg=0.0569541683375;
n1=2000;
n2=1000;
rcycle=2;
ncycle=20;
Up=(E0^2)/(4*omg^2);
Ip=0.90357;
T=2*pi/omg;
dt=rcycle*T/n2;
a=zeros(1,6);
counter=zeros(n1,1);
% 获取当前坐标轴的限制
xlims = get(gca, 'XLim');
ylims = get(gca, 'YLim');
x0_range =40;
cmap = [linspace(0, 1, length(x0_range))', ...
    linspace(0, 1, length(x0_range))', ...
    linspace(0, 1, length(x0_range))'];
hold on;
figure
for k = 1:length(x0_range)
    for j=1:n1
        xarr=zeros(n1,6);
        t0=j/n1*4*T+8*T;
        return_num=0;
        x0 = x0_range(k);
        x1=[-x0,0.2,0];
        x2=[x0,0.2,0];
        x=[x1,x2];
        v0=-0.4;
        v1=[0  0 -v0];
        v2=[0 0 v0];
        r1=sqrt(x1*x1');
        r2=sqrt(x2*x2');
        r12=(x1-x2)*(x1-x2)';
        v=[v1,v2];
        for i=1:n2
            t1=t0+i*dt;
            a1=[0,0,-E0*sin(omg*t1)*(sin((pi*t1)./(ncycle*2*pi/omg))).^2]-x1*2*(r1^2.0)^(-1.5)+(x1-x2)./(r12.^2);
            a2=[0,0,-E0*sin(omg*t1)*(sin((pi*t1)./(ncycle*2*pi/omg))).^2]-x2*2*(r2^2.0)^(-1.5)-(x1-x2)./(r12.^2);
            a=[a1,a2];
            v=v+a*dt;
            x_iter1=x(3);
            x_iter2=x(6);
            x=x+v*dt;
            x1=x(1:3);
            x2=x(4:6);
            r12=(x1-x2)*(x1-x2)';
            r1=sqrt(x1*x1');
            r2=sqrt(x2*x2');
            E1=(v(1:3)*v(1:3)')/2.0-2/r1;
            E2=(v(4:6)*v(4:6)')/2.0-2/r2;
            E=(v(1:6)*v(1:6)')/2.0-2/r1-2/r2+1/r12;
            ts=t1-t0;
            xarr(i,:)=x;
            if ((x_iter2*x(6)<0)||(x_iter1*x(3)<0))%
                if(return_num==0)
                    plot(t1/T,(E+Ip)/omg,'.','Color', [1,0,0])
                    hold on;
                elseif(return_num==1)
                    plot(t1/T,(E+Ip)/omg,'.','Color', [0,1,0])
                    hold on;
                else
                    plot(t1/T,(E+Ip)/omg,'.','Color', [0,0,1])
                    hold on;
                end
                return_num=return_num+1;
                counter(j)=1;

            end
            %         if counter(j)==1
            %         plot3(xarr(:,1),xarr(:,2),xarr(:,3),'b.');
            %         hold on;
            %         end
        end
    end
    
hold on;
end
xlabel("Times(o.c.)")
ylabel("Harmonic Order")