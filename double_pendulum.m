function y = double_pendulum(ivp, duration, fps, movie)

nframes=duration*fps;
h = 1/fps ;
global fixed
t = linspace(0,duration,nframes);
y = zeros(4,nframes-1) ;
y(:,1) = ivp ;
for i = 1:nframes-1
    X = y(:,i) ;
    h1 = RKfun(X) ;
    h2 = RKfun(X+h1*h/2) ;
    h3 = RKfun(X+h2*h/2) ;
    h4 = RKfun(X+h3*h) ;
    y(:,i+1) = X+(h1+h2*2+h3*2+h4)*h/6 ;
end

if movie==true
    phi1=y(1,:)'; dtphi1=y(2,:)';
    phi2=y(3,:)'; dtphi2=y(4,:)';
    l1=fixed(3); l2=fixed(4);
    
    h=plot(0,0,'MarkerSize',30,'Marker','.','LineWidth',2);
    range=1.1*(l1+l2); axis([-range range -range range]); axis square;
    set(gca,'nextplot','replacechildren');
    v = VideoWriter('../results/doublePendulum.mp4','MPEG-4');
    open(v)
    tt = 1 ;
    for i=1:length(phi1)-1
        if (ishandle(h)==1)
            Xcoord=[0,l1*sin(phi1(i)),l1*sin(phi1(i))+l2*sin(phi2(i))];
            Ycoord=[0,-l1*cos(phi1(i)),-l1*cos(phi1(i))-l2*cos(phi2(i))];
            set(h,'XData',Xcoord,'YData',Ycoord);
            drawnow;
            title(['frame = ',num2str(i),'(',num2str(fps),' Hz)']) ;
            mov(tt) = getframe(gcf);
            tt = tt + 1 ;
            if 0
                pause(t(i+1)-t(i));
            end
        end
    end
    writeVideo(v,mov)
    close(v)
end
%% subfunction
    function Y = RKfun(x)
        g=fixed(1); m1=fixed(2); m2=fixed(3); l1=fixed(4); l2=fixed(5);
        m12 = m1+m2 ;
        S = sin(x(1)-x(3)) ;
        C = cos(x(1)-x(3)) ;
        s1 = sin(x(1)) ;
        s2 = sin(x(3)) ;
        dt_th1 = x(2) ;
        dt_th2 = x(4) ;
        
        Y = zeros(4,1);
        Y(1) = dt_th1 ;
        Y(2) = (-g*m12*s1 + m2*g*s2*C - m2*l2*dt_th2^2*S - m2*l1*dt_th1^2*C*S )...
            / (l1*(m12-m2*C^2)) ; 
        Y(3) = dt_th2 ;
        Y(4) = (g*m12*C*s1 - g*m12*s2 + m12*l1*dt_th1^2*S + m2*l2*dt_th2^2*C*S )...
            / (l2*(m12-m2*C^2)); 
        if  m1 == m2 && l1 == l2
            l = l1 ;
            Y(2) = ( -2*g*s1 + g*s2*C - l*dt_th2^2*S - l*dt_th1^2*C*S )...
                / (l*(2-C^2)); 
            Y(4) = (2*g*C*s1 - 2*g*s2 + 2*l*dt_th1^2*S + l*dt_th2^2*C*S )...
                / (l*(2-C^2)); 
        end
    end
end