function [ ] = threeBodySun(  )
%Simulates the movement of the eart, moon and sun through solving the
%three-body-problem.


%Preparing the system
clc;
hold off;

%Animation-parameter
fps=25; %Frames per second
days= 8;%n days in realtime are equivalent to ...
seconds=1;%... m seconds in animationtime

%masses of the bodys (Unit: kg)
mEarth=5.974*10^24; %mass of the earth
mMoon=7.349*10^22; %mass of the moon
mSun=1.989*10^30; %mass of the sun

%Time for the calculation (Unit: s)
TSTART=0;           %Time at the beginning of the calculation
TEND=365.256*24*60*60; %Time at the end of the calculation
tvalues=TSTART:days*((60*24*60)/(fps*seconds)):TEND; %Array which contains 
%the points of time on which the solution of the calculation is needed

%Initial conditions 
%the values of xPos and yPos are the coordinates of the plantets. The unit
%is meter. The Values of xVel and yVel are the velocities of the planets in
%the x and y direction. The unit is meters per second.
xPosEarth= -149600000000;  %Index No.1 in Array y of the solution
xVelEarth= 0;              %Index No.2 in Array y of the solution
xPosMoon= -149194500000;   %Index No.3 in Array y of the solution
xVelMoon= 0;               %Index No.4 in Array y of the solution
xPosSun= 0;                %Index No.5 in Array y of the solution
xVelSun= 0;                %Index No.6 in Array y of the solution
yPosEarth= 0;              %Index No.7 in Array y of the solution
yVelEarth= 29780;          %Index No.8 in Array y of the solution
yPosMoon= 0;               %Index No.9 in Array y of the solution
yVelMoon= 29780+1023;      %Index No.10 in Array y of the solution
yPosSun= 0;                %Index No.11 in Array y of the solution
yVelSun= 0;                %Index No.12 in Array y of the solution

%Settings of the precision of the calculation
opt=odeset('AbsTol', 1e-6, 'RelTol', 1e-6, 'MaxStep', 20000);

%Solving of the first-order system of differential equations
[t,y]=ode45(@(t,y) sysDgl3body(t,y),tvalues,[xPosEarth; xVelEarth; xPosMoon; 
    xVelMoon; xPosSun; xVelSun; yPosEarth; yVelEarth; yPosMoon; yVelMoon; yPosSun; yVelSun],opt);

%Painting of the orbits of the planets
%Orbit of the earth
hdl=plot(y(:,1),y(:,7));
set(hdl,'Color','blue');
hold on
%Orbit of the moon
hdl=plot(y(:,3),y(:,9));
set(hdl,'Color','red');
%Orbit of the sun
hdl=plot(y(:,5),y(:,11));
set(hdl,'Color','green');
%Settings of the rendering
grid on;
title('3 body problem for Earth(blue);Moon(red);Sun(green)');
xmin=min([min(y(:,1)),min(y(:,3)),min(y(:,5))]);
xmax=max([max(y(:,1)),max(y(:,3)),max(y(:,5))]);
ymin=min([min(y(:,7)),min(y(:,9)),min(y(:,11))]);
ymax=max([max(y(:,7)),max(y(:,9)),max(y(:,11))]);
axis([xmin,xmax,ymin,ymax]);
axis equal;

%Counts the number of times the moon orbits the earth
count=0;    %contains the number of orbits
closer=true;    %flag determines which planet (moon or earth) is closer to 
%the sun: true: moon is closer to the sun
for i=1:size(t,1)
    distSunMoon=norm([y(i,3)-y(i,5),y(i,9)-y(i,11)]); %distance between 
    %sun and moon in meters
    distSunEarth=norm([y(i,1)-y(i,5),y(i,7)-y(i,11)]); %distance between 
    %sun and earth in meters
    if(~closer&&distSunMoon<distSunEarth)
        closer=true;
        count=count+1;
    else
        if(closer&&distSunMoon>distSunEarth)
            closer=false;
        end
    end
end

%Printing the amount of orbits
out=['The moon has orbited the earth ', num2str(count), ' times!'];
text(0.25*10^11,0.25*10^11,out);

%Animation of the movement of the planets
dothdl=plot([xPosEarth,xPosMoon,xPosSun],[yPosEarth,yPosMoon,yPosSun],'o');
animationTime=tic;
for i=1:size(t,1)
    %Paints the position of the planets at the point of time t(i)
    set(dothdl,'XData', [y(i,1),y(i,3),y(i,5)],'YData', [y(i,7),y(i,9),y(i,11)]);
    waitForSysTime(t(i,1));
    drawnow;
end

    function waitForSysTime(realTime)
        %Function, which forces the programm to wait for the right time to
        %paint the next frame.
        while realTime>toc(animationTime)*days*((60*60*24)/seconds)
        end
    end

    function dy = sysDgl3body( t,y )
        %Function which represents the first-order system of
        %differential equations.
        
        %Gravitational constant (in m^3/(kg*s^2))
        G=6.67384*10^(-11);
        
        %Calculation of radii (distance between two planets, Unit: m/s^2)
        rEarthMoon=sqrt((y(3)-y(1))^2+(y(9)-y(7))^2);
        rMoonSun=sqrt((y(3)-y(5))^2+(y(9)-y(11))^2);
        rEarthSun=sqrt((y(5)-y(1))^2+(y(7)-y(11))^2);
        
        %Prevents the System from the occurance of division by zero by
        %setting the acceleration instantly to zero whithout calculating.
        if rEarthMoon==0
            x1a1=0;     %1.acceleration of the earth in X-direction
            x2a1=0;     %1.acceleration of the moon in X-direction
            y1a1=0;     %1.acceleration of the earth in Y-direction
            y2a1=0;     %1.acceleration of the moon in Y-direction
        else
            x1a1=G*(mMoon/(rEarthMoon^3))*(y(3)-y(1));     %1.acceleration of the earth in X-direction
            x2a1=G*(mEarth/(rEarthMoon^3))*(y(1)-y(3));     %1.acceleration of the moon in X-direction
            y1a1=G*(mMoon/(rEarthMoon^3))*(y(9)-y(7));     %1.acceleration of the earth in Y-direction
            y2a1=G*(mEarth/(rEarthMoon^3))*(y(7)-y(9));     %1.acceleration of the moon in Y-direction
            
        end
        if rMoonSun==0
            x2a2=0;     %2.acceleration of the moon in X-direction
            x3a2=0;     %2.acceleration of the sun in X-direction
            y2a2=0;     %2.acceleration of the moon in Y-direction
            y3a2=0;     %2.acceleration of the sun in Y-direction
        else
            x2a2=G*(mSun/(rMoonSun^3))*(y(5)-y(3));     %2.acceleration of the moon in X-direction
            x3a2=G*(mMoon/(rMoonSun^3))*(y(3)-y(5));     %2.acceleration of the sun in X-direction
            y2a2=G*(mSun/(rMoonSun^3))*(y(11)-y(9));    %2.acceleration of the moon in Y-direction
            y3a2=G*(mMoon/(rMoonSun^3))*(y(9)-y(11));    %2.acceleration of the sun in Y-direction
        end
        if rEarthSun==0
            x1a2=0;     %2.acceleration of the earth in X-direction
            x3a1=0;     %1.acceleration of the sun in X-direction
            y1a2=0;     %2.acceleration of the earth in Y-direction
            y3a1=0;     %1.acceleration of the sun in Y-direction
        else
            x1a2=G*(mSun/(rEarthSun^3))*(y(5)-y(1));    %2.acceleration of the earth in X-direction
            x3a1=G*(mEarth/(rEarthSun^3))*(y(1)-y(5));    %1.acceleration of the sun in X-direction
            y1a2=G*(mSun/(rEarthSun^3))*(y(11)-y(7));   %2.acceleration of the earth in Y-direction
            y3a1=G*(mEarth/(rEarthSun^3))*(y(7)-y(11));   %1.acceleration of the sun in Y-direction
        end
        
        %The first-order system of differential equations
        dy=zeros(12,1);
        dy(1)=y(2);         %v(t) of the earth in X-direction
        dy(2)=x1a1+x1a2;    %a(t) of the earth in X-direction
        dy(3)=y(4);         %v(t) of the moon in X-direction
        dy(4)=x2a1+x2a2;    %a(t) of the moon in X-direction
        dy(5)=y(6);         %v(t) of the sun in X-direction
        dy(6)=x3a1+x3a2;    %a(t) of the sun in X-direction
        dy(7)=y(8);         %v(t) of the earth in Y-direction
        dy(8)=y1a1+y1a2;    %a(t) of the earth in Y-direction
        dy(9)=y(10);        %v(t) of the moon in Y-direction
        dy(10)=y2a1+y2a2;   %a(t) of the moon in Y-direction
        dy(11)=y(12);       %v(t) of the sun in Y-direction
        dy(12)=y3a1+y3a2;   %a(t) of the sun in Y-direction
    end
end