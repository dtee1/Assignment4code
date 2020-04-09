% David Talson
% ELEC4700 Assignment 3
clear all
clc

nx=200;
ny=100;
CurrDen = [];
G=sparse(nx*ny,nx*ny);
B=zeros(nx*ny,1);
x2=120;
x1=80;
y2=60;
y1=40;
condRegion=1;
condBottle=0.01;
cond=1.*ones(nx,ny);
count=1;

for i = 1:nx
    for j = 1:ny
        if(((i>=x1)&&(i<=x2)&&(j<=y1))||((i>=x1)&&(i<=x2)&&(j>=y2)))
            cond(i,j) = condBottle;
        end
    end
end
for k=0.1:10
    for i=1:nx
        for j=1:ny

            n = j + (i-1)*ny;
            nxm = (i-2)*ny + j;
            nxp = i*ny + j;
            nym = (i-1)*ny + j-1;
            nyp = (i-1)*ny + j+1;

            if i==1
                B(n,1)=k;
                G(n,n)=1;
            elseif i==nx
                B(n,1)=0;
                G(n,n)=1;
            elseif j==1
                B(n,1)=0;

                condxm = (cond(i,j) + cond(i-1,j))/2;
                condxp = (cond(i,j) + cond(i+1,j))/2;
                condyp = (cond(i,j) + cond(i,j+1))/2;

                G(n,n) = -(condxm+condxp+condyp);
                G(n,nxm) = condxm;
                G(n,nxp) = condxp;
                G(n,nyp) = condyp;

            elseif j==ny
                B(n,1)=0;

                condxm = (cond(i,j) + cond(i-1,j))/2;
                condxp = (cond(i,j) + cond(i+1,j))/2;
                condym = (cond(i,j) + cond(i,j-1))/2;

                G(n,n) = -(condxm+condxp+condym);
                G(n,nxm) = condxm;
                G(n,nxp) = condxp;
                G(n,nym) = condym;
            else
                B(n,1) = 0;
                condxm = (cond(i,j) + cond(i-1,j))/2;
                condxp = (cond(i,j) + cond(i+1,j))/2;
                condyp = (cond(i,j) + cond(i,j+1))/2;
                condym = (cond(i,j) + cond(i,j-1))/2;

                G(n,n) = -(condxm+condxp+condyp+condym);
                G(n,nxm) = condxm;
                G(n,nxp) = condxp;
                G(n,nym) = condym;
                G(n,nyp) = condyp;

            end
        end  
    end
    V=G\B;
    for i = 1:nx
        for j = 1:ny
            n = j+(i-1)*ny;
            map(i,j) = V(n,1);
        end
    end

    [Ex,Ey] = gradient(map',1,1);

    global C

    addpath ../geom2d/geom2d

    C.q_0 = 1.60217653e-19;             % electron charge
    C.hb = 1.054571596e-34;             % Dirac constant
    C.h = C.hb * 2 * pi;                    % Planck constant
    C.m_0 = 9.10938215e-31;             % electron mass
    C.kb = 1.3806504e-23;               % Boltzmann constant
    C.eps_0 = 8.854187817e-12;          % vacuum permittivity
    C.mu_0 = 1.2566370614e-6;           % vacuum permeability
    C.c = 299792458;                    % speed of light
    C.g = 9.80665;                      % metres (32.1740 ft) per s²
    C.m_eff = 0.26*C.m_0;               % Effective mass of electrons 

    regionWidth = 2e-7;                 % Nominal width of the region
    regionLength = 1e-7;                % Nominal length of the region
    T = 300;                            % Assumed temperature

    vth = sqrt((2*C.kb*T)/C.m_eff); % Thermal velocity as mean of magnitude of velocity
    tmn = 2e-13;                        % Mean time between collisions
    freepath = vth*2e-13;

    numElectrons = 1000;
    electronXpos = rand(1,numElectrons).*regionWidth;
    electronYpos = rand(1,numElectrons).*regionLength;

    % Check if there are electrons in the box
    outbox = 1;
    for chk = 1:numElectrons

             while (electronXpos(chk) >= 0.8e-7 & electronXpos(chk) <= 1.2e-7) & ((electronYpos(chk) >= 0 & electronYpos(chk) <= 0.4e-7)|(electronYpos(chk) >= 0.6e-7 & electronYpos(chk) <= 1e-7))
                  electronXpos(chk) = rand(1, 1).* regionWidth;
                  electronYpos(chk) = rand(1, 1).* regionLength;
             end

    end

    angle = rand(1,numElectrons).*2*pi;

    sig = sqrt(C.kb*T/C.m_eff)/4;
    MBdist = makedist('Normal',vth,sig);
    electronVel = random(MBdist,1,numElectrons);
    angle = rand(1,numElectrons).*2*pi;
    electronXvel = electronVel.*cos(angle);
    electronYvel = electronVel.*sin(angle);
    deltaT = 1e-9/vth;
    probScat = 1 - exp(-deltaT/tmn);
    electronVel = sqrt(electronXvel.^2 + electronYvel.^2);
    electronAccX =0;
    electronAccY =0;
    xposNew=electronXpos;
    yposNew=electronYpos;
    

        for i=1:1000
              for j = 1:numElectrons 
                  if probScat > rand
                        angleNew = rand(1).*2*pi;
                        electronXvel(1,j) = random(MBdist,1).*cos(angleNew);
                            electronYvel(1,j) = random(MBdist,1).*sin(angleNew);
                  end

                  if electronYpos(1,j) + electronYvel(1,j).*deltaT >=1e-7|| electronYpos(1,j) + electronYvel(1,j).*deltaT <= 0||((electronYpos(1,j) + electronYvel(1,j).*deltaT >= 0.6e-7||electronYpos(1,j) + electronYvel(1,j).*deltaT <= 0.4e-7) & (electronXpos(1,j) + electronXvel(1,j)*deltaT>=0.8e-7 & electronXpos(1,j) + electronXvel(1,j)*deltaT<=1.2e-7))
                        electronYvel(1,j) = electronYvel(1,j)*-1;          
                  end              
                  yposNew(i,j) = electronYpos(1,j) + electronYvel(1,j).*deltaT;
                  if electronXpos(1,j) + electronXvel(1,j)*deltaT >= regionWidth || electronXpos(1,j) + electronXvel(1,j)*deltaT <= 0||((electronYpos(1,j) + electronYvel(1,j).*deltaT >= 0.6e-7||electronYpos(1,j) + electronYvel(1,j).*deltaT <= 0.4e-7) & (electronXpos(1,j) + electronXvel(1,j)*deltaT>=0.8e-7 & electronXpos(1,j) + electronXvel(1,j)*deltaT<=1.2e-7))
                    if electronXpos(1,j) + electronXvel(1,j)*deltaT >= regionWidth ||  (electronXpos(1,j) + electronXvel(1,j)*deltaT >= 0.8e-7 & electronXpos(1,j) + electronXvel(1,j)*deltaT >= 1e-7)
                           xposNew(i,j) = 0;
                    else
                           xposNew(i,j) = regionWidth;
                    end
                  else
                         electronXvel(1,j) = electronXvel(1,j) + electronAccX*(deltaT);
                         xposNew(i,j) = electronXpos(1,j) + electronXvel(1,j).*deltaT + 0.5*electronAccX*(deltaT^2);
                  end
                        electronXvel(1,j) = electronXvel(1,j) + electronAccX*(deltaT);
                        xposNew(i,j) = electronXpos(1,j) + electronXvel(1,j).*deltaT + 0.5*electronAccX*(deltaT^2);

                      xii = ceil((1e7)*xposNew(i,j));
                      if xii<=0
                          xii=1;
                      end
                      if xii>200
                          xii=200;
                      end

                      yii=ceil((1e7)*yposNew(i,j));
                      if yii<=0
                          yii=1;
                      end
                      if yii>100
                          yii=100;
                      end
                      EX=(10^9)*Ex(yii,xii);
                      EY=(10^9)*Ey(yii,xii);
                      electronAccX = ((EX)*C.q_0)/C.m_0;
                      electronAccY = ((EY)*C.q_0)/C.m_0;      

                end

            electronXpos = xposNew(i,:);
            electronYpos = yposNew(i,:);
            t1 = sqrt(electronXvel(1,:).^2);
            t2(i) = 1000*(mean(t1)*C.q_0*(10e15)*regionWidth);
        end
   current(count) = mean(t2(i));
   
   voltage(count)=k;
   count=count+1;
end

y1 = polyfit(voltage,current,1);
resistance = 1/y1(1)
x1=linspace(0.1,10);
p=polyval(y1,x1);
str = sprintf('y = %f * x + %f', y1(1), y1(2));

figure(1)
plot(current,voltage,'o')
hold on
plot(p,x1,'LineWidth',2)
hold off
title('Current vs Voltage plot')
xlabel('Current')
ylabel('Voltage')
grid on

figure(2)
plot(t2,'LineWidth',2)
grid on 
xlabel('Simulation time')
ylabel('Current')
title('Current in x direction')