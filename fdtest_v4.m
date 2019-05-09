% fdtest_v4 Exp 
% a script to try out a 2D explicit method, in factored form
% March 7, 2019
% Stablility dependant on epsilon
% R. Vestrum
% 


close all;

N=100; % number of grid points in each direction
Nm1 = N-1; % useful constants
Nm2 = N-2; 
tmax=1; %Max time to run
x=linspace(-2,2,N);
y=linspace(-2,2,N);
[X,Y] = meshgrid(x,y);

v = 3; % wavespped
dx = x(2)-x(1);

dy = y(2)-y(1);

dt_factor=[1,2,4,8,12,16];
nm=size(dt_factor,2); %Number of times to run code
for m=1:nm
    if m==1
        dt(m)=floor(dt_factor(m)*min(dx,dy)/(16*v)*10000)/10000;
    else
        dt(m)=dt(1)*dt_factor(m);
    end
    ex2(m) = (v*dt(m)/dx)^2;
    ey2(m) = (v*dt(m)/dy)^2;
    thetax(m) = asin(-2*ex2(m)/(1+2*ex2(m)))/2;
    a(m) = sqrt(1 + 2*ex2(m))*cos(thetax(m));
    b(m) = sqrt(1 + 2*ex2(m))*sin(thetax(m));
    thetay(m) = asin(-2*ey2(m)/(1+2*ey2(m)))/2;
    c(m) = sqrt(1 + 2*ey2(m))*cos(thetay(m));
    d(m) = sqrt(1 + 2*ey2(m))*sin(thetay(m));
end

% lets make u to hold k-1, k and k+1 time indices
u = zeros(N,N,3,nm);
lhs = zeros(N,N); % left hand side
mhs = zeros(N,N); % middle hand side
rhs = zeros(N,N); % right hand side
z = zeros(N,N,nm); % for plotting

% a one point initial condition
u(N/2,N/2,:,:) = 1;

%%
% a Gaussian initial condition
u(:,:,1,1) = exp(-(X.*X + Y.*Y)/.01);
for i=1:3
    for j=1:nm
        u(:,:,i,j) = u(:,:,1,1);
    end
end
u_e=u;

%%
% we need a cyclic counter for indexing into u time steps. 
cyc = zeros(nm);

t=0;

%% Creates plots for later

if nm == 4
    i=2;
    j=2;
else
    i=2;
    j=3;
end
for p=1:2
    pl(p)=figure;
    for m=1:nm
        pax(p,m)=subplot(i,j,m);
    end
end


%% Starts running from t=0 to t=tmax

k = zeros(nm);
km1 = zeros(nm);
kp1 = zeros(nm);

step=0;

while t < tmax  % start the timing loop 
t = t + dt(1);
step=step+1;
tic


%
% Each dt is a multiple of the first dt. This runs over each dt_factor and
% checks if the current step is a multiple of the dt_factor. If so it
% calculates for that dt.
%
for m=1:nm %for each dt_factor
    if mod(step,dt_factor(m))==0 %If current step is a multiple of dt_factor 
        update(m)=1;
        cyc(m) = mod(cyc(m) + 1,3);
        k(m) = 1 + mod(cyc(m),3);
        km1(m) = 1 + mod(cyc(m)-1,3);
        kp1(m) = 1 + mod(cyc(m)+1,3);

        rhs = 0*rhs; % reset to zero
        rhs(2:Nm1,2:Nm1) = 4*u(2:Nm1,2:Nm1,k(m),m) ...
          + ex2(m)*(u(3:N,2:Nm1,km1(m),m) + u(1:Nm2,2:Nm1,km1(m),m)) ...
          + ey2(m)*(u(2:Nm1,3:N,km1(m),m) + u(2:Nm1,1:Nm2,km1(m),m)) ...
          - 2*(1+ex2(m)+ey2(m))*u(2:Nm1,2:Nm1,km1(m),m);

        % solve [(a+bR) + 1i(c+dU)]mhs = rhs
        % mhs = 
        mhs = 0*mhs; % make sure edges have zeros
        for i=2:N
            for j=2:N
                mhs(i,j) = (rhs(i,j) - b(m)*mhs(i-1,j) - 1i*d(m)*mhs(i,j-1) )/(a(m)+1i*c(m));
            end
        end

        %  
        % solve [(a+bL) - 1i(c+dD)]lhs = mhs
        %
        lhs = 0*lhs;
        for i=N-1:-1:1
            for j=N-1:-1:1
                lhs(i,j) = (mhs(i,j) - b(m)*lhs(i+1,j) + 1i*d(m)*lhs(i,j+1) )/(a(m)-1i*c(m));
            end
        end

        %
        % now stick into u
        %
        u(2:Nm1,2:Nm1,kp1(m),m) = real(lhs(2:Nm1,2:Nm1));
        z(:,:,m) = (u(:,:,kp1(m),m));
    else
        update(m)=0;
    end
end

% Updates plots                             
for m=1:nm
    for p=1:2
        %
        % Updates first 4 plots if the 4th plot has changed. The
        % intermediate updates for the first 3 plots was small and this
        % speeds up the plotting. The last two plots are only updated when
        % the plot changes.
        %
        if (update(4)==1 && m<4) || (update(m)==1 && m>3)
            surf(pax(p,m),x,y,z(:,:,m),'linestyle','none')

            title_text=sprintf("time step: %1.1f ms   time: %1.1f ms\n $\\epsilon_x^2 = %5.5f$",dt(m)*1000,t*1000,ex2(m));
            title(pax(p,m),title_text,'interpreter','latex')
            if p==2
                view(pax(p,m),2)
                pax(p,m).NextPlot = 'add';
                contour3(pax(p,m),x,y,z(:,:,m),3,'k')
                pax(p,m).NextPlot = 'replace';

            end

        end
                
    end
end

pause(0.01) %Pauses so the computer can update the scaling, images, 
            %before resuming the computation.
  
end
 
