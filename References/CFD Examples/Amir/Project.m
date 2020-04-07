% CV, staggered mesh
clc

%% inputs
%lengths
w_up=.2; %width of upper inlet
w_left=.2; %width of left inlet
l_in=.2; %input length
h_in=.2; %input height
l_out=.6; %output length
l=l_in+w_up+l_out; %total length
h=h_in+w_left; %total height

%mass inlets
u_in_left=.2; %inlet from left boundary
v_in_left=0;
v_in_up=-.2; %inlet from upper boundary (negative because downward)
u_in_up=0;
m_in=w_left*u_in_left+w_up*(-v_in_up);

%concentration inlets
C_in_left=0; %concentration from left boundary
C_in_up=1; %concentration from upper boundary

%fluid properties
rho=1000; %density
nu=9e-4; %kinematic viscosity
D=2e-4; %diffusivity
P_0=0; %atmospheric pressure

%% discretization
m=30; %number of cells in a row (small for fast solving, in the report 60 and 90 are used for this geometry)
n=round(m*h/l); %number of cells in a column (to get dx=dy, can set it to any value if dx~=dy)
dx=l/m;
dy=h/n;
dt=0.01 * min( (min(dx,dy))^2/nu , min(dx,dy)/(u_in_left-v_in_up) );
dt=round(dt,4);
t_end=1000;

%sigma
sigma_nu=dt/((min(dx,dy))^2/nu);
sigma_D=dt/((min(dx,dy))^2/D);
%check stability
if max(sigma_nu,sigma_D)>0.1
    'dt too high (sigma)'
end
%to check steady state
epsilon_u=1e-14;
epsilon_v=1e-14;
epsilon_P=1e-10;
epsilon_C=1e-12;

%% problem matrices

%unknowns
P=P_0*ones(m,n); %pressure
C=zeros(m,n); %concentration
u=zeros(m+1,n); %u-velocity
du=zeros(size(u));
v=zeros(m,n+1); %v-velocity
dv=zeros(size(v));
P_prime=zeros(size(P));
u_prime=zeros(size(u));
v_prime=zeros(size(v));

%out of channel
P([1:round(l_in/dx),round((l_in+w_up)/dx)+1:end],round(w_left/dy)+1:end)=NaN;
C([1:round(l_in/dx),round((l_in+w_up)/dx)+1:end],round(w_left/dy)+1:end)=NaN;
u([1:round(l_in/dx),round((l_in+w_up)/dx)+2:end],round(w_left/dy)+1:end)=NaN;
v([1:round(l_in/dx),round((l_in+w_up)/dx)+1:end],round(w_left/dy)+2:end)=NaN;

%% determine type of cell (mask)
%P
cell_P=zeros(size(P));
cell_P([1:round(l_in/dx),round((l_in+w_up)/dx)+1:end],round(w_left/dy)+1:end)=NaN;
cell_P(2:end-1,1)=1; %lower wall
cell_P([2:round(l_in/dx),round((l_in+w_up)/dx)+1:end-1],round(w_left/dy))=2; %upper wall
cell_P(round(l_in/dx)+1,round(w_left/dy)+1:end-1)=3; %left wall
cell_P(round((l_in+w_up)/dx),round(w_left/dy)+1:end-1)=4; %right wall
cell_P(1,2:round(w_left/dy)-1)=5; %left inlet
cell_P(round(l_in/dx)+2:round((l_in+w_up)/dx)-1,end)=6; %upper inlet
cell_P(end,2:round(w_left/dy)-1)=7; %outlet
cell_P(1,1)=8; cell_P(end,1)=9; cell_P(1,round(w_left/dy))=10; cell_P(end,round(w_left/dy))=11;
cell_P(round(l_in/dx)+1,end)=12; cell_P(round((l_in+w_up)/dx),end)=13; %corners

%u
cell_u=zeros(size(u));
cell_u([1:round(l_in/dx),round((l_in+w_up)/dx)+2:end],round(w_left/dy)+1:end)=NaN;
cell_u(2:end-1,1)=1; %lower wall
cell_u([2:round(l_in/dx),round((l_in+w_up)/dx)+2:end-1],round(w_left/dy))=2; %upper wall
cell_u(round(l_in/dx)+1,round(w_left/dy)+1:end-1)=3; %left wall
cell_u(round((l_in+w_up)/dx)+1,round(w_left/dy)+1:end-1)=4; %right wall
cell_u(1,2:round(w_left/dy)-1)=5; %left inlet
cell_u(round(l_in/dx)+2:round((l_in+w_up)/dx),end)=6; %upper inlet
cell_u(end,2:round(w_left/dy)-1)=7; %outlet
cell_u(1,1)=8; cell_u(end,1)=9; cell_u(1,round(w_left/dy))=10; cell_u(end,round(w_left/dy))=11;
cell_u(round(l_in/dx)+1,end)=12; cell_u(round((l_in+w_up)/dx)+1,end)=13; %corners

%v
cell_v=zeros(size(v));
cell_v([1:round(l_in/dx),round((l_in+w_up)/dx)+1:end],round(w_left/dy)+2:end)=NaN;
cell_v(2:end-1,1)=1; %lower wall
cell_v([2:round(l_in/dx),round((l_in+w_up)/dx)+1:end-1],round(w_left/dy)+1)=2; %upper wall
cell_v(round(l_in/dx)+1,round(w_left/dy)+2:end-1)=3; %left wall
cell_v(round((l_in+w_up)/dx),round(w_left/dy)+2:end-1)=4; %right wall
cell_v(1,2:round(w_left/dy))=5; %left inlet
cell_v(round(l_in/dx)+2:round((l_in+w_up)/dx)-1,end)=6; %upper inlet
cell_v(end,2:round(w_left/dy))=7; %outlet
cell_v(1,1)=8; cell_v(end,1)=9; cell_v(1,round(w_left/dy)+1)=10; cell_v(end,round(w_left/dy)+1)=11;
cell_v(round(l_in/dx)+1,end)=12; cell_v(round((l_in+w_up)/dx),end)=13; %corners

%% for solving P correction system
cell_number=zeros(size(P));

%numbering P_cells
k=0;
for i=1:m
    for j=1:n
        if ~isnan(cell_P(i,j))
            cell_number(i,j)=k+1;
            k=k+1;
        end
    end
end
A=zeros(k);
b=zeros(k,1);

%setting A matrix for P_prime solving (constant)
for i=1:m
    for j=1:n
        switch cell_P(i,j)
            case 0 %general
                A( cell_number(i,j) , cell_number(i,j) ) = -(2/dx^2+2/dy^2);
                A( cell_number(i,j) , cell_number(i+1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i-1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j+1) ) = 1/dy^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 1 %lower wall
                A( cell_number(i,j) , cell_number(i,j) ) = -(2/dx^2+1/dy^2);
                A( cell_number(i,j) , cell_number(i+1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i-1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j+1) ) = 1/dy^2;
            case 2 %upper wall
                A( cell_number(i,j) , cell_number(i,j) ) = -(2/dx^2+1/dy^2);
                A( cell_number(i,j) , cell_number(i+1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i-1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 3 %left wall
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+2/dy^2);
                A( cell_number(i,j) , cell_number(i+1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j+1) ) = 1/dy^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 4 %right wall
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+2/dy^2);
                A( cell_number(i,j) , cell_number(i-1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j+1) ) = 1/dy^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 5 %left inlet
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+2/dy^2);
                A( cell_number(i,j) , cell_number(i+1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j+1) ) = 1/dy^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 6 %upper inlet
                A( cell_number(i,j) , cell_number(i,j) ) = -(2/dx^2+1/dy^2);
                A( cell_number(i,j) , cell_number(i+1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i-1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 7 %outlet
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+2/dy^2);
                A( cell_number(i,j) , cell_number(i-1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j+1) ) = 1/dy^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 8 %corner 1
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+1/dy^2);
                A( cell_number(i,j) , cell_number(i+1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j+1) ) = 1/dy^2;
            case 9 %corner 2
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+1/dy^2);
                A( cell_number(i,j) , cell_number(i-1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j+1) ) = 1/dy^2;
            case 10 %corner 3
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+1/dy^2);
                A( cell_number(i,j) , cell_number(i+1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 11 %corner 4
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+1/dy^2);
                A( cell_number(i,j) , cell_number(i-1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 12 %corner 5
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+1/dy^2);
                A( cell_number(i,j) , cell_number(i+1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
            case 13 %corner 6
                A( cell_number(i,j) , cell_number(i,j) ) = -(1/dx^2+1/dy^2);
                A( cell_number(i,j) , cell_number(i-1,j) ) = 1/dx^2;
                A( cell_number(i,j) , cell_number(i,j-1) ) = 1/dy^2;
        end
    end
end
A = A(1:end-1,1:end-1); %since last P doesn't need correction

%% IC
u(cell_u==5)=u_in_left;
u(cell_u==8)=u_in_left;
u(cell_u==10)=u_in_left;
v(cell_v==6)=v_in_up;
v(cell_v==12)=v_in_up;
v(cell_v==13)=v_in_up;
% for i=1:m
%     for j=1:n
%         if ~isnan(cell_P(i,j)) && (j-1/2)*dy>w_left
%             C(i,j)=C_in_up;
%         end
%         if cell_u(i,j)==0 && (i-1)*dx<l_in
%             u(i,j)=u_in_left;
%         end
%         if cell_v(i,j)==0 && (j-1)*dy>w_left
%             v(i,j)=v_in_up;
%         end
%     end
% end

%% for saving
u_save=u;
v_save=v;
C_save=C;
P_save=P;
t_save=0;

%% iterations
for t=dt:dt:t_end
    %% u momentum
    du=zeros(size(u));
    for i=1:m+1 %on x
        for j=1:n %on y
            if ~isnan(cell_u(i,j))
                Fc=0; Fd=0; Qp=0;
                switch cell_u(i,j)
                    case 0 %general
                        Fc = dy/4*( (u(i+1,j)+u(i,j)) * (u(i+1,j)+u(i,j)) - (u(i-1,j)+u(i,j)) * (u(i-1,j)+u(i,j)) ) ...
                            +dx/4*( (v(i-1,j+1)+v(i,j+1)) * (u(i,j+1)+u(i,j)) - (v(i-1,j)+v(i,j)) * (u(i,j-1)+u(i,j)) );
                        Fd = nu*( dy/dx*( (u(i+1,j)-u(i,j)) - (u(i,j)-u(i-1,j)) ) + dx/dy*( (u(i,j+1)-u(i,j)) - (u(i,j)-u(i,j-1)) ) );
                        Qp = 1/rho*dy*(P(i,j)-P(i-1,j));
                    case 1 %lower wall
                        Fc = dy/4*( (u(i+1,j)+u(i,j)) * (u(i+1,j)+u(i,j)) - (u(i-1,j)+u(i,j)) * (u(i-1,j)+u(i,j)) ) ...
                            +dx/4*( (v(i-1,j+1)+v(i,j+1)) * (u(i,j+1)+u(i,j)) - 0 );
                        Fd = nu*( dy/dx*( (u(i+1,j)-u(i,j)) - (u(i,j)-u(i-1,j)) ) + dx/dy*( (u(i,j+1)-u(i,j)) - 1/3*(9*u(i,j)-u(i,j+1)) ) );
                        Qp = 1/rho*dy*(P(i,j)-P(i-1,j));
                    case 2 %upper wall
                        Fc = dy/4*( (u(i+1,j)+u(i,j)) * (u(i+1,j)+u(i,j)) - (u(i-1,j)+u(i,j)) * (u(i-1,j)+u(i,j)) ) ...
                            +dx/4*( 0 - (v(i-1,j)+v(i,j)) * (u(i,j-1)+u(i,j)) );
                        Fd = nu*( dy/dx*( (u(i+1,j)-u(i,j)) - (u(i,j)-u(i-1,j)) ) + dx/dy*( 1/3*(-9*u(i,j)+u(i,j-1)) - (u(i,j)-u(i,j-1)) ) );
                        Qp = 1/rho*dy*(P(i,j)-P(i-1,j));
                    case 6 %upper inlet
                        Fc = dy/4*( (u(i+1,j)+u(i,j)) * (u(i+1,j)+u(i,j)) - (u(i-1,j)+u(i,j)) * (u(i-1,j)+u(i,j)) ) ...
                            +dx/4*( (4*v_in_up*u_in_up) - (v(i-1,j)+v(i,j)) * (u(i,j-1)+u(i,j)) );
                        Fd = nu*( dy/dx*( (u(i+1,j)-u(i,j)) - (u(i,j)-u(i-1,j)) ) + dx/dy*( 1/3*(8*u_in_up-9*u(i,j)+u(i,j-1)) - (u(i,j)-u(i,j-1)) ) );
                        Qp = 1/rho*dy*(P(i,j)-P(i-1,j));    
                end
                du(i,j)=-Fc+Fd-Qp; %not units of u
            end
        end
    end
    
    %% v momentum
    dv=zeros(size(v));
    for i=1:m
        for j=1:n+1
            if ~isnan(cell_v(i,j))
                Fc=0; Fd=0; Qp=0;
                switch cell_v(i,j)
                    case 0 %general
                        Fc = dy/4*( (u(i+1,j)+u(i+1,j-1)) * (v(i+1,j)+v(i,j)) - (u(i,j)+u(i,j-1)) * (v(i-1,j)+v(i,j)) )...
                            +dx/4*( (v(i,j+1)+v(i,j)) * (v(i,j+1)+v(i,j)) - (v(i,j-1)+v(i,j)) * (v(i,j-1)+v(i,j)) );
                        Fd = nu*( dy/dx*( (v(i+1,j)-v(i,j)) - (v(i,j)-v(i-1,j)) ) + dx/dy*( (v(i,j+1)-v(i,j)) - (v(i,j)-v(i,j-1)) ) );
                        Qp=1/rho*dx*(P(i,j)-P(i,j-1));
                    case 3 %left wall
                        Fc = dy/4*( (u(i+1,j)+u(i+1,j-1)) * (v(i+1,j)+v(i,j)) - 0 )...
                            +dx/4*( (v(i,j+1)+v(i,j)) * (v(i,j+1)+v(i,j)) - (v(i,j-1)+v(i,j)) * (v(i,j-1)+v(i,j)) );
                        Fd = nu*( dy/dx*( (v(i+1,j)-v(i,j)) - 1/3*(9*v(i,j)-v(i+1,j)) ) + dx/dy*( (v(i,j+1)-v(i,j)) - (v(i,j)-v(i,j-1)) ) );
                        Qp=1/rho*dx*(P(i,j)-P(i,j-1));
                    case 4 %right wall
                        Fc = dy/4*( 0 - (u(i,j)+u(i,j-1)) * (v(i-1,j)+v(i,j)) )...
                            +dx/4*( (v(i,j+1)+v(i,j)) * (v(i,j+1)+v(i,j)) - (v(i,j-1)+v(i,j)) * (v(i,j-1)+v(i,j)) );
                        Fd = nu*( dy/dx*( 1/3*(-9*v(i,j)+v(i-1,j)) - (v(i,j)-v(i-1,j)) ) + dx/dy*( (v(i,j+1)-v(i,j)) - (v(i,j)-v(i,j-1)) ) );
                        Qp=1/rho*dx*(P(i,j)-P(i,j-1));
                    case 5 %left inlet
                        Fc = dy/4*( (u(i+1,j)+u(i+1,j-1)) * (v(i+1,j)+v(i,j)) - 4*u_in_left*v_in_left )...
                            +dx/4*( (v(i,j+1)+v(i,j)) * (v(i,j+1)+v(i,j)) - (v(i,j-1)+v(i,j)) * (v(i,j-1)+v(i,j)) );
                        Fd = nu*( dy/dx*( (v(i+1,j)-v(i,j)) - 1/3*(-8*v_in_left+9*v(i,j)-v(i+1,j)) ) + dx/dy*( (v(i,j+1)-v(i,j)) - (v(i,j)-v(i,j-1)) ) );
                        Qp=1/rho*dx*(P(i,j)-P(i,j-1));
                    case 7 %outlet
                        Fc = dy/4*( (u(i+1,j)+u(i+1,j-1)) * (2*v(i,j)) - (u(i,j)+u(i,j-1)) * (v(i-1,j)+v(i,j)) )...
                            +dx/4*( (v(i,j+1)+v(i,j)) * (v(i,j+1)+v(i,j)) - (v(i,j-1)+v(i,j)) * (v(i,j-1)+v(i,j)) );
                        Fd = nu*( dy/dx*( 0 - (v(i,j)-v(i-1,j)) ) + dx/dy*( (v(i,j+1)-v(i,j)) - (v(i,j)-v(i,j-1)) ) );
                        Qp=1/rho*dx*(P(i,j)-P(i,j-1));
                end
                dv(i,j)=-Fc+Fd-Qp; %not units of v
            end
        end
    end
    
    %% concentration equation
    dC=zeros(size(C));
    for i=1:m
        for j=1:n
            if ~isnan(cell_P(i,j))
                Fc=0; Fd=0;
                switch cell_P(i,j)
                    case 0 %general
                        Fc = dy/2*( u(i+1,j) * (C(i+1,j)+C(i,j)) - u(i,j) * (C(i,j)+C(i-1,j)) )...
                            +dx/2*( v(i,j+1) * (C(i,j+1)+C(i,j)) - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( (C(i+1,j)-C(i,j)) - (C(i,j)-C(i-1,j)) ) + dx/dy*( (C(i,j+1)-C(i,j)) - (C(i,j)-C(i,j-1)) ) );
                    case 1 %lower wall
                        Fc = dy/2*( u(i+1,j) * (C(i+1,j)+C(i,j)) - u(i,j) * (C(i,j)+C(i-1,j)) )...
                            +dx/2*( v(i,j+1) * (C(i,j+1)+C(i,j)) - 0 );
                        Fd = D*( dy/dx*( (C(i+1,j)-C(i,j)) - (C(i,j)-C(i-1,j)) ) + dx/dy*( (C(i,j+1)-C(i,j)) - 0 ) );
                    case 2 %upper wall
                        Fc = dy/2*( u(i+1,j) * (C(i+1,j)+C(i,j)) - u(i,j) * (C(i,j)+C(i-1,j)) )...
                            +dx/2*( 0 - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( (C(i+1,j)-C(i,j)) - (C(i,j)-C(i-1,j)) ) + dx/dy*( 0 - (C(i,j)-C(i,j-1)) ) );
                    case 3 %left wall
                        Fc = dy/2*( u(i+1,j) * (C(i+1,j)+C(i,j)) - 0 )...
                            +dx/2*( v(i,j+1) * (C(i,j+1)+C(i,j)) - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( (C(i+1,j)-C(i,j)) - 0 ) + dx/dy*( (C(i,j+1)-C(i,j)) - (C(i,j)-C(i,j-1)) ) );
                    case 4 %right wall
                        Fc = dy/2*( 0 - u(i,j) * (C(i,j)+C(i-1,j)) )...
                            +dx/2*( v(i,j+1) * (C(i,j+1)+C(i,j)) - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( 0 - (C(i,j)-C(i-1,j)) ) + dx/dy*( (C(i,j+1)-C(i,j)) - (C(i,j)-C(i,j-1)) ) );
                    case 5 %left inlet
                        Fc = dy/2*( u(i+1,j) * (C(i+1,j)+C(i,j)) - u_in_left * 2*C_in_left )...
                            +dx/2*( v(i,j+1) * (C(i,j+1)+C(i,j)) - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( (C(i+1,j)-C(i,j)) - 1/3*(-8*C_in_left+9*C(i,j)-C(i+1,j)) ) + dx/dy*( (C(i,j+1)-C(i,j)) - (C(i,j)-C(i,j-1)) ) );
                    case 6 %upper inlet
                        Fc = dy/2*( u(i+1,j) * (C(i+1,j)+C(i,j)) - u(i,j) * (C(i,j)+C(i-1,j)) )...
                            +dx/2*( v_in_up * 2*C_in_up - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( (C(i+1,j)-C(i,j)) - (C(i,j)-C(i-1,j)) ) + dx/dy*( 1/3*(8*C_in_up-9*C(i,j)-C(i,j-1)) - (C(i,j)-C(i,j-1)) ) );
                    case 7 %outlet
                        Fc = dy/2*( u(i+1,j) * (2*C(i,j)) - u(i,j) * (C(i,j)+C(i-1,j)) )...
                            +dx/2*( v(i,j+1) * (C(i,j+1)+C(i,j)) - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( 0 - (C(i,j)-C(i-1,j)) ) + dx/dy*( (C(i,j+1)-C(i,j)) - (C(i,j)-C(i,j-1)) ) );
                    case 8 %corner 1
                        Fc = dy/2*( u(i+1,j) * (C(i+1,j)+C(i,j)) - u_in_left * 2*C_in_left )...
                            +dx/2*( v(i,j+1) * (C(i,j+1)+C(i,j)) - 0 );
                        Fd = D*( dy/dx*( (C(i+1,j)-C(i,j)) - 1/3*(-8*C_in_left+9*C(i,j)-C(i+1,j)) ) + dx/dy*( (C(i,j+1)-C(i,j)) - 0 ) );
                    case 9 %corner 2
                        Fc = dy/2*( u(i+1,j) * (2*C(i,j)) - u(i,j) * (C(i,j)+C(i-1,j)) )...
                            +dx/2*( v(i,j+1) * (C(i,j+1)+C(i,j)) - 0 );
                        Fd = D*( dy/dx*( 0 - (C(i,j)-C(i-1,j)) ) + dx/dy*( (C(i,j+1)-C(i,j)) - 0 ) );
                    case 10 %corner 3
                        Fc = dy/2*( u(i+1,j) * (C(i+1,j)+C(i,j)) - u_in_left * 2*C_in_left )...
                            +dx/2*( 0 - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( (C(i+1,j)-C(i,j)) - 1/3*(-8*C_in_left+9*C(i,j)-C(i+1,j)) ) + dx/dy*( 0 - (C(i,j)-C(i,j-1)) ) );
                    case 11 %corner 4
                        Fc = dy/2*( u(i+1,j) * (2*C(i,j)) - u(i,j) * (C(i,j)+C(i-1,j)) )...
                            +dx/2*( 0 - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( 0 - (C(i,j)-C(i-1,j)) ) + dx/dy*( 0 - (C(i,j)-C(i,j-1)) ) );
                    case 12 %corner 5
                        Fc = dy/2*( u(i+1,j) * (C(i+1,j)+C(i,j)) - 0 )...
                            +dx/2*( v_in_up * 2*C_in_up - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( (C(i+1,j)-C(i,j)) - 0 ) + dx/dy*( 1/3*(8*C_in_up-9*C(i,j)-C(i,j-1)) - (C(i,j)-C(i,j-1)) ) );
                    case 13 %corner 6
                        Fc = dy/2*( 0 - u(i,j) * (C(i,j)+C(i-1,j)) )...
                            +dx/2*( v_in_up * 2*C_in_up - v(i,j) * (C(i,j)+C(i,j-1)) );
                        Fd = D*( dy/dx*( 0 - (C(i,j)-C(i-1,j)) ) + dx/dy*( 1/3*(8*C_in_up-9*C(i,j)-C(i,j-1)) - (C(i,j)-C(i,j-1)) ) );
                end
                dC(i,j)=-Fc+Fd; %not units of C
            end
        end
    end
    
    %% update uncorrected u and v
    u=u+du*dt/dx/dy;
    v=v+dv*dt/dx/dy;
    C=C+dC*dt/dx/dy;
    
    %set BCs
    %u
    u(cell_u==5|cell_u==8|cell_u==10) = u_in_left; %inlet and corners
    u(cell_u==3|cell_u==4|cell_u==12|cell_u==13) = 0; %left/right walls and conrners
    u(end,1:round(w_left/dy))=u(end-1,1:round(w_left/dy)); %outlet
    
    %check outlet
    m_out=sum(u(end,1:round(w_left/dy))*dy); %mass outlet
    if m_out>1e-8
        u(end,1:round(w_left/dy))=u(end,1:round(w_left/dy))*m_in/m_out; %set m_out to m_in (no mass source)
    else
        u(end,1:round(w_left/dy))=u(end,1:round(w_left/dy))+m_in/w_left; %if m_out=0
    end
    
    %v
    v(cell_v==6|cell_v==12|cell_v==13) = v_in_up; %inlet and corners
    v(cell_v==1|cell_v==2|cell_v==8|cell_v==9|cell_v==10|cell_v==11) = 0; %lower/upper walls and corners
    
    %% pressure correction
    %set up vector of b (based on delta_m_dot(i,j))
    for i=1:m
        for j=1:n
            if ~isnan(cell_P(i,j))
                b(cell_number(i,j)) = ( (u(i+1,j)-u(i,j))/dx + (v(i,j+1)-v(i,j))/dy ) *rho/dt;
            end
        end
    end
    
    %set p_prime(corner 4) to 0 (is known=P_0) (doesn't need correction)
    %doesn't need anything to do except reducing matrix size
    
    %size of A already reduced
    P_prime_vector=A\b(1:end-1);
    P_prime_vector=[P_prime_vector;0];
    for i=1:m
        for j=1:n
            if cell_number(i,j)~=0
                P_prime(i,j) = P_prime_vector(cell_number(i,j));
            end
        end
    end
    
    P=P+P_prime;
    
    %% velocity correction
    %u & v since boundary elements don't need correction can use same loop (m*n)
    for i=1:m
        for j=1:n
            %for middle, low/upper walls, and upper inlet
            if cell_u(i,j)==0 || cell_u(i,j)==1 || cell_u(i,j)==2 || cell_u(i,j)==6
                u_prime(i,j) = -dt/rho/dx*( P_prime(i,j)-P_prime(i-1,j) );
            end
            %for middle, left/right walls, left inlet, and outlet
            if cell_v(i,j)==0 || cell_v(i,j)==3 || cell_v(i,j)==4 || cell_v(i,j)==5 || cell_v(i,j)==7
                v_prime(i,j) = -dt/rho/dy*( P_prime(i,j)-P_prime(i,j-1) );
            end
        end
    end
    
    u=u+u_prime;
    v=v+v_prime;
    
    %% save parameters
    if mod(round(t,5),0.1)==0
        t
        if mod(round(t,5),1)==0
            save_time=t
            t_save(end+1)=t;
            u_save(:,:,end+1)=u;
            v_save(:,:,end+1)=v;
            C_save(:,:,end+1)=C;
            P_save(:,:,end+1)=P;
            save('Data');
            
            %% check steady state
            if max(max(abs(u_save(:,:,end)-u_save(:,:,end-1))))<epsilon_u
                if max(max(abs(v_save(:,:,end)-v_save(:,:,end-1))))<epsilon_v
                    if max(max(abs(P_save(:,:,end)-P_save(:,:,end-1))))<epsilon_P
                        if max(max(abs(C_save(:,:,end)-C_save(:,:,end-1))))<epsilon_C
                            t_end=t
                            break
                        end
                    end
                end
            end
        end
    end
    
    %% check Re and Courant #
    %Re
    Re_u=max(max(abs(u)))*w_left/nu;
    Re_v=max(max(abs(v)))*w_up/nu;
    %check laminar
    if max(Re_u,Re_v)>1000
        'Re # too high'
        break
    end
    %CFL
    CFL_u=dt/min(min(dx./abs(u)));
    CFL_v=dt/min(min(dy./abs(v)));
    %check stability
    if max(CFL_u,CFL_v)>0.1
        'dt too high (CFL)'
        break
    end
    
end