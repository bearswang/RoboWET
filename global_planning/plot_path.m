load('example_-60.mat');

%% generate the map
mon=1;
hold on;
% plot starting point
p0=plot(UGV_location_all(1,1,mon),UGV_location_all(2,1,mon),'ob','MarkerSize',15);
% plot all vertices
p1=plot(UGV_location_all(1,:,mon),UGV_location_all(2,:,mon),'square','MarkerSize',8);
% plot IoT users
p2=plot(user_location_all(1,:,mon),user_location_all(2,:,mon),'+','MarkerSize',8);

D=D_all(:,:,mon);
v_opt= v_all(:,3,mon);

%% plot the shortest path visiting all vertices
% solve the TSP using Gurobi
cvx_begin quiet
cvx_solver Mosek
% original variables
variable x(M,M) binary
variable u(M)

minimize trace(transpose(D)*x) % total path length
subject to

% divergence constraint
for m=1:M
    sum(x(m,:))==1;
    sum(x(:,m))==1;
    x(m,m)==0;
end

% subtour elimination constraint
for m=2:M
    for j=2:M
        if m~=j
            u(m)-u(j)+(M-1)*x(m,j)+(M-3)*x(j,m)<=M-2;
        end
    end
end

for m=2:M
    u(m)>=1;
    u(m)<=M-1;
end

cvx_end
% 
% for m=1:M
%     for j=1:M
%         if x(m,j)==1
%             p3=plot(UGV_location_all(1,[m,j],mon),UGV_location_all(2,[m,j],mon),'-k','LineWidth',2);
%         end
%     end
% end

%% plot the proposed path at -60 dBm
MM=sum(v_opt);
if MM==1
    Sm=0;
    Em=0;
    x=zeros(M,M);
    p4=plot(UGV_location_all(1,1,mon),UGV_location_all(2,1,mon),'-k','LineWidth',2);
else
    epsilon=1e3;
    % solve the TSP using Gurobi
    cvx_begin quiet
    cvx_solver Mosek
    % original variables
    variable x(M,M) binary
    variable u(M)
    
    minimize trace(transpose(D)*x) % total path length
    subject to
    
    % divergence constraint
    for m=1:M
        if v_opt(m)==1;
            sum(x(m,:))==1;
            sum(x(:,m))==1;
        else
            sum(x(m,:))==0;
            sum(x(:,m))==0;
        end
        x(m,m)==0;
    end

    % subtour elimination constraint
    for m=2:M
        for j=2:M
            if m~=j && v_opt(m)==1 && v_opt(j)==1
                u(m)-u(j)+(MM-1)*x(m,j)+(MM-3)*x(j,m)<=MM-2;
            end
        end
    end
    
    for m=2:M
        u(m)>=v_opt(m);
        u(m)<=(MM-1)*v_opt(m);
    end
    
    cvx_end
    
    for m=1:M
        for j=1:M
            if x(m,j)>0.5
                p4=plot(UGV_location_all(1,[m,j],mon),UGV_location_all(2,[m,j],mon),'-k','LineWidth',2);
            end
        end
    end
    
end

load('example_-90.mat');
%% plot the proposed path at -80 dBm
MM=sum(v_opt);
if MM==1
    Sm=0;
    Em=0;
    x=zeros(M,M);
    p5=plot(UGV_location_all(1,1,mon),UGV_location_all(2,1,mon),'--r','LineWidth',2);
else
    epsilon=1e3;
    % solve the TSP using Gurobi
    cvx_begin quiet
    cvx_solver Mosek
    % original variables
    variable x(M,M) binary
    variable u(M)
    
    minimize trace(transpose(D)*x) % total path length
    subject to
    
    % divergence constraint
    for m=1:M
        if v_opt(m)==1;
            sum(x(m,:))==1;
            sum(x(:,m))==1;
        else
            sum(x(m,:))==0;
            sum(x(:,m))==0;
        end
        x(m,m)==0;
    end

    % subtour elimination constraint
    for m=2:M
        for j=2:M
            if m~=j && v_opt(m)==1 && v_opt(j)==1
                u(m)-u(j)+(MM-1)*x(m,j)+(MM-3)*x(j,m)<=MM-2;
            end
        end
    end
    
    for m=2:M
        u(m)>=v_opt(m);
        u(m)<=(MM-1)*v_opt(m);
    end
    
    cvx_end
    
    for m=1:M
        for j=1:M
            if x(m,j)>0.5
                p5=plot(UGV_location_all(1,[m,j],mon),UGV_location_all(2,[m,j],mon),'--r','LineWidth',2);
            end
        end
    end
    
end

legend([p0,p1,p2,p4,p5],'Starting point','Vertices','IoT users','Optimal path1','Optimal path2');

