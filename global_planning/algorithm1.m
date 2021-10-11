function [E_opt, v_opt, history, Em_opt, Et_opt]= algorithm1(K, M, h, g, D, T, gamma, N0, miu, E_current, iter_ls)
beta=0.5; % performance loss due to modulation
eta=1/1.1; % backscatter efficiency

A=zeros(K,M);
for k=1:K
    for m=1:M
        A(k,m)=beta*eta*norm(h(k,m))^2*norm(g(k,m))^2/N0; % equation (7)
    end
end

Pool=cell(2,1); % initial living pool
Pool{1}=[1,0];
Pool{2}=[1,1];

iter=iter_ls; % add the number of iterations of Algorithm 2
history.size(1:iter)=2^(M-1);

while ~isempty(Pool)
    % compute the number of candidate solutions in the pool
    number=0;
    for i=1:length(Pool)
        number=number+2^(M-length(Pool{i}));
    end
    history.size(iter)=number;
    
    % select a node
    node = Pool{1};
    [Em, Et, E_total] = bounding(K, M, h, g, D, T, gamma, N0, miu, node); % compute the lower bound using (22)
    if length(node)==M
        if E_total<=E_current
            E_current=E_total;
            v_current=[node];
            Em_current=Em;
            Et_current=Et;
            Pool(1)=[];
            continue;
        end
    end
    lower_bound = E_total;
    if lower_bound > E_current % discard
        Pool(1)=[];
    else % branch
        Pool(1)=[];
        if length(node)<M
            children = branch(node);
            Pool=[children{2}, Pool];
            Pool=[children{1}, Pool];
        end
    end

    if mod(iter,100)==0
        disp('100 iterations!');
    end
    iter=iter+1; % increase the iteration counter
    
end
disp('Total number of iterations:');
disp(iter);
E_opt=E_current;
v_opt=v_current;
Em_opt=Em_current;
Et_opt=Et_current;

end

