function childern = branch( parent )
childern=cell(2,1);
coin=randi([0,1],1); % randonly select a child node
childern{1}=[parent,coin];
childern{2}=[parent,1-coin];
end

