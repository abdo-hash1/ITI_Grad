function [I,Q] = QPSK(data)

maped_I= [-1,-1,1,1]./sqrt(2);
maped_Q= [-1,1,-1,1]./sqrt(2);


data = reshape(data,2,[])' ;

for i=1:size(data,1)
    I(i) = maped_I(bit2int(data(i,:)',2)+1);
    Q(i) = maped_Q(bit2int(data(i,:)',2)+1);
end

end

