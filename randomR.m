function outIndex = randomR(weight)
L=length(weight);
 
outIndex=zeros(1,L);

 
u=unifrnd(0,1,1,L);
u=sort(u);

cdf=cumsum(weight);

i=1;
for j=1:L
  
    while (i<=L) && (u(i)<=cdf(j))
 
        outIndex(i)=j;
 
        i=i+1;
    end
end