function[x] = initialize_curtailment(xmin,xmax,J,T,S,it,CM)
  x=zeros(J,T,S,it);
for s=1:S
for j=1:J
        M=0;
        while M~=1
                 for t=1:T
                 x(j,t,s,1) = xmin(j) + (xmax(j)-xmin(j))*rand(1,1);
                 end
    
                 cond2=0;

                 for t=1:T
                 cond2=cond2+x(j,t,s,1);
                 end

                 if (cond2)<=CM(j)
                 M=M+1;
                 
                 end    
         end
end
end
end