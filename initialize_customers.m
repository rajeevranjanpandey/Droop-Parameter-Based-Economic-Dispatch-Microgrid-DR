function[y] = initialize_customers(J,T,S,K1,K2,x,th,UB)
y=zeros(J,T,S,1);
for s=1:S
    %
    M=0;
    while M~=1
        %
        for j=1:J
            for t=1:T
                y(j,t,s,1) = K1(j)*x(j,t,s,1).^2+K2(j)*x(j,t,s,1)*(1-th(j));
            end
        end
        
        
        cond1=0;
        for t=1:T
            for j=1:J
                cond1=cond1+y(j,t,s,1);
            end
        end
        
        if (cond1<=UB)
            M=M+1;
            
        end
    end
end
end

