function[ff,x1]=objective_function(x,y,a,b,ld,ww,T,J,I,s,ii,pdg)
% optimization function
x1=zeros(1,1);
x3=zeros(1,1);
for t=1:T
    for i=1:I
        x1=x1+(a(i)*((pdg(i,t)*1000).^2)+(b(i)*pdg(i,t)*1000));
    end
end


ff(s,1)= ww*x1;
end
