function[mp,nq] = initialize_droop(mpmin,mpmax,nqmin,nqmax,T,I,S,it)

mp=zeros(I,T,S,it);
nq=zeros(I,T,S,it);
% power generation by generators
for s=1:S
    for i=1:I
            for t=1:T
                mp(i,t,s,1) = mpmin + (mpmax-mpmin)*rand(1,1);
                nq(i,t,s,1) = nqmin + (nqmax-nqmin)*rand(1,1);
            end

        end
    end
end