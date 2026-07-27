function[pdg0,qdg0,w,Pl,Ql,V,del,PLoss,QLoss]=droop_load_flow(ii,x,s,mp,nq,t,ndg,St,Pl00,Ql00)
nbuss=38;                    % no. of bus
nbus=nbuss-(5-ndg);
BMva=1000;                  % base MVA = 1000 kw
KVb=12.66;                  % base kv = 12.66 kv
Zb=((KVb^2)*1000)/BMva;     % base impedence
%Y = ybusppg(nbus,Zb);       % Calling ybusppg.m to get Y-Bus Matrix..
busdd=bus_data(nbuss);      % Calling busdatas..
busd = busdd(1:nbus,:);
% bus = busd(:,1);            % Bus Number..
type = busd(:,2);           % Type of Bus 1-Slack, 2-PV, 3-PQ..
V = busd(:,3);              % Specified Voltage..
del = busd(:,4);            % Voltage Angle..
%Vm = pol2rect(V,del);       % Converting polar to rectangular..
Pg = busd(:,5)/BMva;        % PGi..
Qg = busd(:,6)/BMva;        % QGi..
% Pl0 = busd(:,7)/BMva;        % PLi..
% Ql0 = busd(:,8)/BMva;        % QLi..
Pl0=Pl00(:,t);
Ql0=Ql00(:,t);
Pl0(24)=Pl00(24,t)-St(t)/BMva;

Qmin = busd(:,9)/BMva;      % Minimum Reactive Power Limit..
Qmax = busd(:,10)/BMva;     % Maximum Reactive Power Limit..
% P = Pg - Pl;                % Pi = PGi - PLi..
% Q = Qg - Ql;                % Qi = QGi - QLi..
% Psp = P;                    % P Specified..
% Qsp = Q;                    % Q Specified..
table1= [  0.92   4.04
    0.18   6.00
    1.51   3.40
    0      0   ];
alpha=table1(:,1);
beta=table1(:,2);
cus=busd(:,11);
kpf=0;
kqf=-2;
V0=1.01;
w0=1;

%    DG No. Location mp nq w0 v0 Qmax_dg

dg_data=[ 1 34 5.102*10^(-3) .02      1   1.01   0.9
    2 35 1.502*10^(-3) .03333   1   1.01   0.9
    3 36 4.506*10^(-3) .02      1   1.01   0.9
    4 37 2.253*10^(-3) .05      1   1.01   0.9
    5 38 2.253*10^(-3) .05      1   1.01   0.3];
dg=dg_data(:,1);
loc=dg_data(:,2);
mppp=mp(:,t,s,ii+1);
nqqq=nq(:,t,s,ii+1);
% w0=dg_data(:,5);
% v0=dg_data(:,6);
Qmax_dg=dg_data(:,7);
pv = find(type == 2 | type == 1);   % PV Buses..
pq = find(type == 3);               % PQ Buses..
PQ=sort(find(type==3|type==6));
NPQ=length(PQ);
vf=find(type==6);
npv = length(pv);                   % No. of PV buses..
npq = length(pq);                   % No. of PQ buses.
nvf=length(vf);
Tol = 1;
Iter = 1;
w=1;

% --- safety cap so a bad candidate droop combination (e.g. from TFO's
%     random exploration far outside the PSO neighborhood) cannot make
%     this loop hang forever. Same NR body, just bounded. ---
maxIter = 100;

while (Tol > 1e-5) && (Iter <= maxIter)
    
    
    Vm = polar_to_rect(V,del);       % Converting polar to rectangular..
    Y = build_ybus(nbus,Zb,w,nbuss);       % Calling ybusppg.m to get Y-Bus Matrix..
    theta=angle(Y);
    Ym=abs(Y);
    G = real(Y);                % Conductance matrix..
    B = imag(Y);                % Susceptance matrix..
    Zk=1./Y;
    Rk=real(Zk);
    Xk=imag(Zk);
    % calculation of P and Q
    %  load model
    Pl=zeros(nbus,1);
    Ql=zeros(nbus,1);
    for i=1:nbus
        Pl(i)=Pl0(i)*((V(i)/V0)^alpha(cus(i)))*(1+kpf*(w-w0));
        Ql(i)=Ql0(i)*((V(i)/V0)^beta(cus(i)))*(1+kqf*(w-w0));
    end
    
    % P=zeros(nbus,1);
    % Q=zeros(nbus,1);
    % for VF bus
    n=1;
    for i=1:nbus
        if type(i)==6
            Pg(i)=(1/mppp(n))*(w0-w);
            Qg(i)= (1/nqqq(n))*(V0-V(i));
            % Dg reactive power limit
                        if Qg(i)>Qmax_dg(n)
                            Qg(i)=Qmax_dg(n);
                        end
            n=n+1;
        end
    end
    PP = Pg - Pl;                % Pi = PGi - PLi..
    QQ = Qg - Ql;                % Qi = QGi - QLi..
    Psp = PP;                    % P Specified..
    Qsp = QQ;                    % Q Specified..
    
    % Calculate P and Q for other bus
    P=zeros(nbus,1);
    Q=zeros(nbus,1);
    
    for i = 1:nbus
        for k = 1:nbus
            P(i) = P(i) + V(i)* V(k)*(G(i,k)*cos(del(i)-del(k)) + B(i,k)*sin(del(i)-del(k)));
            Q(i) = Q(i) + V(i)* V(k)*(G(i,k)*sin(del(i)-del(k)) - B(i,k)*cos(del(i)-del(k)));
        end
    end
    
    
    % loss formula
    
    Ploss=0;
    Qloss=0;
    for i=1:nbus
        for j=1:nbus
            Ploss=Ploss+real (Y(i,j)* (conj(Vm(i))*Vm(j) + conj(Vm(j))* Vm(i)));
            Qloss=Qloss+imag (Y(i,j)* (conj(Vm(i)) *Vm(j) + conj(Vm(j))* Vm(i)));
        end
    end
    PLoss=0.5*Ploss;
    QLoss=-0.5*Qloss;
    
    
    Ptot=sum(Pl)+PLoss;
    Qtot=sum(Ql)+QLoss;
     Psys=sum(Pg);
    Qsys=sum(Qg);
%     for i=1:length(dg)
%         Psys=Psys+((1/mp(i))*(w0-w));
%         Qsys=Qsys+((1/nq(i))*(V0-V(loc(i))));
%     end
    dPa=Psp-P;
    dQa=Qsp-Q;
    dP=dPa(2:nbus);
    dQ=dQa(2:nbus);
    dx=Ptot-Psys;
    dy=Qtot-Qsys;
    M=[dP;dQ;dx;dy];
    % Jacobian
    % J1 - Derivative of Real Power Injections with Angles..
    J1 = zeros(nbus-1,nbus-1);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J1(i,k) = J1(i,k) + V(m)* V(n)*(-G(m,n)*sin(del(m)-del(n)) + B(m,n)*cos(del(m)-del(n)));
                end
                J1(i,k) = J1(i,k) - V(m)^2*B(m,m);
            else
                J1(i,k) = V(m)* V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
            end
        end
    end
    
    % J2 - Derivative of Real Power Injections with V..
    J2 = zeros(nbus-1,NPQ);
    for i = 1:(nbus-1)
        m = i+1;
        for k = 1:NPQ
            n = PQ(k);
            if n == m
                for n = 1:nbus
                    J2(i,k) = J2(i,k) + V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J2(i,k) = J2(i,k) + V(m)*G(m,m);
            else
                J2(i,k) = V(m)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J3 - Derivative of Reactive Power Injections with Angles..
    J3 = zeros(NPQ,nbus-1);
    for i = 1:NPQ
        m = PQ(i);
        for k = 1:(nbus-1)
            n = k+1;
            if n == m
                for n = 1:nbus
                    J3(i,k) = J3(i,k) + V(m)* V(n)*(G(m,n)*cos(del(m)-del(n)) + B(m,n)*sin(del(m)-del(n)));
                end
                J3(i,k) = J3(i,k) - V(m)^2*G(m,m);
            else
                J3(i,k) = V(m)* V(n)*(-G(m,n)*cos(del(m)-del(n)) - B(m,n)*sin(del(m)-del(n)));
            end
        end
    end
    
    % J4 - Derivative of Reactive Power Injections with V..
    
    J4 = zeros(NPQ,NPQ);
    for i = 1:NPQ
        m = PQ(i);
        for k =1: NPQ
            n = PQ(k);
            if n == m
                for n = 1:nbus
                    J4(i,k) = J4(i,k) + V(n)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                    %                     if type(m)==6
                    %                         J4(i,k)=(-1)/nq(kk);
                    %                         kk=kk+1;
                    %                     end
                end
                J4(i,k) = J4(i,k) - V(m)*B(m,m);
            else
                J4(i,k) = V(m)*(G(m,n)*sin(del(m)-del(n)) - B(m,n)*cos(del(m)-del(n)));
                %                 if type(m)==6
                %                         J4(i,k)=0;
                %                     end
            end
        end
    end
    ll=1;
    for kk=33:32+ndg
        J4(kk,kk)=J4(kk,kk)+(1)/nqqq(ll);
        ll=ll+1;
    end
    %
    J13=zeros(NPQ,1);
    for i=1:NPQ
        m=PQ(i);
        for n=1:nbus
            J13(i)=J13(i)+ ( (-1)*  (   (Xk(m,n)^2/w)   /   (Rk(m,n)^2+Xk(m,n)^2 )^1.5    )     * V(m)* V(n)* cos(del(m)-del(n)-theta(m,n))   +  (-1)*  ((  Xk(m,n)/(w*Rk(m,n)))   /   (1+(Xk(m,n)/Rk(m,n))^2 ) ) *   V(m)*V(n)* Ym(m,n)*sin(del(m)-del(n)-theta(m,n))  );
        end
    end
    
    J23=zeros(NPQ,1);
    for i=1:NPQ
        m=PQ(i);
        for n=1:nbus
            J23(i)=J23(i)+  ( (-1)*  (   (Xk(m,n)^2/w)   /   (Rk(m,n)^2+Xk(m,n)^2 )^1.5    )     * V(m)* V(n)* sin(del(m)-del(n)-theta(m,n))   -  (-1)*  ((  Xk(m,n)/(w*Rk(m,n)))   /   (1+(Xk(m,n)/Rk(m,n))^2 ) ) *   V(m)*V(n)* Ym(m,n)*cos(del(m)-del(n)-theta(m,n))  );
        end
    end
    
    J14=zeros(NPQ,1);
    for i=1:NPQ
        m=PQ(i);
        J14(i)=J14(i)+ V(m)*Ym(m,1)  *cos(del(m)-del(1)-theta(m,1));
    end
    
    J24=zeros(NPQ,1);
    for i=1:NPQ
        m=PQ(i);
        J24(i)=J24(i)+ V(m)*Ym(m,1)  * sin(del(m)-del(1)-theta(m,1));
    end
    J31=zeros(1,NPQ);
    J32=zeros(1,NPQ);
    J41=zeros(1,NPQ);
    
    J42=zeros(1,NPQ);
    n=1;
    for i=1:NPQ
        m=PQ(i);
        
        if type(m)==6
            J42(1,i)=(-1)/nqqq(n);
            n=n+1;
            
            %       else J42(1,i)=0;
        end
        
    end
    
    J34=zeros(1,1);
    J43=zeros(1,1);
    
    J33=zeros(1,1);
    
    for i=1:ndg
        J33=J33+((-1)/mppp(i));
    end
    
    J44=zeros(1,1);
    n=1;
    for i=1:NPQ
        if type(1)==6
            J44=(-1)/nqqq(n);
            n=n+1;
            
        end
    end
    
    J= [J1 J2 J13 J14 ; J3 J4 J23 J24; J31 J32 J33 J34; J41 J42 J43 J44];

    if rcond(J) < 1e-14 || any(~isfinite(J(:)))
        % ill-conditioned / singular Jacobian for this droop combination —
        % stop iterating instead of hanging or throwing from inv(J).
        break;
    end

    X=inv(J)*M;
    ddel=X(1:NPQ);
    dV=X((NPQ+1):(2*nbus-2));
    dw=X(2*nbus-1);
    dV1= X(2*nbus);
    w=w+dw;
    V(1)=V(1)+dV1;
    for i=1:NPQ
        m=PQ(i);
        V(m)=V(m)+dV(i);
        del(m)=del(m)+ddel(i);
    end
    Iter = Iter + 1;
    Tol = max(abs(M));                  % Tolerance..
    
end

pdg0=Pg(34:(33+ndg));
qdg0=Qg(34:(33+ndg));

end
