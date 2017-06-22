clearvars
x0=0;
xn=10;
kvZahl=25;
RB1Wert=50;
RB2Wert=100;
diffKoeff=4;
konvKoeff=2;
stromrichtung=true;
if konvKoeff<0
    stromrichtung=false; %UDS
end
dx=(xn-x0)/(kvZahl);
x= x0+dx/2:dx:xn-dx/2;
syms k T_i T_ip1 T_im1 T_wp T_ep T_w T_e EQ
f(k) = 0*k/k;
for i=1:1:length(x)
    if i==1 && i==length(x)
        T_w = RB1Wert;
        T_e = RB2Wert;
        if ~stromrichtung 
            T_w = RB2Wert;
            T_e = RB1Wert;
        end
        T_wp = (T_i - RB1Wert)/(x(i)-x0);
        T_ep = (RB2Wert-T_i)/(xn-x(i));
        F = 0.5*(f(x0)+f(xn))*(xn-x0);
    elseif i==1
        T_w = RB1Wert;
        T_e = T_i;
        if ~stromrichtung 
            T_w = RB1Wert;
            T_e = T_ip1;
        end
        T_wp = (T_i - RB1Wert)/(x(i)-x0);
        T_ep = (T_ip1 - T_i)/(x(i+1)-x(i));
        F = 0.5*(f(x0)+f(x(i)+dx/2))*(x(i)+dx/2-x0);
    elseif i==length(x)
        T_w = T_im1;
        T_e = RB2Wert;
        if ~stromrichtung 
            T_w = T_i;
            T_e = RB2Wert;
        end
        T_ep = (RB2Wert-T_i)/(xn-x(i));
        T_wp = (T_i - T_im1)/(x(i)-x(i-1));
        F = 0.5*(f(x(i)-dx/2)+f(xn))*(xn-(x(i)-dx/2));
    else
        T_w = T_im1;    %UDS O(x²)
        T_e = T_i;
        if ~stromrichtung 
            T_w = T_i;
            T_e = T_ip1;
        end
        T_ep = (T_ip1 - T_i)/(x(i+1)-x(i)); %ZDS
        T_wp = (T_i - T_im1)/(x(i)-x(i-1));
        F = 0.5*(f(x(i)-dx/2)+f(x(i)+dx/2))*dx;
    end
    t = konvKoeff*T_e - konvKoeff*T_w + diffKoeff*T_ep - diffKoeff*T_wp == -F;
    EQ(i)=t;    
end
[A, b] = equationsToMatrix(EQ, [T_im1 T_i T_ip1]); %Gesamtmatrix
C = zeros(kvZahl,kvZahl-3);
A = [A C];

cnt=1;      %umformvektor
for j=1:kvZahl
    K(j)=cnt;   
    cnt=cnt-1;
end
[m, n] = size(A);    %umformer
ind = bsxfun(@plus,m*mod(bsxfun(@plus,1:n,K(:)-1),n),(1:m)');
B = zeros(size(A));
B(:) = A(ind);
if kvZahl == 1  %schneider
    B(:,2:end)=[];
elseif kvZahl == 2
    B(:,3)=[];
end
T = linsolve(B,b);  %löser
double(T)