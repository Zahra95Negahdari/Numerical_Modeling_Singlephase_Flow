clear
clc
close all
format short
k=.1;
phi=.2;
mo=5;
B0=1.2;
ip=5000;
cl=1e-6;
dx=1000;
dy=dx;
dz=50;
dt=0.5;
rw=.5;
s1=2;
s2=-2;
q=-100;
alpha=5.615;
Bc=1.127;
B=B0;
T=(Bc*k*dy*dz)/(mo*B*dx);
T1=T;
T2=T;
A=(dx*dy*dz*phi*cl)/(alpha*B*dt);
p=zeros(1,12);
p(1,:)=ip;
pwf2(1)=ip;
pwf4(1)=ip;
pwf8(1)=ip;
pwf22(1)=ip;
pwf44(1)=ip;
pwf88(1)=ip;
q2(1)=0;
q4(1)=0;
q8(1)=0;
req=.2*dx;
eps=1;
for n=2:361
    for i=2:10
        B=B0/(1+cl*(p(n-1,i)-ip));
        B1=B0/(1+cl*(p(n-1,i-1)-ip));
        B2=B0/(1+cl*(p(n-1,i+1)-ip));
        B11=(B1+B)/2;
        B22=(B2+B)/2;
        T1=(Bc*k*dy*dz)/(mo*B11*dx);
        T2=(Bc*k*dy*dz)/(mo*B22*dx);
        j=(2*pi*Bc*k*dz)/(mo*B*(log(req/rw)+s1));
        switch i
            case 2
                m(i-1,:)=fully(i,0,T2,p(n-1,i),p(n-1,i),p(n-1,i-1),p(n-1,i+1),A,eps,0,0,0);
            case 3
                m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n-1,i),p(n-1,i-1),p(n-1,i+1),A,eps,-100,0,0);
            case 5
                m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n-1,i),p(n-1,i-1),p(n-1,i+1),A,eps,0,j,2200);
            case 9
                m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n-1,i),p(n-1,i-1),p(n-1,i+1),A,eps,-100,0,0);
            case 10
                m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n-1,i),p(n-1,i-1),p(n-1,i+1),A,eps,0,0,0);
            otherwise
                m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n-1,i),p(n-1,i-1),p(n-1,i+1),A,eps,0,0,0);
        end
    end
    x=m(:,2:10);
    y=m(:,12);
    dp=x\y;
    p(n,2:10)=p(n-1,2:10)-(dp');
    p(:,11)=ip;
    p(n,1)=p(n,2);
    for jj=1:9
        while abs(dp(jj))>1
            for i=2:10
                B=B0/(1+cl*(p(n,i)-ip));
                B1=B0/(1+cl*(p(n,i-1)-ip));
                B2=B0/(1+cl*(p(n,i+1)-ip));
                B11=(B1+B)/2;
                B22=(B2+B)/2;
                T1=(Bc*k*dy*dz)/(mo*B11*dx);
                T2=(Bc*k*dy*dz)/(mo*B22*dx);
                j=(2*pi*Bc*k*dz)/(mo*B*(log(req/rw)+s1));
                switch i
                    case 2
                        m(i-1,:)=fully(i,0,T2,p(n-1,i),p(n,i),p(n,i-1),p(n,i+1),A,eps,0,0,0);
                    case 3
                        m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n,i),p(n,i-1),p(n,i+1),A,eps,-100,0,0);
                    case 5
                        m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n,i),p(n,i-1),p(n,i+1),A,eps,0,j,2200);
                    case 9
                        m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n,i),p(n,i-1),p(n,i+1),A,eps,-100,0,0);
                    case 10
                        m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n,i),p(n,i-1),p(n-1,i+1),A,eps,0,0,0);
                    otherwise
                        m(i-1,:)=fully(i,T1,T2,p(n-1,i),p(n,i),p(n,i-1),p(n,i+1),A,eps,0,0,0);
                end
            end
            x=m(:,2:10);
            y=m(:,12);
            dp=x\y;
            p(n,2:10)=p(n,2:10)-(dp)';
            p(:,11)=ip;
            p(n,1)=p(n,2);
        end
    end
end
for n=2:361
    acc=0;
    for i=2:11
        B=B0/(1+cl*(p(n-1,i)-ip));
        B1=B0/(1+cl*(p(n-1,i-1)-ip));
        B2=B0/(1+cl*(p(n-1,i+1)-ip));
        B11=(B1+B)/2;
        B22=(B2+B)/2;
        T1=(Bc*k*dy*dz)/(mo*B11*dx);
        T2=(Bc*k*dy*dz)/(mo*B22*dx);
        j1=(2*pi*Bc*k*dz)/(mo*B*(log(req/rw)+s1));
        j2=(2*pi*Bc*k*dz)/(mo*B*(log(req/rw)+s2));
        switch i
            case 3
                flow(n)=T1*(p(n,2)-p(n,3));
                pwf2(n)=p(n,i)-q/(-j1);
                pwf22(n)=p(n,i)-q/(-j2);
                q2(n)=100;
                delta_p1(i)=(dx*dy*dz*cl*phi/(dt*alpha))*(p(n-1,i)-p(n,i));
                delta_p2(i)=p(n-1,i)-p(n,i);
            case 5
                q4(n)=j1*(p(n,i)-2200);
                q44(n)=j2*(p(n,i)-2200);
                pwf4(n)=2200;
                pwf44(n)=2200;
                delta_p1(i)=(dx*dy*dz*cl*phi/(dt*alpha))*(p(n-1,i)-p(n,i));
                delta_p2(i)=p(n-1,i)-p(n,i);
            case 9
                pwf8(n)=p(n,i)-q/(-j1);
                pwf88(n)=p(n,i)-q/(-j2);
                q8(n)=100;
                delta_p1(i)=(dx*dy*dz*cl*phi/(dt*alpha))*(p(n-1,i)-p(n,i));
                delta_p2(i)=p(n-1,i)-p(n,i);
            case 11
                delta_p1(i)=(dx*dy*dz*cl*phi/(dt*alpha))*(p(n-1,i)-p(n,i));
                in(n)=T1*abs(p(n,i)-p(n,i-1))+sum(delta_p1);
                entering=T1*abs(p(n,i)-p(n,i-1));
                delta_p2(i)=p(n-1,i)-p(n,i);
            otherwise
                delta_p1(i)=(dx*dy*dz*cl*phi/(dt*alpha))*(p(n-1,i)-p(n,i));
                delta_p2(i)=p(n-1,i)-p(n,i);
        end
        acc=(dx*dy*dz/alpha)*(phi*cl/B0)*delta_p2(i)+acc;
    end
    mb(n-1)=(abs(q4(n)+200-in(n))/(q4(n)+200))*100;
    IMB(n-1)=acc/(dt*(200+q4(n)-entering));
end
for n=1:361
    figure(1)
    plot(p(n,2:11))
    hold on
    xlabel('block number')
    ylabel('pressure')
end
p(91,2:11)
p(181,2:11)
p(271,2:11)
p(361,2:11)
figure(2)
plot(flow)
grid on
xlabel('time steps')
ylabel('flow between grid block 1&2')
%=============================
figure(3)
plot(q4)
hold on
plot(q2)
grid on
xlabel('time steps')
ylabel('well flow rate VS time ')
legend('well flow for Grid block4','well flow for Grid blocks 2&8')
%================================
figure(4)
plot(pwf2)
hold on
plot(pwf4)
hold on
plot(pwf8)
grid on
xlabel('time steps')
ylabel('BHP for Grid blocks')
legend('BHP2','BHP4','BHP8')
%=======================
figure(5)
plot(pwf2)
hold on
plot(pwf22)
grid on
xlabel('time steps')
ylabel('BHP for Grid block2')
legend('s=2','s=-2')
%=====================
figure(6)
plot(pwf4)
hold on
plot(pwf44)
grid on
xlabel('time steps')
ylabel('BHP for Grid block4')
legend('s=2','s=-2')
%===================
figure(7)
plot(pwf8)
hold on
plot(pwf88)
grid on
xlabel('time steps')
ylabel('BHP for Grid block8')
legend('s=2','s=-2')
%======================
figure(8)
plot(q4)
hold on
plot(q44)
grid on
xlabel('time steps')
ylabel('flow rate for grid block4')
legend('s=2','s=-2')
%======================
figure(9)
plot(mb)
grid on
title('MB with eq-21')
%========================
figure(10)
plot(IMB)
grid on
title('MB with eq-23')