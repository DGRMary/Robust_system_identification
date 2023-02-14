clear 
close all 
clc 

load data_exam_2A_hammer_ver2.mat

u = Input; 
y = Output;

Ts

Dh = 0.02;
N = length(u);
n = 1;

PUI = zeros(4,2);

for c_PUI = 1:2
    if c_PUI == 1
        min_th = 1;
    else 
        min_th= -1;
    end
    
    for r_PUI = 1:4
        %obj function 
        objPoly.typeCone = 1;
        objPoly.dimVar = 4+2*N;
        objPoly.degree = 1;
        objPoly.noTerms = 1;

        support_mtx = zeros(1,4+2*N);
        support_mtx(r_PUI) = 1;
        coef_vect = min_th;

        objPoly.supports = sparse(support_mtx);
        objPoly.coef = sparse(coef_vect);

        ineqPolySys = cell(1,2*N-1);

        for t = 2:N
            ineqPolySys{t-n}.typeCone = -1;
            ineqPolySys{t-n}.dimVar = 4+2*N;
            ineqPolySys{t-n}.degree = 2;
            ineqPolySys{t-n}.noTerms = 6;

            support_mtx = zeros(6,4+2*N);
            support_mtx(2,4+t) = 1;
            support_mtx(3,1) = 1;
            support_mtx(4,1) = 1;support_mtx(4,4+t-1) = 1;
            support_mtx(5,2) = 1;support_mtx(5,4+N+t) = 1;
            support_mtx(6,3) = 1;support_mtx(6,4+N+t-1) = 1;

            coef_vect = [y(t) -1 y(t-1) -1 -1 -1]';

            ineqPolySys{t-n}.supports = sparse(support_mtx);
            ineqPolySys{t-n}.coef = sparse(coef_vect);
        end

        for t = 1:N
            ineqPolySys{N-1+t}.typeCone = -1;
            ineqPolySys{N-1+t}.dimVar = 4+2*N;
            ineqPolySys{N-1+t}.degree = 1;
            ineqPolySys{N-1+t}.noTerms = 3;

            support_mtx = zeros(3,4+2*N);
            support_mtx(1,4+N+t) = 1;
            support_mtx(3,4) = 1;
           
            coef_vect = [1 -0.25*u(t) -(u(t))^3]';

            ineqPolySys{N-1+t}.supports = sparse(support_mtx);
            ineqPolySys{N-1+t}.coef = sparse(coef_vect);
        end

            ubd = [1e10*ones(1,4) Dh*ones(1,N) 1e10*ones(1,N)];
            ubd(1) = 0.999;
            lbd = -ubd;

            param.POPsolver = 'active-set';
            param.relaxOrder = 2;

            [param,SDPobjValue,POP,cpuTime,SDPsolverInfo,SDPinfo] =...
                sparsePOP(objPoly,ineqPolySys,lbd,ubd,param);
            
            PUI(r_PUI,c_PUI) = POP.xVectL(r_PUI);
    end
end
theta_min = PUI(:,1);
theta_max = PUI(:,2);
theta_c = (theta_max+theta_min)/2;
acc = 100*(theta_max-theta_min)./theta_c;

table(theta_min,theta_c,theta_max,acc)

%%
theta_c = theta_c';
G = tf(theta_c(2:3),[1 theta_c(1)],Ts);
x = 0.25*u + theta_c(4)*u.^3;

y_sim = lsim(G,x);

figure
plot(y,'b')
hold on
plot(y_sim,'r--')
grid on


