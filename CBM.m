function C= CBMWX8(Df, Dp, N)

digits(10);

cj = 1;
cy = 2;
cz = 8;
L = 52;

galpha=1;
gbeta=0.5;
betaPowAlpha=gbeta^galpha;%gbeta^galpha
gammaAlpha=gamma(galpha);%gamma(galpha)
gammaPara=betaPowAlpha/gammaAlpha;%gbeta^galpha/gamma(galpha)

dividedDp = [];
for i=0:N
    dividedDp=[dividedDp,i*Dp/N];
end

alpha=[];%役龄回退因子
for i=0:L
    alpha=[alpha,i / (50 * i + 5)];
%     alpha=[alpha,0];
end

pj = [1];
py = [0];
pz = [0];

syms x y;
omega={0};
omegashu={@(x)0*x,@(x)0*x,@(x)0*x};

C=0;

for i=1:L
    otshu1=0;%0<x<Dp
    otshu2=0;%Dp<x<Df
    otshu3=0;%Df<x
    if i==1
        ot=vpa(piecewise(x<=0,0,x>0,gammaPara*x^(galpha-1)*exp(-gbeta*x)));
        childs=children(ot);
        childsSize=size(childs);
        childsSize=childsSize(1);
        otshu1=childs(2,1);
        otshu1=otshu1{1,1};
        otshu1=matlabFunction(otshu1);
        if childsSize==2
            otshu2=otshu1;
            otshu3=otshu1;
        end
        omega=[omega,ot];
    else
        omegay=subs(omega(i),x,y);
        
        f1=vpa(simplify(omegay*gammaPara*x^(galpha-1)*exp(-gbeta*(x-y))));
        f2=vpa(simplify(omegay*gammaPara*x^(galpha-1)*exp(-gbeta*(x))));
        
        ot1=vpa(simplify(piecewise(x<=0,0,(0<x)&(x<Dp),pj(i)*int(f1,y,0,x),x>=Dp,pj(i)*int(f1,y,0,Dp))));
        ot2=vpa(simplify(piecewise(x<=0,0,x>0,pj(i)*int(f2,y,Dp,Inf))));
        ot3=vpa(simplify(piecewise(x<=0,0,(0<x)&(x<Df),(1-pj(i))*int(f1,y,0,x),x>=Df,(1-pj(i))*int(f1,y,0,Df))));
        ot4=vpa(simplify(piecewise(x<=0,0,x>0,(1-pj(i))*int(f2,y,Df,Inf))));
        
        ot=vpa(simplify(ot1+ot2+ot3+ot4));
        childs=children(ot);
        childsSize=size(childs);
        childsSize=childsSize(1);
        otshu1=childs{2,1};
        otshu1=matlabFunction(otshu1);
        if childsSize==2
            otshu2=otshu1;
            otshu3=otshu1;
        elseif childsSize==3
            otshu2=otshu1;
            otshu3=childs{3,1};
            otshu3=matlabFunction(otshu3);
        elseif childsSize==4
            otshu2=childs{3,1};
            otshu2=matlabFunction(otshu2);
            otshu3=childs{4,1};
            otshu3=matlabFunction(otshu3);
        elseif childsSize==5
            otshu2=childs{4,1};
            otshu2=matlabFunction(otshu2);
            otshu3=childs{5,1};
            otshu3=matlabFunction(otshu3);
        end
        omega=[omega,ot];
%         omega(i+1)
    end
    omegashu{i+1,1}=otshu1;
    omegashu{i+1,2}=otshu2;
    omegashu{i+1,3}=otshu3;

%     otshu1
%     otshu2
%     otshu3
    
    intoh1=integral(otshu1,0,Dp);
    intoh2=integral(otshu2,Dp,Df);
    intoh3=integral(otshu3,Df,Inf);
    disp(['t:',num2str(i)])
    disp(['0<x<Dp:',num2str(intoh1),'Dp<=x<Df:',num2str(intoh2),'Df<=x:',num2str(intoh3),'sum:',num2str(intoh1+intoh2+intoh3)])
    
    if i<N
        pj=[pj,0];
        py=[py,0];
        pz=[pz,integral(otshu3,Df,Inf)];
    else
        pj1=pj(i-N+1)*(integral(omegashu{i-N+1,2},Dp,Df)+pz(i-N+1));
        for j=i-N+2:i
            pj1=pj1*(1-pj(j))*(1-pz(j));
%             pj1=pj1*(1-pz(j));
        end
        
        if i==N
            pj2=pj(i-N+1);
        else
            pj2=pj(i-N+1)*integral(omegashu{i-N+1,1},dividedDp(1),dividedDp(2));
        end
        for j=i-N+2:i
            pj2=pj2*(1-pj(j))*(1-pz(j));
%             pj2=pj2*(1-pz(j));
        end
        
        pj3=(1-pj(i-N+1))*pz(i-N+1);
        for j=i-N+2:i
            pj3=pj3*(1-pj(j))*(1-pz(j));
%             pj3=pj3*(1-pz(j));
        end
        
        pj4=0;
        for j=1:N-2
            pj4temp=pj(i-N+j+1)*integral(omegashu{i-N+1,1},dividedDp(j+1),dividedDp(j+2));
            for k=i-N+j+2:i
                pj4temp=pj4temp*(1-pj(k))*(1-pz(k));
%                 pj4temp=pj4temp*(1-pz(k));
            end
            pj4=pj4+pj4temp;
        end
        
        pj5=pj(i)*integral(omegashu{i,1},dividedDp(N),dividedDp(N+1));
        pjtemp=pj1+pj2+pj3+pj4+pj5;
        pj=[pj,pjtemp];
        py=[py,pjtemp*integral(otshu2,Dp,Df)];
        pz=[pz,integral(otshu3,Df,Inf)];
    end
    
    C=C+cj*pj(i+1)+cy*py(i+1)+cz*pz(i+1);
    disp(['检测概率:',num2str(pj(i+1))])
    disp(['预防维修概率:',num2str(py(i+1))])
    disp(['故障后维修概率',num2str(pz(i+1))])
    disp(['总成本:',num2str(C)])
    
    
        
    
%     pj1=0;
%     pj2=0;
%     pj3=0;
%     
%     for j=1:min(i,N)
%         if j==N
%             pj1=1;
%             pj2=1;
% %             pj1=pj1*pj(i-N+1)*(py(i-N+1)+pz(i-N+1));
%             pj1=pj1*pj(i-N+1)*(integral(omegashu{i-N+1,2},Dp,Df)+integral(omegashu{i-N+1,3},Df,Inf));
%             pj2=pj2*(1-pj(i-N+1))*pz(i-N+1);
%             for k=1:(N-1)
%                 pj1=pj1*(1-pj(i-k+1));
%                 pj1=pj1*(1-pz(i-k+1));
%                 pj2=pj2*(1-pj(i-k+1));
%                 pj2=pj2*(1-pz(i-k+1));
%             end
%         end
%         pj3temp=1;
%         if j==N&&i==N
%             pj3temp=pj3temp*pj(i-j+1);
%         else
%             pj3temp=pj3temp*pj(i-j+1)*integral(omegashu{i-j+1,1},dividedDp(N-j+1),dividedDp(N-j+2));
%         end
%         if j>1
%             for k=1:(j-1)
%                 pj3temp=pj3temp*(1-pj(i-k+1));
%                 pj3temp=pj3temp*(1-pz(i-k+1));
%             end
%         end
%         pj3=pj3+pj3temp;
%     end
% 
%     pj=[pj,pj1+pj2+pj3];
%     py=[py,pj(i+1)*integral(otshu2,Dp,Df)];
%     pz=[pz,integral(otshu3,Df,Inf)];
    
%     C=C+cj*pj(i+1)+cy*py(i+1)+cz*pz(i+1);
%     disp(['检测概率:',num2str(pj(i+1))])
%     disp(['预防维修概率:',num2str(py(i+1))])
%     disp(['故障后维修概率',num2str(pz(i+1))])
%     disp(['总成本:',num2str(C)])
end