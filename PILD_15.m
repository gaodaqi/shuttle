t0=clock;

load D:\highimbalance\shuttle\trainbpca15.txt;
load D:\highimbalance\shuttle\ntr.txt;

n1=size(trainbpca15);
dim=n1(2)-2;
m=mean(mean(trainbpca15(:,1:dim)));
if abs(m)<0.1
   m=0.1;
end
m=1;
k0=1;
c=size(ntr,1);

for r=1:c-1,
    for t=r+1:c,
        ntr1=[ntr(r,1),ntr(r,2),ntr(t,1),ntr(t,2)];
        p0=ntr1(2)-ntr1(1)+1;
        q0=ntr1(4)-ntr1(3)+1;


        p100=0;
        q100=0;

        d00=ones(p0,1);
        d00(p0+1:p0+q0,1)=-ones(q0,1);
        d0=d00;        

        d=d0;
        a0=[];
        a0(:,2:dim+1)=trainbpca15(ntr1(1):ntr1(2),1:dim);
        a0(p0+1:p0+q0,2:dim+1)=trainbpca15(ntr1(3):ntr1(4),1:dim);
        a0(:,1)=m*ones(p0+q0,1);
        a2=a0;


        p=p0;
        q=q0;

        sign=1;
        n=1;
        n0=1;

        n000(k0)=0;

        nmax=150;

        while n<nmax&sign==1,

              if p==1&q==1
                 minnumerror(k0)=0;
                 w1(k0,1:dim)=a2(1,2:dim+1)-a2(2,2:dim+1);
                 mu=a2(1,2:dim+1)*w1';
                 mu(2)=a2(2,2:dim+1)*w1';
                 theta1(k0)=-(mu(1)+mu(2))/2;
                 break;
              end


              n000(k0)=n000(k0)+1;
              if n==1
                 atrans=a2'*a2;
                 arank=rank(atrans);
                 adet=det(atrans);
                 nn00=1;

                 if abs(adet)<1.0e-10|arank<dim+1
                    atrans=atrans+0.00005*eye(dim+1);
                    nn00=-1;
%                  fprintf('det(%d %d)=%g, rank(%d %d)=%g, Matrix(%d,%d) is close to singular\n',r,t,adet,r,t,arank, r,t);
                 end
                 invmat=inv(atrans)*a2';
              end

              w=invmat*d;
%              theta(n)=w(1)-(p-q)/(2*(p+q));

              project=a2(:,2:dim+1)*w(2:dim+1);
              mu(1)=mean(project(1:p));
              mu(2)=mean(project(p+1:p+q));

%//////////////////////////////////////////////////

              p1=0;
              re1=[];
              for i=1:p,
                  if (project(i)<mu(1))&(project(i)>mu(2))
                     p1=p1+1;
                     re1(p1)=project(i);
                  end
              end

              q1=0;
              re2=[];
              for i=p+1:p+q,
                  if (project(i)<mu(1))&(project(i)>mu(2))
                     q1=q1+1;
                     re2(q1)=project(i);
                  end
              end

%             mmaxmin=min(re1);
              if p1==0
                 mmaxmin(1)=mu(1);
              else
                 mmaxmin(1)=max(re1);
              end
              if q1==0
                 mmaxmin(2)=mu(2);
              else
                 mmaxmin(2)=max(re2);
              end

%              maverage=mean(re1);
              if p1==0
                 maverage(1)=mu(1);
              else
                 maverage(1)=mean(re1);
              end 
              if q1==0
                 maverage(2)=mu(2);
              else
                 maverage(2)=mean(re2);
              end              

              theta(n)=-(mmaxmin(1)+mmaxmin(2))/(m*2);

%//////////////////////////////////////////////////

%              theta(n)=-(mu(1)+mu(2))/(m*2);
%              theta(n)=-(p*mu(1)+q*mu(2))/(m*(p+q));
%              theta(n)=w(1)-(p-q)/(m*(p+q));
%              theta(n)=w(1);
              theta(n)=-(maverage(1)+maverage(2))/(m*2);

% Reassigned the desired outputs according to theta

              numerror1=0;
              numerror2=0;
              project0=a0(:,2:dim+1)*w(2:dim+1);
              result0=project0+m*theta(n);

              for i=1:p0,
                  if result0(i)<0
                     numerror1=numerror1+1;
                  end
              end
              for i=p0+1:p0+q0,
                  if result0(i)>0
                     numerror2=numerror2+1;
                  end
              end

              numerror(n)=numerror1+numerror2;

              fprintf('numerror(%d)=%d %d %d %d\n',n,numerror(n),numerror1,numerror2,nn00);

              result1=project+m*theta(n);
              mmin=10000000;


              for i=1:p0,
                  if result0(i)>0
                     if result0(i)<mmin
                        mmin=result0(i);
                     end
                  end
              end

              if numerror1==p0
                 mmin=-((3*mu(1)+mu(2))/4+theta(n));
              end

              for i=1:p,
                  d(i,1)=result1(i);
                  if result1(i)<0
                     d(i,1)=mmin/2;
                  end
              end

              mmax=-10000000;

              for i=p0+1:p0+q0,
                  if result0(i)<0
                     if mmax<result0(i)
                        mmax=result0(i);
                     end
                  end
              end

              for i=p+1:p+q,
                  d(i,1)=result1(i);
                  if result1(i)>0
                     d(i,1)=mmax/2;
                  end
              end

              dminmax(n,1)=mmin;
              dminmax(n,2)=mmax;

% list 5 continuous weights, thresholds and errors:

              if numerror(n)==0
                 minnumerror(k0)=0;
                 w1(k0,:)=w(2:dim+1);
                 theta1(k0)=theta(n);
                 break;
              end
              
              if n==1
                 w0=zeros(5,dim);
                 theta0=zeros(5);
                 numerror0=zeros(5,1);
              end

              if n<=5
                 w0(n,:)=w(2:dim+1);
                 theta0(n)=theta(n);
                 numerror0(n)=numerror(n);
              else
                 w0(1:4,:)=w0(2:5,:);
                 w0(5,:)=w(2:dim+1);
                 theta0(1:4)=theta0(2:5);
                 theta0(5)=theta(n);
                 numerror0(1:4)=numerror0(2:5);
                 numerror0(5)=numerror(n);
              end

              sign1=0;
              if (n>=5)&(numerror0(1)<min(numerror0(2:5)))
                 sign1=1;
              else
                 n=n+1;
              end

% Set up the stop condition
              
              if (n0==1)&(sign1==1)
                  sign=0;
                  minnumerror(k0)=numerror0(1);
                  w1(k0,:)=w0(1,:);
                  theta1(k0)=theta0(1);
              end

%///////////////////////////////////////////////////////////////
              if (n0>1)&(sign1==1)
                 sign=0;
                 if numerror0(1)<minnumerror(k0)
                    minnumerror(k0)=numerror0(1);
                    w1(k0,:)=w0(1,:);
                    theta1(k0)=theta0(1);
                 end
              end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
              if (n0==1)&(sign1==0)&(n==nmax)
                 sign=0;
                 sign1=1;

                 for ii=1:5,
                     numerror_re(ii)=numerror0(6-ii);
                     w_re(ii,:)=w0(6-ii,:);
                     theta_re(ii)=theta0(6-ii);
                 end

                 [minvalue,minindex]=min(numerror_re);
                  minnumerror(k0)=numerror_re(minindex);
                  w1(k0,:)=w_re(minindex,:);
                  theta1(k0)=theta_re(minindex);

              end              
%///////////////////////////////////////////////////////////////
              if (n0>1)&(sign1==0)&(n==nmax)
                  sign=0;
                  sign1=1;

                 for ii=1:5,
                     numerror_re(ii)=numerror0(6-ii);
                     w_re(ii,:)=w0(6-ii,:);
                     theta_re(ii)=theta0(6-ii);
                 end

                 [minvalue,minindex]=min(numerror_re);
                 if numerror_re(minindex)<minnumerror(k0)
                    minnumerror(k0)=numerror_re(minindex);
                    w1(k0,:)=w_re(minindex,:);
                    theta1(k0)=theta_re(minindex);
                 end

              end
%//////////////////////////////////////////////////////////////// 
              
% Sample decomposition

              if ((sign1==1)&((p/dim>=8)|(q/dim>=8)))
%              if ((sign1==1)&((p/dim>=8)|(q/dim>=8)))|((n==nmax)&((p/dim>=8)|(q/dim>=8)))
                  sign=1;
                  project00=[];
                  if n==nmax
                     project00=a2(:,2:dim+1)*(w_re(minindex,:))';
                  else
                     project00=a2(:,2:dim+1)*(w0(1,:))';
                  end
                  mu0=mean(project00(1:p));
                  mu0(2)=mean(project00(p+1:p+q));
                  a20=[];
                  d0=[];

                  if p/dim>8
                     p10=0;
                     for i=1:p,
                         if (project00(i)<mu0(1))&(project00(i)>mu0(2))
                             p10=p10+1;
                             a20(p10,:)=a2(i,:);
                             d0(p10,1)=d(i,1);
                         end
                     end
                  else
                     p10=p;
                     a20=a2(1:p10,:);
                     d0=d(1:p10,1);
                  end

                  if q/dim>8
                     q10=0;
                     for i=p+1:p+q,
                         if (project00(i)<mu0(1))&(project00(i)>mu0(2))
                             q10=q10+1;
                             a20(p10+q10,:)=a2(i,:);
                             d0(p10+q10,1)=d(i,1);
                         end
                     end
                  else
                     q10=q;
                     a20(p10+1:p10+q10,:)=a2(p+1:p+q10,:);
                     d0(p10+1:p10+q10,1)=d(p+1:p+q10,1);
                  end

                  if p==p10&q==q10
                     break;
                  end
                  
                  a2=a20;
                  d=d0;
                  if p10>0
                     p=p10;
                  else
                      break;
                  end

                  if q10>0
                     q=q10;
                  else
                     break;
                  end

                  fprintf('(p, q)=%d, %d\n',p,q);
                  n=1;
                  n0=n0+1;

                  if p10==p100&q10==q100
                     sign=0;
                  else
                     p100=p10;
                     q100=q10;
                  end

                  if n0>=log(p0+q0)/log(2)+1
                     sign=0;
                  end
                  numerror=[];
                  numerror0=[];
             end
        end

        fprintf('n(%d,%d)=%d\n\n',r,t,n);
        minnumerror
        k0=k0+1;
    end
end
t=etime(clock,t0)
n000
sum(n000)
%n0

fid=fopen('D:\highimbalance\shuttle\invtca_d\w.txt','w');
    fprintf(fid,'%g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n',w1');
fclose(fid);
fid=fopen('D:\highimbalance\shuttle\invtca_d\theta.txt','w');
    fprintf(fid,'%g\n',theta1');
fclose(fid);
fid=fopen('D:\highimbalance\shuttle\invtca_d\m.txt','w');
    fprintf(fid,'%g\n',m');
fclose(fid);

t0=clock;

load D:\highimbalance\shuttle\invtca_d\w.txt;
load D:\highimbalance\shuttle\invtca_d\theta.txt;
load D:\highimbalance\shuttle\invtca_d\m.txt;

% Total, average, geometric accuracy; AUC:
% For the trainbpca15ing set.

k=0;
for r=1:c-1,
    for t=r+1:c,
        re=[];
        a=[];
        k=k+1;
        p0=ntr(r,3);
        q0=ntr(t,3);
        a=trainbpca15(ntr(r,1):ntr(r,2),1:dim);
        a(p0+1:p0+q0,:)=trainbpca15(ntr(t,1):ntr(t,2),1:dim);
        re=a*(w(k,:))';

% AUC
        y=[];
        [y(:,1),y(:,2)]=sort(re,'descend');

        b=[0 0];
        j10=0;
        j00=0;
        j1=0;
        j0=0;

        for i=1:p0+q0,
            if y(i,2)<=p0
               j10=j10+1;
               j1=j10/p0;
            else
               j00=j00+1;
               j0=j00/q0;
            end
            b(i+1,1)=j0;
            b(i+1,2)=j1;
        end

        figure(1);
        plot(b(:,1),b(:,2));
        area_trainbpca15(r,t)=0;
        for i=2:p0+q0+1,
            area_trainbpca15(r,t)=area_trainbpca15(r,t)+b(i,2)*(b(i,1)-b(i-1,1));
        end
        area_trainbpca15(r,t);

% Total, average, geometric accuracy

        for j=1:size(theta,2),
            numtr0=[0 0];
            for i=1:p0,
                if re(i)+m*theta(k,j)<0
                   numtr0(1)=numtr0(1)+1;
                end
            end
            for i=p0+1:p0+q0,
                if re(i)+m*theta(k,j)>0
                   numtr0(2)=numtr0(2)+1;
                end
            end
            num_trainbpca15(c*(j-1)+r,t)=numtr0(1);
            num_trainbpca15(c*(j-1)+t,r)=numtr0(2);

            err_trainbpca15(c*(j-1)+r,t-1)=100*(1-(numtr0(1)+numtr0(2))/(p0+q0));
            err_trainbpca15(c*(j-1)+r,c+t-2)=100*((p0-numtr0(1))/p0+(q0-numtr0(2))/q0)/2;
            err_trainbpca15(c*(j-1)+r,2*c+t-3)=100*sqrt(((p0-numtr0(1))/p0)*((q0-numtr0(2))/q0));
            err_trainbpca15(c*(j-1)+r,3*c+t-4)=100*area_trainbpca15(r,t);
            err_trainbpca15(c*(j-1)+r,4*c+t-5)=(err_trainbpca15(c*(j-1)+r,t-1)+err_trainbpca15(c*(j-1)+r,c+t-2)+err_trainbpca15(c*(j-1)+r,2*c+t-3)+area_trainbpca15(r,t)*100)/4;        

        end
    end
end

err_trainbpca15

t=etime(clock,t0)

% For testbpca15 set.

load D:\highimbalance\shuttle\testbpca15.txt;
load D:\highimbalance\shuttle\nte.txt;

% For testbpca15 set.

k=0;
for r=1:c-1,
    for t=r+1:c,
        re=[];
        a=[];
        k=k+1;
        p0=nte(r,3);
        q0=nte(t,3);
        a=testbpca15(nte(r,1):nte(r,2),1:dim);
        a(p0+1:p0+q0,:)=testbpca15(nte(t,1):nte(t,2),1:dim);
        re=a*(w(k,:))';

% AUC

        y=[];
        [y(:,1),y(:,2)]=sort(re,'descend');

        b=[0 0];
        j10=0;
        j00=0;
        j1=0;
        j0=0;

        for i=1:p0+q0,
            if y(i,2)<=p0
               j10=j10+1;
               j1=j10/p0;
            else
               j00=j00+1;
               j0=j00/q0;
            end
            b(i+1,1)=j0;
            b(i+1,2)=j1;
        end

        figure(2);
        plot(b(:,1),b(:,2));
        area_testbpca15(r,t)=0;
        for i=2:p0+q0+1,
            area_testbpca15(r,t)=area_testbpca15(r,t)+b(i,2)*(b(i,1)-b(i-1,1));
        end
        area_testbpca15(r,t);

% Total, average, geometric accuracy

        for j=1:size(theta,2),
            numtestbpca150=[0 0];
            for i=1:p0,
                if re(i)+m*theta(k,j)<0
                   numtestbpca150(1)=numtestbpca150(1)+1;
                end
            end
            for i=p0+1:p0+q0,
                if re(i)+m*theta(k,j)>0
                   numtestbpca150(2)=numtestbpca150(2)+1;
                end
            end
            num_testbpca15(c*(j-1)+r,t)=numtestbpca150(1);
            num_testbpca15(c*(j-1)+t,r)=numtestbpca150(2);

            err_testbpca15(c*(j-1)+r,t-1)=100*(1-(numtestbpca150(1)+numtestbpca150(2))/(p0+q0));
            err_testbpca15(c*(j-1)+r,c+t-2)=100*((p0-numtestbpca150(1))/p0+(q0-numtestbpca150(2))/q0)/2;
            err_testbpca15(c*(j-1)+r,2*c+t-3)=100*sqrt(((p0-numtestbpca150(1))/p0)*((q0-numtestbpca150(2))/q0));
            err_testbpca15(c*(j-1)+r,3*c+t-4)=100*area_testbpca15(r,t);
            err_testbpca15(c*(j-1)+r,4*c+t-5)=(err_testbpca15(c*(j-1)+r,t-1)+err_testbpca15(c*(j-1)+r,c+t-2)+err_testbpca15(c*(j-1)+r,2*c+t-3)+area_testbpca15(r,t)*100)/4;
        end
    end
end

err_testbpca15

fid=fopen('D:\highimbalance\shuttle\invtca_d\err_trainbpca15.txt','w');
    fprintf(fid,'%g %g %g %g %g\n',err_trainbpca15');
fclose(fid);
fid=fopen('D:\highimbalance\shuttle\invtca_d\err_testbpca15.txt','w');
    fprintf(fid,'%g %g %g %g %g\n',err_testbpca15');
fclose(fid);

t=etime(clock,t0)

t0=clock;

%///////////////////////////////////////////////////////////////////////////%

%clear all;

% Majority votestbpca15
% trainbpca15ing set

t0=clock;

load D:\highimbalance\shuttle\trainbpca15.txt;
load D:\highimbalance\shuttle\ntr.txt;

load D:\highimbalance\shuttle\invtca_d\w.txt;
load D:\highimbalance\shuttle\invtca_d\theta.txt;
load D:\highimbalance\shuttle\invtca_d\m.txt;


% Majority votestbpca15
% trainbpca15 set

re=trainbpca15(:,1:size(w,2))*w';
nerr=[];

n0=size(re,1);
k1=size(ntr,1);
k2=k1-1;

for t=1:size(theta,2),
    re1=re+m*ones(n0,1)*(theta(:,t))';

    for i=1:n0,
        k3=1;
        k4=k2;
        b=0;
        for j=1:k2,
            a(j:k2,j)=(re1(i,k3:k4))';
            k3=k4+1;
            k4=k4+k2-j;
        end

        for j=1:k2,
            for k=j+1:k1,
                a(j,k:k)=-a(k-1,j);
            end
        end

        for q=1:k1,
            k=k2;
            for j=1:k2,
                if a(j,q)<0
                   k=k-1;
                end
            end
            restestbpca151(i,q)=k;
        end
    end

    nerr1=zeros(k1,k1);
    kk0=0;
    for r=1:k1,
        nerr(t,r)=0;
        for i=ntr(r,1):ntr(r,2),
            b0=restestbpca151(i,r);
            b1=restestbpca151(i,:);
            b1(r)=-1;
            [amax,nmax]=max(b1);
            if amax>b0
               nerr(t,r)=nerr(t,r)+1;
               nerr1(r,nmax)=nerr1(r,nmax)+1;
               kk0=kk0+1;
            end
        end
    end
end

nerr
nerr1
sum(nerr')
1-sum(nerr')/ntr(k1,2)

t=etime(clock,t0)

t0=clock;

% testbpca15 set
load D:\highimbalance\shuttle\testbpca15.txt;
load D:\highimbalance\shuttle\nte.txt;


nerr=[];
re=testbpca15(:,1:size(w,2))*w';

n0=size(re,1);
k1=size(nte,1);
k2=k1-1;
re1=[];
a=[];
for t=1:size(theta,2),
    re1=re+m*ones(n0,1)*(theta(:,t))';

    for i=1:n0,
        k3=1;
        k4=k2;
        b=0;
        for j=1:k2,
            a(j:k2,j)=(re1(i,k3:k4))';
            k3=k4+1;
            k4=k4+k2-j;
        end

        for j=1:k2,
            for k=j+1:k1,
                a(j,k:k)=-a(k-1,j);
            end
        end

        for q=1:k1,
            k=k2;
            for j=1:k2,
                if a(j,q)<0
                   k=k-1;
                end
            end
            restestbpca151(i,q)=k;
        end
    end

    kk0=0;
    nerr1=zeros(k1,k1);
    for r=1:k1,
        nerr(t,r)=0;
        for i=nte(r,1):nte(r,2),
            b0=restestbpca151(i,r);
            b1=restestbpca151(i,:);
            b1(r)=-1;
            [amax,nmax]=max(b1);
            if amax>b0
               nerr(t,r)=nerr(t,r)+1;
               nerr1(r,nmax)=nerr1(r,nmax)+1;
            end
        end
    end
end

nerr
nerr1
sum(nerr')
1-sum(nerr')/nte(k1,2)

t=etime(clock,t0)
