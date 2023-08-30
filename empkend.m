function [Kn] = empkend(x,tq)
n=size(x,1);
n_tq=length(tq);
w=zeros(n,1);
wu(:,1)=sort(x(:,1));
wu(:,2)=sort(x(:,2));
hi_ind=zeros(n,1);
hj_ind=zeros(n,1);
Kn=zeros(n_tq,1);
for k=1:n_tq
for j=1:n
for i=1:n
    hi_ind(i)=(x(i,1)<wu(j,1))&(x(i,2)<wu(j,2));
end
    w(j)=sum(hi_ind)/(n+1);
    hj_ind(j)=(w(j)<tq(k));
end
    Kn(k)=sum(hj_ind)/(n);
end
