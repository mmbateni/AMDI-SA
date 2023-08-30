function [ cop_empri ] = empcop( w )
    %% Calculating Empirical Copula for GOFs; Empirical calculated based on
    %  Aghakouchk's work on Non-Parametric Approach. input must be a (n*2)
    %  matrix
    nn=size(w,1);
    bp=zeros(nn,1);
    for i=1:nn
    td=zeros(nn,3);
    td(w(:,1)<=w(i,1),1)=1;
    td(w(:,2)<=w(i,2),2)=1;
    td(:,3)=td(:,1).*td(:,2);
    bp(i)=sum(td(:,3));
    
    end
    % Can be replaced with Gringorten Or California plotting position
    % cop_empri(:,ui)=(bp-044)./(nn+0.12);
    % cop_empri(:,ui)=(bp)./(nn);
    cop_empri=(bp)./(nn+1);

end

