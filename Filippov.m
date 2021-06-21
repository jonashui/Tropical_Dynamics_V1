function a = Filippov(Line,x,L,xd,yd,r)
%Computes the Filippov flow (and looks at the next largest in the opposite
%direction if the Filippov flow is equal 0. Input is a line in the format 
%v=a*u+b. x is the point on the line and r is the direction vector of the line.

%xd and yd are the data given from the structs in extractData. 

syms u v
H1=coeffs(Line-v,u);
H2=coeffs(Line-v,v);

%Computing the normal vector.
if(has(Line,u))
    H_x=double([H1(2),H2(2)]);
elseif(Line==0)
    H_x=[0,double(H2)];
else
    H_x=[0,double(H2(2))];
end

%Extracting the arrows from both sides of the equation
xplus=subs(L,{u,v},{x+0.1*H_x});
m=max(xplus);
ii=find(xplus==m);
F1=[0,0];
for j=1:length(ii)
    if(ii(j) > xd.l)
        F1=F1+[0,yd.signs(ii(j)-xd.l)];
    else
        F1=F1+[xd.signs(ii(j)),0];
    end    
end


xminus=subs(L,{u,v},{x-0.1*H_x});
m=max(xminus);
ii=find(xminus==m);
F2=[0,0];
for j=1:length(ii)
    if(ii(j) > xd.l)
        F2=F2+[0,yd.signs(ii(j)-xd.l)];
    else
        F2=F2+[xd.signs(ii(j)),0];
    end    
end


%If the vectors are opposite, we have to look in the other direction. 
if(F1(2)==0 && F2(2)==0 && sign(F1(1)) ~= sign(F2(1)) )
    vec=subs(L(xd.l+1:end),{u,v},{x});
    m=max(vec);
    i=find(vec==m);
    temp=[0,0];
    for j=1:length(i)
        yd.signs(i(j));
        temp=temp+[0,yd.signs(i(j))];
    end
    a=dot(temp,r)/(norm(r)^2)*r;
        return   
elseif(F1(1)==0 && F2(1)==0 && sign(F1(2)) ~= sign(F2(2)))
    vec=subs(L(1:xd.l),{u,v},{x});
    m=max(vec);
    i=find(vec==m);
    temp=[0,0];
    for j=1:length(i)
        yd.signs(i(j));
        temp=temp+[xd.signs(i(j)),0];
    end
    a=dot(temp,r)/(norm(r)^2)*r;
    return    
end




%Standard Filippov flow
d=dot(H_x,(F2-F1));

if(d ==0)
    a=[0,0];
else
    lambda=dot(H_x,F2)/d;
    if(lambda <= 0 || lambda ==1)
        a=[0,0];
    else
         a=lambda*F1+(1-lambda)*F2;
    end
end







