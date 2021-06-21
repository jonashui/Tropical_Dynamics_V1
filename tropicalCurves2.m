function [out, out2] = tropicalCurves2(xm,ym, uint, vint)
%Plots the tropical Curves of the system of differential equations given by
%xm and ym. xm and ym are symbolic differential equations with tropicalised
%constants. uint and vint are the plotting intervals.

%Output is just for additional testing. 

%Extracting data from the differential equations
syms x y 
[xd,yd]=extractData(xm,ym);
syms u v 


%Computing the exponents that we compare against eachother. 
L1=ones(1, xd.l)*u;
L2=ones(1, yd.l)*u;
for i=1:xd.l
    L1(i)=xd.alphas(i)+(xd.degrees(1,i)-1)*u+xd.degrees(2,i)*v;
end

for i=1:yd.l
    L2(i)=yd.alphas(i)+yd.degrees(1,i)*u+(yd.degrees(2,i)-1)*v;
end
%Turn it into one vector. 
L=[L1,L2];
L

Lsigns=[xd.signs, yd.signs];


Lines={};
Vert={};
L1={};
L2={};

V1={};
V2={};
%Setting all exponents equal to eachother to generate all lines
%Also saving which exponents are used for each line
for i=1:length(L)
    for j=i+1:length(L)
        if(  contains(char(simplify(L(i)==L(j))),'v') )
            Lines=[Lines, solve(L(i)==L(j),v)];
            L1=[L1,i];
            L2=[L2,j];
            
        elseif(contains(char(simplify(L(i)==L(j))),'u') && ~contains(char(simplify(L(i)==L(j))),'true'))
            Vert=[Vert, solve(L(i)==L(j),u)];
            V1=[V1,i];
            V2=[V2,j];
        end
    end
end    


L1=cell2mat(L1);
L2=cell2mat(L2);
V1=cell2mat(V1);
V2=cell2mat(V2);

hold on

ineq=ones(1,length(L)-2)*u;

%Computing which lines are active by checking if the two exponents set
%equal to eachother are the largest ones. 
ALines=zeros(length(Lines),5)*u;
for i=1:length(Lines)
    tempL=L;
    tempL([L1(i),L2(i)])=[];
    
    
    ustart=uint(1);
    uend=uint(2);
    
    v=Lines(i);
   
    for j=1:length(tempL)
        ineq(j)=L(L1(i)) > tempL(j); 
        temp=L(L1(i)) == tempL(j);
       
        if(contains(char(temp),'v') && simplify(subs(temp,v)))
             ineq(j)=symtrue;
        elseif(simplify(L(L1(i)))== simplify(tempL(j)))
             ineq(j)=symtrue; 
        end
    end
  
    ineq=subs(ineq,v);
   
    ineq=simplify(ineq);
    ineq(ineq==in(u, 'real'))=symtrue;
        
   
    
    if(ismember(symfalse,ineq))
        continue
    end  
    
    for j=1:length(tempL)
        if(ineq(j)==symtrue)
            continue
        end
        u0=double(solve(lhs(ineq(j))==rhs(ineq(j)),u));
       
        
        c1=double(subs(ineq(j),u, ustart));
        c2=double(subs(ineq(j),u, uend));
        c3=(u0 > ustart && u0 < uend);
       
        
        if(~c1 && ~c2)
            uend=0;
            ustart=0;
            break
        elseif(c1  && c2)
            continue
        elseif(c1 && c3)
           uend=u0;
        elseif(c2 && c3)
            ustart=u0;
        end
        
        
    end
    
   
  
    ALines(i,:)=[Lines(i), ustart, uend, subs(Lines(i),u,ustart), subs(Lines(i),u,uend)];  
    
    if(ALines(i,4) < vint(1))
        ALines(i,2)=solve(Lines(i)==vint(1),u);
        ALines(i,4)=subs(Lines(i),u,ALines(i,2));
    elseif(ALines(i,4) > vint(2))
        ALines(i,2)=solve(Lines(i)==vint(2),u);
        ALines(i,4)=subs(Lines(i),u,ALines(i,2));
    elseif(ALines(i,5) < vint(1))
        ALines(i,3)=solve(Lines(i)==vint(1),u);
        ALines(i,3)=subs(Lines(i),u,ALines(i,3));
    elseif(ALines(i,5) > vint(2))
        ALines(i,3)=solve(Lines(i)==vint(2),u);
        ALines(i,5)=subs(Lines(i),u,ALines(i,3));
        
    end
     %Plotting the active lines
     plot([ALines(i,2) ALines(i,3)], [subs(Lines(i),u,ALines(i,2)) subs(Lines(i),u,ALines(i,3))])
end



TF2=double(subs(simplify(ALines(:,2)==0),u,-100));
TF3=double(subs(simplify(ALines(:,3)==0),u,-100));

ALines(TF2 & TF3,:)=[];

syms u v

%Again which lines are active, now for the vertical ones. 
AVert=zeros(length(Vert),3)*u;
for i=1:length(Vert)
    tempL=L;
    tempL([V1(i),V2(i)])=[];
    
    vstart=vint(1);
    vend=vint(2);
    
    u=double(Vert(i));
    
    for j=1:length(tempL)
        ineq(j)=L(V1(i)) > tempL(j);
        temp=L(V1(i)) == tempL(j);
        if(contains(char(temp),'u') && simplify(subs(temp)))
             ineq(j)=symtrue;
        elseif(simplify(L(V1(i)))== simplify(tempL(j)))
             ineq(j)=symtrue; 
        end   
    end
    
    ineq=subs(ineq);
    ineq=simplify(ineq);
    
    ineq(ineq==in(v, 'real'))=symtrue;
    
   
    
    if(ismember(symfalse,ineq))
        continue
    end  
    
    for j=1:length(tempL)
        if(ineq(j)==symtrue)
            continue
        end
        v0=double(solve(lhs(ineq(j))==rhs(ineq(j)),v));
        
        c1=double(subs(ineq(j),v, vstart));
        c2=double(subs(ineq(j),v, vend));
        c3=(v0 > vstart && v0 < vend);
        
        if(~c1 && ~c2)
            vend=0;
            vstart=0;
            break
        elseif(c1  && c2)
            continue
        elseif(c1 && c3)
           vend=v0;
        elseif(c2 && c3)
           vstart=v0;
        end
    end
    %Plotting the active lines
    plot([Vert(i) Vert(i)], [vstart vend])
    AVert(i,:)=[Vert(i), vstart, vend];
end
syms u v

TF2=double(simplify(AVert(:,2)==0));
TF3=double(simplify(AVert(:,3)==0));

AVert(TF2 & TF3,:)=[];



%Creating a grid to plot the arrows
X=[uint(1):(uint(2)-uint(1))/20:uint(2)];
Y=[vint(1):(vint(2)-vint(1))/20:vint(2)];

Grid=subs(L,u,X);
Grid=subs(Grid,v,Y);

Grid=double(reshape(Grid,[length(X),length(Y),length(L)]));

%Checking which exponent is max for each point in the grid. 
[m,~]=max(Grid, [], 3);

U=zeros(length(Y),length(X));
V=zeros(length(Y),length(X));
%Looping over all points and setting the arrows to the sign(s) of the
%largest one(s). 
for i=1:length(X)
    for j=1:length(Y)
        ii=find(m(j,i)==subs(L,{u,v},{X(i),Y(j)}));
        
        for k=1:length(ii)
            if(ii(k) > xd.l)
               V(j,i)=Lsigns(ii(k));
            else
               U(j,i)=Lsigns(ii(k));
            end
        end
    end
end


% 
% U=Lsigns(ii);
% U(ii > xd.l)=0;
% 
% V=Lsigns(ii);
% V(ii <= xd.l)=0;

% for i=1:length(X)
%     for j=1:length(Y)
%         for k=1:size(ALines,1)
%             if(point_to_line(X(i),Y(j),ALines(k,2),ALines(k,4),ALines(k,3),ALines(k,5)) < 0.2)
%                 U(j,i)=0;
%                 V(j,i)=0;
%                 break
%             end
%         end 
%             
%          for k=1:size(AVert,1)
%             if(point_to_line(X(i),Y(j),AVert(k,1),AVert(k,2),AVert(k,1),AVert(k,3)) < 0.2)
%                 U(j,i)=0;
%                 V(j,i)=0;
%                 break
%             end
%         end
%     end
% end




%%%%%

%Storing all intersections to create polygons. 
xp=ALines(:,2:3);
xp=[xp(:)];

xp2=AVert(:,1);
xp2=[xp2(:)]';
xp2=repelem(xp2, 2);

xp=double([xp', xp2]);

yp=ALines(:,4:5);
yp=yp(:);
yp2=AVert(:,2:3);
yp2=yp2';
yp2=[yp2(:)];

yp=double([yp', yp2']);

p=[xp;yp];

p=[p, [uint(1);vint(1)], [uint(1);vint(2)], [uint(2);vint(1)], [uint(2);vint(2)]];
    
p=unique(p','rows')';
size(p);


polyx=cell(length(L),1);
polyy=cell(length(L),1);    

%Adding the intersections to the polygons where the exponents are largest.
for i=1:size(p,2)
    test=subs(L,{u,v},{p(1,i),p(2,i)});
    m=max(test);
    ind=find(test==m);
    
    for j=1:length(ind)
       polyx{ind(j)}=[polyx{ind(j)}, p(1,i)];
       polyy{(ind(j))}=[polyy{ind(j)}, p(2,i)];
    end 
end

poly=cell(length(L),1);

%Creating the polygons
for i=1:length(L)
    if(length(polyx{i}) < 3)
        continue
    end
    x=polyx{i};
    y=polyy{i};
    cx = mean(x);
    cy = mean(y);
    
    a = atan2(y - cy, x - cx);
    
    [~, order] = sort(a);
    
    x=x(order);
    y=y(order);
    
    %Making the polygons slightly smaller.
    for j=1:length(x)
        x(j)=x(j)-0.12*(x(j)-cx);
        y(j)=y(j)-0.12*(y(j)-cy);
    end
    
    polyx{i} = x;
    polyy{i} = y;
    
    
    poly{i}=polyshape(polyx{i},polyy{i});
end


poly=poly(~cellfun('isempty',poly));

% figure
% hold on
% for i=1:length(poly)
%     plot(poly{i})
% end
% 
% figure
% hold on

[XG,YG]=meshgrid(X,Y);
XG=reshape(XG,[size(XG,1)^2,1]);
YG=reshape(YG,[size(YG,1)^2,1]);

%Only plot the arrows inside the polygons, to avoid arrow crossing lines.
in2=cell(length(poly),1);
for i=1:size(poly,1)
     in2{i}=isinterior(poly{i},XG,YG);
end


in2=[in2{:}];
in2=all(in2==0,2);
in2=reshape(in2, size(U));

U(in2==1)=0;
V(in2==1)=0;

%Plotting the arrows
quiver(X,Y,U,V, 0.4, 'Color', '#6A5ACD')

xlabel('u')
ylabel('v')

axis([uint(1) uint(2) vint(1) vint(2)])


% Filippov flow for the lines
for i=1:size(ALines,1)
   for j=0.4:0.46:1.9  
    r=1/2*([ALines(i,3),ALines(i,5)]-[ALines(i,2),ALines(i,4)]);
    x=[ALines(i,2),ALines(i,4)]+j*r;
    a=Filippov(ALines(i,1),x, L, xd, yd,r);
    
    quiver(x(1),x(2),a(1),a(2),0.3,'LineWidth', 2, 'Color', 'k', 'MaxHeadSize', 2);
   
   end
end

% Filippov for the vertical lines
for i=1:size(AVert,1)
   for j=0.4:0.4:1.6 
    r=1/2*([AVert(i,1),AVert(i,3)]-[AVert(i,1),AVert(i,2)]);
    x=[AVert(i,1),AVert(i,2)]+r*j;
    a=FilippovV(x, L, xd, yd,r);
   
    quiver(x(1),x(2),a(1),a(2),0.3,'LineWidth', 2, 'Color', 'k', 'MaxHeadSize', 2);
    
   end
end

%Optional output for additional testing. 
out=ALines;
out2=L;

