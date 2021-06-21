function [xd,yd] = extractData(xm,ym)
% extractData takes two symbolic differential equations with variables x 
% and y and tropical coefficients and converts them to two structueres that
%contain important Data for tropicalization. 

%The struct xd contains:
%xd.l=The number of terms in the differential equation.
%xd.signs= A vector containing the non-tropicalized coefficients/signs
%of the terms in the xm differential equation.
%xd.alphas=A vector containing the tropicalized coefficients of the terms. 
%xd.degrees=A 2 \times xd.l matrix that contains the x and y degree of each
%term. 

%yd is identical to xd, but for the other differential equation. 


syms x y 

temp=char(xm);
%Expanding the terms (if there is more than 1) and putting them in an array
if(contains(temp,'+') || count(temp(2:end),'-')-count(temp,'exp(-') > 0 )
    tempx=children(expand(xm));
    tempx=simplify([tempx{:}]);
else
    tempx=xm;
end

temp=char(ym);
if(contains(temp,'+') || count(temp(2:end),'-')-count(temp,'exp(-') > 0 )
    tempy=children(expand(ym));
    tempy=simplify([tempy{:}]);
else
    tempy=ym;
end


%Inserting the length
xd.l=length(tempx);
yd.l=length(tempy);


%Extracting the signs
xd.signs=ones(1,xd.l);
yd.signs=ones(1,yd.l);
for i=1:xd.l
    xd.signs(i)=coeffs(tempx(i));
   
end

for i=1:yd.l
    yd.signs(i)=coeffs(tempy(i));
end

%Extracting the alphas
xd.alphas=ones(xd.l,1);
yd.alphas=ones(yd.l,1);
for i=1:xd.l
    temp=char(tempx(i));
    if(contains(temp,'exp'))
        temps=temp(strfind(temp,'exp')+4:strfind(temp,'eps')-2);
        xd.alphas(i)=double(sym(erase(temps,'(')));
    else
        xd.alphas(i)=0;
    end
end


for i=1:yd.l
    temp=char(tempy(i));
    if(contains(temp,'exp'))
        temps=temp(strfind(temp,'exp')+4:strfind(temp,'eps')-2);
        yd.alphas(i)=double(sym(erase(temps,'(')));
    else
        yd.alphas(i)=0;
    end
end

%Extracting the degrees
xd.degrees=[polynomialDegree(tempx,x); polynomialDegree(tempx,y)];
yd.degrees=[polynomialDegree(tempy,x); polynomialDegree(tempy,y)];


end




