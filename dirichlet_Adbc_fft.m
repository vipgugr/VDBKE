function [Ad,b,atAa,Aa,flipx] = dirichlet_Adbc_fft(x,y,m1,m2,lambda_C,C)
%Compute Ad, b, and function hanldes for XtX
if ~exist('lambda_C','var')
    lambda_C=0;
end
L=length(x);
Ad=zeros(m1,m2);
b=Ad;
flipx=cell(L,1);
Xr=cell(L,1);
X=cell(L,1);
for t=1:length(x)    
    [M1,M2]=size(x{t});
    X{t}=fft2(x{t});
    flipx{t}=rot90(x{t},2);
    Xr{t}=fft2(flipx{t});
    %b=b-2*conv2(flipx{t},y{t},'valid');
    b=b-2*valid_conv_by_fft(Xr{t},y{t});
    xx=x{t}.^2;
    for i=m1:-1:1
        for j=m2:-1:1
            Ad(m1+1-i,m2+1-j)=Ad(m1+1-i,m2+1-j)+sum(sum(xx(i:i+M1-m1,j:j+M2-m2)));
        end
    end
end
fprintf('Ad_min=%f,Ad_max=%f\n',min(Ad(:)),max(Ad(:)));
if lambda_C>0
    CtC = conv2(C,C,'same');
    Cd = CtC(2,2);
    Ad=Ad+lambda_C*Cd;
end
atAa=@(alpha)xtAx(alpha,X,x,L,lambda_C,C);
Aa=@(Xalpha)Ax(Xalpha,Xr,L,lambda_C,C);
end
function [y, Xalpha]=xtAx(alpha,X,x,L,lambda_C,C)
%y=alphat*XtX*alpha and Xalpha=X*alpha
y=0;
Xalpha=cell(2,1);
N=numel(alpha);
    for i=1:L
        if N>600
            Xalpha{i}=valid_conv_by_fft(X{i},alpha);%conv2(x{i},alpha,'valid');
        else
            Xalpha{i}=conv2(x{i},alpha,'valid');
        end
        y=y+sum(Xalpha{i}(:).^2);
    end
    if lambda_C>0
       Xalpha{L+1}=conv2(alpha,C,'same');
       y=y+lambda_C*norm(Xalpha{L+1},'fro')^2; 
    end
end
function [y]=Ax(Xalpha,Xr,L,lambda_C,C)
%XtX*alpha
y=0;
for i=1:L
    y=y+valid_conv_by_fft(Xr{i},Xalpha{i});%conv2(flipx{i},Xalpha{i},'valid');
end
if lambda_C>0
    y=y+lambda_C*conv2(Xalpha{L+1},C,'same');
end
end
function [y]=valid_conv_by_fft(X,h)
[M1,M2]=size(X);
[s1,s2]=size(h);
H=fft2(padarray(h,[M1-s1,M2-s2],0,'post'));
temp=real(ifft2(X.*H));
y=temp(s1:M1,s2:M2);
end

