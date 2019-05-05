function OutImage=mirror(InImage,n)
%´¹Ö±¾µÏñ±ä»»
I=InImage;
[M,N,G]=size(I);
J=I;
for i=1:M
    for j=1:N
        J(i,j,:)=I(i,N-j+1,:);
    end
end
OutImage=J;