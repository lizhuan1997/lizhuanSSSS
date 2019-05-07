%%%contraction2
%%%2019.5.7

F=free_energy(1);
disp("F=" + F)
n=free_energy(10)*60-10*66;
disp("number of degeneration = " +round(exp(n)))


function F=free_energy(beta)
J=1;
B=[exp(-J*beta),exp(J*beta);exp(J*beta),exp(-J*beta)];
k=norm(B,'fro');
% k=1;
B=B/k;
id=zeros(60,3);
indx=[5,6,2;1,7,3;2,8,4;3,9,5;4,10,1;...
      1,20,11;2,12,13;3,14,15;4,16,17;5,18,19;...
      6,21,12;7,11,22;7,23,14;8,13,24;8,25,16;9,15,26;9,27,18;10,17,28;10,29,20;6,19,30;...
      11,30,31;12,32,23;13,22,33;14,34,25;15,24,35;16,36,27;17,26,37;18,38,29;19,28,39;20,40,21;...
      21,41,32;22,31,42;23,34,43;24,33,44;25,45,36;26,46,35;27,47,38;28,48,37;29,49,40;30,50,39;...
      31,50,51;32,51,43;33,42,52;34,52,45;35,44,53;36,53,47;37,46,54;38,54,49;39,48,55;40,55,41;...
      41,42,56;43,44,57;45,58,46;47,48,59;49,50,60;...
      60,51,57;56,52,58;57,53,59;54,58,60;55,59,56];
  
for i =1:60
    id(i,:)=[min(i,indx(i,1))*100+max(i,indx(i,1)),min(i,indx(i,2))*100+max(i,indx(i,2)),min(i,indx(i,3))*100+max(i,indx(i,3))];
    
end


E=zeros(2,2,2);
E(1,1,1)=1;
E(2,2,2)=1;
t1=E;
idt1=id(1,:);
for i =[6,11,12,21,22,31,32,41,42,51,56]
    [t1,idt1]=tri_tensor_product(t1,idt1,E,id(i,:),B);
end

MT=permute(t1,[1,3,5,7,9,2,4,6,8,10]);
for i =1:5
    MT=tensor_product(MT,[1,2,3,4,5,6,7,8,9,10],B,[6,11]);
    
end
MT=reshape(MT,32,32);
F=(log(trace(MT^5))+log(k)*90)/60;


end
function [t,idt]=tri_tensor_product(T1,idT1,T2,idT2,B)
        id1=idT1;
for i =1:length(id1)
    if sum(idT2==id1(i))>0
         idT2(idT2==id1(i))=id1(i)*10;
        [T1,idT1]=tensor_product(T1,idT1,B,[id1(i),id1(i)*10]);
       
    end
    
    
end
[t,idt]=tensor_product(T1,idT1,T2,idT2);
end