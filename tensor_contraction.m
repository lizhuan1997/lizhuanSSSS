%%%tensor contraction
%%%2019.5.6
global J
J=1;

disp("F= "+free_energy(beta))
n=free_energy(10)*60-10*66;

disp("number of degeneration = " +round(exp(n)))


function F=free_energy(beta)
global J
tensor1=zeros(2,2,2);
tensor1(1,1,1)=1;
tensor1(2,2,2)=1;
B=[exp(-J*beta),exp(J*beta);exp(J*beta),exp(-J*beta)];


k=norm(B,'fro');
% k=1;
B=B/k;

tensor2=tensor_product(tensor1,[1,2,3],B,[3,4]);
tensor2=tensor_product(tensor2,[1,2,3],tensor1,[3,4,5]);


t1=tensor_product(tensor1,[1,2,3],B,[3,4]);
t1=permute(t1,[3,1,2]);

t2=tensor_product(tensor2,[1,2,3,4],B,[1,5]);
t2=tensor_product(t2,[1,2,3,4],B,[2,5]);
t2=permute(t2,[4,2,1,3]);

t3=tensor_product(t2,[1,2,3,4],B,[2,5]);
t3=permute(t3,[1,4,2,3]);

T1=tensor_product(t1,[1,2,3],t2,[4,5,6,2]);
T1=permute(T1,[1,3,4,2,5]);

T2=tensor_product(T1,[1,2,3,4,5],t3,[6,7,8,3]);
T2=permute(T2,[1,2,5,6,3,4,7]);

t22=tensor_product(tensor2,[1,2,3,4],B,[1,5]);
t22=tensor_product(t22,[1,2,3,4],B,[1,5]);
t22=permute(t22,[3,4,1,2]);


T12=tensor_product(t1,[1,2,3],t22,[4,5,6,2]);
T12=permute(T12,[1,3,4,2,5]);


T2=tensor_product(T2,[1,2,3,4,5,6,7],T12,[8,9,4,10,11]);
T2=permute(T2,[1,2,3,8,7,4,5,6,10,9]);

t31=permute(t3,[1,4,3,2]);
T11=tensor_product(t1,[1,2,3],t31,[4,5,6,2]);
T11=permute(T11,[1,3,4,2,5]);

T3=tensor_product(T11,[1,2,3,4,5],t1,[11,3,13]);
T3=permute(T3,[1,2,5,3,4,6]);

MT=tensor_product(T2,[1,2,3,4,5,6,7,8,9,10],T3,[7,8,9,11,12,13]);
MT=permute(MT,[1,2,3,4,5,6,8,9,10,7]);




MT=reshape(MT,32,32);
logZ=norm(MT);
MT=MT/logZ;
Z=(trace(MT^5));


logz=log(logZ)*5+log(Z)+log(k)*90;
F=logz/60;

end



% exp(log(z(100))-(betan*E1))

