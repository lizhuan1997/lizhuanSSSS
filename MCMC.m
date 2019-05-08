%%%MCMC buckyball
%%%2019.5.6
beta=1.5;
J=1;
global B

eMC=[];
E=[];
m=20000;
sample1=[];
for beta=1
B=[(J*beta),(-J*beta);(-J*beta),(J*beta)];
samples=rand(60,m);
loops=60;
samples(samples>0.5)=1;
samples(samples<0.5)=2;


indx=[5,6,2;1,7,3;2,8,4;3,9,5;4,10,1;...
      1,20,11;2,12,13;3,14,15;4,16,17;5,18,19;...
      6,21,12;7,11,22;7,23,14;8,13,24;8,25,16;9,15,26;9,27,18;10,17,28;10,29,20;6,19,30;...
      11,30,31;12,32,23;13,22,33;14,34,25;15,24,35;16,36,27;17,26,37;18,38,29;19,28,39;20,40,21;...
      21,41,32;22,31,42;23,34,43;24,33,44;25,45,36;26,46,35;27,47,38;28,48,37;29,49,40;30,50,39;...
      31,50,51;32,51,43;33,42,52;34,52,45;35,44,53;36,53,47;37,46,54;38,54,49;39,48,55;40,55,41;...
      41,42,56;43,44,57;45,58,46;47,48,59;49,50,60;...
      60,51,57;56,52,58;57,53,59;54,58,60;55,59,56];
e=energy(samples,indx);

  for loop=1:loops
       for i =1:60

                j1=samples(indx(i,1),:);
                j2=samples(indx(i,2),:);
                j3=samples(indx(i,3),:);


                in1=(j1-1)*2+samples(i,:);
                in2=(j2-1)*2+samples(i,:);
                in3=(j3-1)*2+samples(i,:);

                e1=B(in1)+B(in2)+B(in3);
                e2=B(in1-samples(i,:)+3-samples(i,:))+B(in2-samples(i,:)-samples(i,:)+3)+B(in3-samples(i,:)-samples(i,:)+3);
                
                p=exp(beta*(e1-e2));
                x=rand(1,m);
                
                
                samples2=3-samples(i,:);
                ek=e+e2-e1;
                
                e(x<p)=ek(x<p);
                samples(i,x<p)=samples2(x<p);
                
                e(e2<e1)=ek(e2<e1);
                samples(i,e2<e1)=samples2(e2<e1);
                if mod(i,5)==0
                    E=[E,e];
                    sample1=[sample1,samples];
                end
                
       end
       
  end
  eMC=[eMC,mean(energy(samples,indx))];
%   
  
end

[~,b]=unique(sample1','row');
n=sum(e(b)==-66);
disp("# of degeneration = "+n);


  function e=energy(samples,indx)
  global B
        e=zeros(1,size(samples,2));
        for i =1:60
            j1=samples(indx(i,1),:);
            j2=samples(indx(i,2),:);
            j3=samples(indx(i,3),:);
            
            
            in1=(j1-1)*2+samples(i,:);
            in2=(j2-1)*2+samples(i,:);
            in3=(j3-1)*2+samples(i,:);
%             
            
            e=e+B(in1);
            e=e+B(in2);
            e=e+B(in3);
            
        end
        
        
        e=e/2;
  end
  
  

