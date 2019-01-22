function result=viterbi(S,hmm,ep)
    
    % conversao de probabilidades em logaritmo
    hmm=log2(hmm);
    ep=log2(ep);
    
    [m,~] = size(hmm);

    s=zeros(1,length(S));
    
    %descodificação do vetor de letras
    for i=1:length(S)
        if(S(i)=='G')
            s(i)=3;
        end
        if(S(i)=='A')
            s(i)=1;
        end
        if(S(i)=='C')
            s(i)=2;
        end
        if(S(i)=='T')
            s(i)=4;
        end
    end
    
    res=zeros(1,length(S));
    score=zeros(m,length(S));
    scorei=zeros(1,m);
    
    %inicialização
    for j=1:m
        score(j,1)=log2(1/m)+ep(s(1),j);
    end
    
    %calculo dos diferentes scores
    for i=2:length(s)
        for j=1:m
            for k=1:m
                scorei(k)=score(k,i-1)+hmm(k,j);
            end
            [~,ptr(j,i)]=max(scorei);
            score(j,i)=ep(s(i),j)+max(scorei);
        end
    end
    
    [~,result(length(s))]=max(score(:,length(s)));
    
    for i=length(s):-1:2
       result(i-1)=ptr(result(i),i); 
    end
       
    disp('Sequência:');
    disp(result);
end