function [midEN, SNR, R]=thre_search(S, pun, iter)
    [m, n]=size(S);
    R =(n-m)/(n-length(pun)); % code rate calculating from protograph size

    leEN=-10;
    riEN=10;
    while( ~pexit(S,riEN,R,pun, iter))
        leEN=riEN;
        riEN=riEN+1;
    end
    for L=1:1:20
        midEN=(riEN+leEN)/2;
        if pexit(S,midEN,R,pun, iter)
            riEN=midEN;
        else
            leEN=midEN;
        end
    end
    midEN =(riEN+leEN)/2;
    SNR = midEN + 10*log10(log2(4)*R);
end