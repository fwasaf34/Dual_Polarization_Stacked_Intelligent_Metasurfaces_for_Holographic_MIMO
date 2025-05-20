

cap=0;
capWF=0;
for s1=1:2*S
    inte=0;
    inteWF=0;
    for s2=1:2*S
        if s2~=s1
            inte=inte+P/(2*S)*abs(alpha*DpH(s1,s2)).^2;
            inteWF=inteWF+Ps(s2)*abs(alpha*DpH(s1,s2)).^2;
        end
    end
    cap=cap+log2(1+P/(2*S)*abs(alpha*DpH(s1,s1)).^2/(inte+noise));
    capWF=capWF+log2(1+Ps(s1)*abs(alpha*DpH(s1,s1)).^2/(inteWF+noise));
end

CDpHmimo(round)=cap;
CDpHmimoWF(round)=capWF;



% alpha*DpH-lam
% norm(alpha*DpH-lam,'fro')^2/norm(lam,'fro')^2