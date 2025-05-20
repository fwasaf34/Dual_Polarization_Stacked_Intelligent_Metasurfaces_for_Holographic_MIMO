% 
% % 
% cap=0;
% for s1=1:S
%     for s2=1:S
%         if s2~=s1
%             inte(s2)=P/S*abs(alpha*H(s1,s2)).^2;
%         end
%     end
%     cap=cap+log2(1+P/S*abs(alpha*H(s1,s1)).^2/(sum(inte)+noise));
% end
% 
% CHmimo(round)=cap;


cap=0;
capWF=0;
for s1=1:S
    inte=0;
    inteWF=0;
    for s2=1:S
        if s2~=s1
            inte=inte+P/S*abs(alpha*H(s1,s2)).^2;
            inteWF=inteWF+Ps(s2)*abs(alpha*H(s1,s2)).^2;
        end
    end
    cap=cap+log2(1+P/S*abs(alpha*H(s1,s1)).^2/(inte+noise));
    capWF=capWF+log2(1+Ps(s1)*abs(alpha*H(s1,s1)).^2/(inteWF+noise));
end

CHmimo(round)=cap;
CHmimoWF(round)=capWF;