function res=woe_corr(res,i)



%for i=1:length(res.CATMoS_VT_pred)
M=zeros(5,7);
W=zeros(5,7);
AD=zeros(5,7);
ADi=zeros(5,7);
Cfi=zeros(5,7);

    if res.CATMoS_VT_pred(i,1)==1
        M(1,1:2)=ones(1,2);
        W(1,1:2)=repelem((res.AD_index_VT(i,1)+res.Conf_index_VT(i,1))*res.AD_VT(i,1)/2,2);
        AD(1,1:2)=repelem(res.AD_VT(i,1),2);
        ADi(1,1:2)=repelem(res.AD_index_VT(i,1),2);
        Cfi(1,1:2)=repelem(res.Conf_index_VT(i,1),2);
    else
        M(1,3:end)=ones(1,5);
        W(1,3:end)=repelem((res.AD_index_VT(i,1)+res.Conf_index_VT(i,1))*res.AD_VT(i,1)/2,5);
        AD(1,3:end)=repelem(res.AD_VT(i,1),5);
        ADi(1,3:end)=repelem(res.AD_index_VT(i,1),5);
        Cfi(1,3:end)=repelem(res.Conf_index_VT(i,1),5);
    end
    if res.CATMoS_NT_pred(i,1)==1
        M(2,6:7)=ones(1,2);
        W(2,6:7)=repelem((res.AD_index_NT(i,1)+res.Conf_index_NT(i,1))*res.AD_NT(i,1)/2,2);
        AD(2,6:7)=repelem(res.AD_NT(i,1),2);
        ADi(2,6:7)=repelem(res.AD_index_NT(i,1),2);
        Cfi(2,6:7)=repelem(res.Conf_index_NT(i,1),2);
    else
        M(2,1:5)=ones(1,5);
        W(2,1:5)=repelem((res.AD_index_NT(i,1)+res.Conf_index_NT(i,1))*res.AD_NT(i,1)/2,5);
        AD(2,1:5)=repelem(res.AD_NT(i,1),5);
        ADi(2,1:5)=repelem(res.AD_index_NT(i,1),5);
        Cfi(2,1:5)=repelem(res.Conf_index_NT(i,1),5);
    end
    switch(res.CATMoS_EPA_pred(i,1))
        case 1
            M(3,1:2)=ones(1,2);
            W(3,1:2)=repelem((res.AD_index_EPA(i,1)+res.Conf_index_EPA(i,1))*res.AD_EPA(i,1)/2,2);
            AD(3,1:2)=repelem(res.AD_EPA(i,1),2);
            ADi(3,1:2)=repelem(res.AD_index_EPA(i,1),2);
            Cfi(3,1:2)=repelem(res.Conf_index_EPA(i,1),2);
        case 2
            M(3,3:4)=ones(1,2);
            W(3,3:4)=repelem((res.AD_index_EPA(i,1)+res.Conf_index_EPA(i,1))*res.AD_EPA(i,1)/2,2);
            AD(3,3:4)=repelem(res.AD_EPA(i,1),2);
            ADi(3,3:4)=repelem(res.AD_index_EPA(i,1),2);
            Cfi(3,3:4)=repelem(res.Conf_index_EPA(i,1),2);
        case 3
            M(3,5:6)=ones(1,2);
            W(3,5:6)=repelem((res.AD_index_EPA(i,1)+res.Conf_index_EPA(i,1))*res.AD_EPA(i,1)/2,2);
            AD(3,5:6)=repelem(res.AD_EPA(i,1),2);
            ADi(3,5:6)=repelem(res.AD_index_EPA(i,1),2);
            Cfi(3,5:6)=repelem(res.Conf_index_EPA(i,1),2);
        case 4
            M(3,7)=1;
            W(3,7)=(res.AD_index_EPA(i,1)+res.Conf_index_EPA(i,1))*res.AD_EPA(i,1)/2;
            AD(3,7)=res.AD_EPA(i,1);
            ADi(3,7)=res.AD_index_EPA(i,1);
            Cfi(3,7)=res.Conf_index_EPA(i,1);
    end
    switch(res.CATMoS_GHS_pred(i,1))
        case 1
            M(4,1)=1;
            W(4,1)=(res.AD_index_GHS(i,1)+res.Conf_index_GHS(i,1))*res.AD_GHS(i,1)/2;
            AD(4,1)=res.AD_GHS(i,1);
            ADi(4,1)=res.AD_index_GHS(i,1);
            Cfi(4,1)=res.Conf_index_GHS(i,1);
        case 2
            M(4,2)=1;
            W(4,2)=(res.AD_index_GHS(i,1)+res.Conf_index_GHS(i,1))*res.AD_GHS(i,1)/2;
            AD(4,2)=res.AD_GHS(i,1);
            ADi(4,2)=res.AD_index_GHS(i,1);
            Cfi(4,2)=res.Conf_index_GHS(i,1);
        case 3
            M(4,3)=1;
            W(4,3)=(res.AD_index_GHS(i,1)+res.Conf_index_GHS(i,1))*res.AD_GHS(i,1)/2;
            AD(4,3)=res.AD_GHS(i,1);
            ADi(4,3)=res.AD_index_GHS(i,1);
            Cfi(4,3)=res.Conf_index_GHS(i,1);
        case 4
            M(4,4:5)=ones(1,2);
            W(4,4:5)=repelem((res.AD_index_GHS(i,1)+res.Conf_index_GHS(i,1))*res.AD_GHS(i,1)/2,2);
            AD(4,4:5)=repelem(res.AD_GHS(i,1),2);
            ADi(4,4:5)=repelem(res.AD_index_GHS(i,1),2);
            Cfi(4,4:5)=repelem(res.Conf_index_GHS(i,1),2);
        case 5
            M(4,6:7)=ones(1,2);
            W(4,6:7)=repelem((res.AD_index_GHS(i,1)+res.Conf_index_GHS(i,1))*res.AD_GHS(i,1)/2,2);
            AD(4,6:7)=repelem(res.AD_GHS(i,1),2);
            ADi(4,6:7)=repelem(res.AD_index_GHS(i,1),2);
            Cfi(4,6:7)=repelem(res.Conf_index_GHS(i,1),2);
    end
    if ~isnan(res.CATMoS_LD50_pred(i,1))&& (res.CATMoS_LD50_pred(i,1)+0.3)<log10(5000)
        A=find(((res.CATMoS_LD50_pred(i,1)-0.3)-log10([5,50,300,500,2000,5000]))<0);
        B=find(((res.CATMoS_LD50_pred(i,1)+0.3)-log10([5,50,300,500,2000,5000]))<0);
        M(5,A(1):B(1))=ones(1,B(1)-A(1)+1);
        W(5,A(1):B(1))=repelem((res.AD_index_LD50(i,1)+res.Conf_index_LD50(i,1))*res.AD_LD50(i,1)/2,B(1)-A(1)+1);
        AD(5,A(1):B(1))=repelem(res.AD_LD50(i,1),B(1)-A(1)+1);
        ADi(5,A(1):B(1))=repelem(res.AD_index_LD50(i,1),B(1)-A(1)+1);
        Cfi(5,A(1):B(1))=repelem(res.Conf_index_LD50(i,1),B(1)-A(1)+1);
    elseif ~isnan(res.CATMoS_LD50_pred(i,1))&& (res.CATMoS_LD50_pred(i,1)-0.3)<log10(5000)
        M(5,6:7)=ones(1,2);
        W(5,6:7)=repelem((res.AD_index_LD50(i,1)+res.Conf_index_LD50(i,1))*res.AD_LD50(i,1)/2,2);
        AD(5,6:7)=repelem(res.AD_LD50(i,1),2);
        ADi(5,6:7)=repelem(res.AD_index_LD50(i,1),2);
        Cfi(5,6:7)=repelem(res.Conf_index_LD50(i,1),2);
    elseif  ~isnan(res.CATMoS_LD50_pred(i,1))&& (res.CATMoS_LD50_pred(i,1)-0.3)>log10(5000)
        M(5,7)=1;
        W(5,7)=(res.AD_index_LD50(i,1)+res.Conf_index_LD50(i,1))*res.AD_LD50(i,1)/2;
        AD(5,7)=res.AD_LD50(i,1);
        ADi(5,7)=res.AD_index_LD50(i,1);
        Cfi(5,7)=res.Conf_index_LD50(i,1);
    end
    
    
    if (sum(M(:,1))>=3 && min(find(sum(M)==max(sum(M))))<=3 && ((res.CATMoS_GHS_pred(i,1)==1 && res.Conf_index_GHS(i,1)>=0.5)  || res.CATMoS_LD50_pred(i,1)<=log10(17.5))) || (sum(M(:,2))>=3 && min(find(sum(M)==max(sum(M))))<=3 && ((res.CATMoS_GHS_pred(i,1)==1 && res.Conf_index_GHS(i,1)>=0.5) || res.CATMoS_LD50_pred(i,1)<=log10(17.5)))
        WoE=1;
    elseif (sum(M(:,2))>=4 && min(find(sum(M)==max(sum(M))))<=3) || (sum(M(:,2))>=2 && min(find(sum(M)==max(sum(M))))<=3 && (res.CATMoS_GHS_pred(i,1)<=2 || res.CATMoS_LD50_pred(i,1)<=log10(75) || res.CATMoS_EPA_pred(i,1)==1 || (res.CATMoS_VT_pred(i,1)==1 && res.Conf_index_VT(i,1)>=0.55)))
        WoE=2;
    elseif sum(M(:,7))>=3 && res.CATMoS_LD50_pred(i,1)>=log10(4000) && res.Conf_index_LD50(i,1)>= res.Conf_index_EPA(i,1)
        WoE=7;
    elseif (sum(M(:,6))>=4 && max(find(sum(M)==max(sum(M))))>=6) && res.CATMoS_NT_pred(i,1)==1 && res.CATMoS_LD50_pred(i,1)>=log10(2500) && res.CATMoS_LD50_pred(i,1)<=log10(5000) && res.CATMoS_EPA_pred(i,1)==3
        WoE=6;
    elseif (res.AD_index_VT(i,1)+ res.AD_index_NT(i,1) + res.AD_index_EPA(i,1) + res.AD_index_GHS(i,1)+ res.AD_index_LD50(i,1)<5) && res.CATMoS_LD50_pred(i,1)<=log10(500) && min(find(sum(M)>=2))<=4
        %WoE=min(find(sum(M)>=2));
        WoE=find([res.CATMoS_LD50_pred(i,1)<=log10(5) res.CATMoS_LD50_pred(i,1)>log10(5)&&res.CATMoS_LD50_pred(i,1)<=log10(50) res.CATMoS_LD50_pred(i,1)>log10(50)&&res.CATMoS_LD50_pred(i,1)<=log10(300)...
            res.CATMoS_LD50_pred(i,1)>log10(300)&&res.CATMoS_LD50_pred(i,1)<=log10(500) res.CATMoS_LD50_pred(i,1)>log10(500)&&res.CATMoS_LD50_pred(i,1)<=log10(2000) res.CATMoS_LD50_pred(i,1)>log10(2000)&&res.CATMoS_LD50_pred(i,1)<=log10(5000)]);
        if WoE > min(find(sum(M)>=2))
            WoE=min(find(sum(M)>=2));
        end
    else
        WoE=find(sum(M)==max(sum(M)));
        if length(WoE)>1
            if min(WoE)<=3 && max(WoE)<=4 && res.CATMoS_LD50_pred(i,1)<=log10(500)%&& abs(sum(M(max(WoE)))-sum(M(min(WoE))))==1
                if min(WoE)>1
                    WoE=min(WoE);
                elseif res.CATMoS_GHS_pred(i,1)==1 || res.CATMoS_LD50_pred(i,1)<=log10(100)
                    WoE=min(WoE);
                elseif res.CATMoS_LD50_pred(i,1)<=log10(300)
                    WoE=round(mean(WoE));
                else
                    WoE=max(WoE);
                end

            elseif max(WoE)>=6 && (res.CATMoS_NT_pred(i,1)==1 || res.CATMoS_LD50_pred(i,1)>=log10(2500))
                WoE=max(WoE);
            else
                W=sum(W)./sum(M);
                if max(W(WoE))>=max(W(W(WoE)~=max(W(WoE))))+max(W(W(WoE)~=max(W(WoE))))/3
                    WoE=WoE(W(WoE)==max(W(WoE)));
                end
                WoE=min(WoE);
                %WoE
            end
        else
            %M=sum(M);
            %if min(find(M==max(M(M~=max(M)))))==(WoE-1) && WoE<=4 && WoE>3 && M(min(find(M==max(M(M~=max(M))))))==max(M)-1 && res.CATMoS_LD50_pred(i,1)<=log10(450)
            if WoE==4 && sum(M(:,4))-sum(M(:,3))<=1 && res.CATMoS_LD50_pred(i,1)<=log10(500)
                WoE=WoE-1;
            end
        end
    end
    WoE=min(WoE);
    %WoE
    switch WoE
        case 1
            res.CATMoS_VT_pred(i,1)=1;
            res.CATMoS_NT_pred(i,1)=0;
            res.CATMoS_EPA_pred(i,1)=1;
            res.CATMoS_GHS_pred(i,1)=1;
            
%             if res.CATMoS_LD50_pred(i,1)-0.3<= log10(5)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(5))/2;
%             elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(5)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
%             elseif res.CATMoS_LD50_pred(i,1)> log10(5)
%                 res.CATMoS_LD50_pred(i,1)=log10(2.5);
%             end
            
            if res.CATMoS_LD50_pred(i,1)> log10(5)
                if res.CATMoS_LD50_pred(i,1)-0.3<= log10(5)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(5))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-1)/0.65<= log10(5)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-1)/0.65;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(5)*res.CATMoS_LD50_pred(i,1)/4;
                end
            end



        case 2
            res.CATMoS_VT_pred(i,1)=1;
            res.CATMoS_NT_pred(i,1)=0;
            res.CATMoS_EPA_pred(i,1)=1;
            res.CATMoS_GHS_pred(i,1)=2;
            
%             if res.CATMoS_LD50_pred(i,1)-0.3<= log10(50) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(50))/2>= log10(5)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(50))/2;
%             elseif res.CATMoS_LD50_pred(i,1)+0.3>= log10(5) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(5))/2<= log10(50)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(5))/2;
%             elseif(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(50) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(5)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
%             elseif res.CATMoS_LD50_pred(i,1)> log10(50) || res.CATMoS_LD50_pred(i,1)< log10(5)
%                 res.CATMoS_LD50_pred(i,1)=log10(27.5); 
%             end
            
            if res.CATMoS_LD50_pred(i,1)> log10(50)
                if res.CATMoS_LD50_pred(i,1)-0.3<= log10(50) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(50))/2>= log10(5)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(50))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-1)/0.65<= log10(50) && (res.CATMoS_LD50_pred(i,1)-1)/0.65>= log10(5)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-1)/0.65;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(5)+(log10(50)-log10(5))*res.CATMoS_LD50_pred(i,1)/4;
                end
            elseif res.CATMoS_LD50_pred(i,1)< log10(5)
                if res.CATMoS_LD50_pred(i,1)+0.3>= log10(5) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(5))/2<= log10(50)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(5))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(5) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(50)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(5)+(log10(50)-log10(5))*res.CATMoS_LD50_pred(i,1)/4;
                end
            end
        case 3
            res.CATMoS_VT_pred(i,1)=0;
            res.CATMoS_NT_pred(i,1)=0;
            res.CATMoS_EPA_pred(i,1)=2;
            res.CATMoS_GHS_pred(i,1)=3;
            
%             if res.CATMoS_LD50_pred(i,1)-0.3<= log10(300) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(300))/2>= log10(50)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(300))/2;
%             elseif res.CATMoS_LD50_pred(i,1)+0.3>= log10(50) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(50))/2<= log10(300)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(50))/2;
%             elseif(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(300) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(50)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
%             elseif res.CATMoS_LD50_pred(i,1)> log10(300) || res.CATMoS_LD50_pred(i,1)< log10(50)
%                 res.CATMoS_LD50_pred(i,1)=log10(175); 
%             end
            
            if res.CATMoS_LD50_pred(i,1)> log10(300)
                if res.CATMoS_LD50_pred(i,1)-0.3<= log10(300) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(300))/2>= log10(50)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(300))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(300) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(50)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(50)+(log10(300)-log10(50))*res.CATMoS_LD50_pred(i,1)/4;
                end
            elseif res.CATMoS_LD50_pred(i,1)< log10(50)
                if res.CATMoS_LD50_pred(i,1)+0.3>= log10(50) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(50))/2<= log10(300)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(50))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(50) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(300)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(50)+(log10(300)-log10(50))*res.CATMoS_LD50_pred(i,1)/4;
                end
            end
        case 4
            res.CATMoS_VT_pred(i,1)=0;
            res.CATMoS_NT_pred(i,1)=0;
            res.CATMoS_EPA_pred(i,1)=2;
            res.CATMoS_GHS_pred(i,1)=4;
            
%             if res.CATMoS_LD50_pred(i,1)-0.3<= log10(500) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(500))/2>= log10(300)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(500))/2;
%             elseif res.CATMoS_LD50_pred(i,1)+0.3>= log10(300) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(300))/2<= log10(500)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(300))/2;
%             elseif(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(500) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(300)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
%             elseif res.CATMoS_LD50_pred(i,1)> log10(500) || res.CATMoS_LD50_pred(i,1)< log10(300)
%                 res.CATMoS_LD50_pred(i,1)=log10(400); 
%             end
            
            if res.CATMoS_LD50_pred(i,1)> log10(500)
                if res.CATMoS_LD50_pred(i,1)-0.3<= log10(500) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(500))/2>= log10(300)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(500))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(500) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(300)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(300)+(log10(500)-log10(300))*res.CATMoS_LD50_pred(i,1)/4;
                end
            elseif res.CATMoS_LD50_pred(i,1)< log10(300)
                if res.CATMoS_LD50_pred(i,1)+0.3>= log10(300) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(300))/2<= log10(500)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(300))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(300) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(500)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(300)+(log10(500)-log10(300))*res.CATMoS_LD50_pred(i,1)/4;
                end
            end
        case 5
            res.CATMoS_VT_pred(i,1)=0;
            res.CATMoS_NT_pred(i,1)=0;
            res.CATMoS_EPA_pred(i,1)=3;
            res.CATMoS_GHS_pred(i,1)=4;
            
%             if res.CATMoS_LD50_pred(i,1)+0.3>= log10(500) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(500))/2<= log10(2000)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(500))/2;
%                 elseif res.CATMoS_LD50_pred(i,1)-0.3<= log10(2000) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(2000))/2>= log10(500)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(2000))/2;
%             elseif(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(2000) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(500)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
%             elseif res.CATMoS_LD50_pred(i,1)> log10(2000) || res.CATMoS_LD50_pred(i,1)< log10(500)
%                 res.CATMoS_LD50_pred(i,1)=log10(1250); 
%             end
            
            if res.CATMoS_LD50_pred(i,1)> log10(2000)
                if res.CATMoS_LD50_pred(i,1)-0.3<= log10(2000) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(2000))/2>= log10(500)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(2000))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(2000) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(500)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(500)+(log10(2000)-log10(500))*res.CATMoS_LD50_pred(i,1)/4;
                end
            elseif res.CATMoS_LD50_pred(i,1)< log10(500)
                if res.CATMoS_LD50_pred(i,1)+0.3>= log10(500) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(500))/2<= log10(2000)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(500))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(500) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(2000)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(500)+(log10(2000)-log10(500))*res.CATMoS_LD50_pred(i,1)/4;
                end
            end
        case 6
            res.CATMoS_VT_pred(i,1)=0;
            res.CATMoS_NT_pred(i,1)=1;
            res.CATMoS_EPA_pred(i,1)=3;
            res.CATMoS_GHS_pred(i,1)=5;
            
%             if res.CATMoS_LD50_pred(i,1)+0.3>= log10(2000) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(2000))/2<= log10(5000)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(2000))/2;
%             elseif(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(5000) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(2000)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
%             elseif res.CATMoS_LD50_pred(i,1)-0.3<= log10(5000) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(5000))/2>= log10(2000)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(5000))/2;
%             elseif res.CATMoS_LD50_pred(i,1)> log10(5000) || res.CATMoS_LD50_pred(i,1)< log10(2000)
%                 res.CATMoS_LD50_pred(i,1)=log10(3500); 
%             end
            
            if res.CATMoS_LD50_pred(i,1)> log10(5000) && (res.CATMoS_LD50_pred(i,1)-0.3+log10(5000))/2>= log10(2000)
                if res.CATMoS_LD50_pred(i,1)-0.3<= log10(5000)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.3+log10(5000))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(5000) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(2000)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(2000)+(log10(5000)-log10(2000))*res.CATMoS_LD50_pred(i,1)/4;
                end
            elseif res.CATMoS_LD50_pred(i,1)< log10(2000)
                if res.CATMoS_LD50_pred(i,1)+0.3>= log10(2000) && (res.CATMoS_LD50_pred(i,1)+0.3+log10(2000))/2<= log10(5000)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(2000))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(2000) && (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369<= log10(5000)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(2000)+(log10(5000)-log10(2000))*res.CATMoS_LD50_pred(i,1)/4;
                end
                
            end
        case 7
            res.CATMoS_VT_pred(i,1)=0;
            res.CATMoS_NT_pred(i,1)=1;
            res.CATMoS_EPA_pred(i,1)=4;
            res.CATMoS_GHS_pred(i,1)=5;
            
%             if res.CATMoS_LD50_pred(i,1)+0.3>= log10(5000)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(5000))/2;
%             elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(5000)
%                 res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
%             elseif res.CATMoS_LD50_pred(i,1)< log10(5000)
%                 res.CATMoS_LD50_pred(i,1)=log10(7500);
%             end
            
            if res.CATMoS_LD50_pred(i,1)< log10(5000)
                if res.CATMoS_LD50_pred(i,1)+0.3>= log10(5000)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)+0.3+log10(5000))/2;
                elseif (res.CATMoS_LD50_pred(i,1)-0.781)/0.7369>= log10(5000)
                    res.CATMoS_LD50_pred(i,1)=(res.CATMoS_LD50_pred(i,1)-0.781)/0.7369;
                else
                    res.CATMoS_LD50_pred(i,1)=log10(5000)+(log10(10000)-log10(5000))*res.CATMoS_LD50_pred(i,1)/4;
                end
            end
    end
    
    res.AD_CATMoS(i,1)=round(mean(AD(find(AD(:,WoE)),WoE)));
    res.AD_index_CATMoS(i,1)=max(mean(ADi(find(ADi(:,WoE)),WoE)),median(ADi(find(ADi(:,WoE)),WoE)));
    %res.Conf_index_CATMoS(i,1)=max(mean(Cfi(find(Cfi(:,WoE)),WoE)),median(Cfi(find(Cfi(:,WoE)),WoE)));
    res.Conf_index_CATMoS(i,1)=median([mean(Cfi(find(Cfi(:,WoE)),WoE)),median(Cfi(find(Cfi(:,WoE)),WoE)),sum(M(:,WoE))/5]);
    
   
%end


% res.AD_VT=[];
% res.AD_index_VT=[];
% res.Conf_index_VT=[];
% res.AD_NT=[];
% res.AD_index_NT=[];
% res.Conf_index_NT=[];
% res.AD_EPA=[];
% res.AD_index_EPA=[];
% res.Conf_index_EPA=[];
% res.AD_GHS=[];
% res.AD_index_GHS=[];
% res.Conf_index_GHS=[];
% res.AD_LD50=[];
% res.AD_index_LD50=[];
% res.Conf_index_LD50=[];


    
    
