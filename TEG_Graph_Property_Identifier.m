%% Identifying Properties 
function TEG_Graph_parameters = TEG_Graph_Property_Identifier(T,Y_input,type)
Y_TEG=Y_input/2;

% R-time : Time to start forming clot, reach 2mm (5-10 minutes)
data_id=find(abs(Y_TEG-2)<.2,1);
if isempty(data_id)==1
    data_id=find(abs(Y_TEG-2)<.5,1);
end
tegR_id=data_id;
TEG_R=T(data_id);

% K time : time to reach a certain strength, after R time, from 2 to 20 mm amplitude (1-3 minutes)
data_id=find(abs(Y_TEG-20)<1,1);
if isempty(data_id)==1
    data_id=find(abs(Y_TEG-20)<2,1);
end
time_20mm=T(data_id);
TEG_K=time_20mm-TEG_R;

% Alpha angle : speed of fibrin accumulation (53-72 degrees)
dydx = gradient(Y_TEG(:)) ./ gradient(T(:));
dydx2 = gradient(dydx(:)) ./ gradient(T(:)); 
smDyDx=smoothdata(smoothdata(dydx2));
aa=find(abs(smDyDx)<0.5);
[~,bb]=min(smDyDx);
data_id=aa(find(aa>bb,1));
TEG_alpha=abs(atand(dydx(data_id)));

% Maximum amplitude : Highest vertical amplitude of TEG (50-70 mm) - clot strength 
if type==1
%Method 1
    [TEG_MA1,data_id]=max(Y_TEG);
    TEG_MAtime1= T(data_id) ;
    tegMA_id1=data_id;
    TEG_MA=TEG_MA1;
    TEG_MAtime=TEG_MAtime1;
    tegMA_id=tegMA_id1;

elseif type==2
%Method 2
    aa=find(abs(dydx)<0.02);
    [~,bb]=max(dydx);
    data_id=aa(find(aa>bb,1));
    TEG_MA2=Y_TEG(data_id);
    TEG_MAtime2=T(data_id);
    tegMA_id2=data_id;
    TEG_MA=TEG_MA2;
    TEG_MAtime=TEG_MAtime2;
    tegMA_id=tegMA_id2;
end

%Degredation (AUC) - Method 2
if tegMA_id+360<=length(T)
    AUC30=trapz(T(tegMA_id:tegMA_id+360),Y_TEG(tegMA_id:tegMA_id+360));
else
    c1=polyfit(T(tegMA_id:end),Y_TEG(tegMA_id:end),1);
    y1 = polyval(c1,T(tegMA_id)+30);
    AUC30=trapz([T(tegMA_id:end);T(tegMA_id)+30],[Y_TEG(tegMA_id:end);y1]);
end
AUC_MA=30*TEG_MA;
Ly30_method2=(AUC_MA-AUC30)/AUC_MA*100;
TEG_Ly30=Ly30_method2;

TEG_Graph_parameters=[TEG_R,TEG_K,TEG_alpha,TEG_MA,TEG_Ly30,TEG_MAtime];
return