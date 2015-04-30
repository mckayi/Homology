function z = EarthMoverFromIntervals(intervalsA,intervalsB,dimension,max_filtration)

%Load barcodes as Mx2 matrices
load(intervalsA);
StructA=HDintervals;
load(intervalsB);
StructB=HDintervals;
A=cell2mat(StructA(dimension+1));
B=cell2mat(StructB(dimension+1));
%Replace infinity with the max filtration time for both A and B
for i=1:size(A,1)
    if strcmp(num2str(A(i,2)),'Inf')
        A(i,2)=max_filtration;
    end
end
for i=1:size(B,1)
    if strcmp(num2str(B(i,2)),'Inf')
        B(i,2)=max_filtration;
    end
end

%Sort A and B
[Aunique,~,AcountIndex]=unique(A,'rows');
[Bunique,~,BcountIndex]=unique(B,'rows');
%Form histograms of points in A and B. The number of times each row is 
%repeated, will serve as the bin height.
%Vectorize everything (ie. use arrayfun)
Ahist = arrayfun(@(i)sum(ismember(A,Aunique(i,:),'rows')),1:size(Aunique,1))';
Bhist = arrayfun(@(i)sum(ismember(B,Bunique(i,:),'rows')),1:size(Bunique,1))';
%Form the projections proj_{[1;1]})[a;b]=1/2[a+b;a+b];
Aproj = cell2mat(arrayfun(@(i)1/2*[sum(Aunique(i,:)) sum(Aunique(i,:))],1:size(Aunique,1),'UniformOutput',false)');
Bproj = cell2mat(arrayfun(@(i)1/2*[sum(Bunique(i,:)) sum(Bunique(i,:))],1:size(Bunique,1),'UniformOutput',false)');
%Form full histograms (counts in set)\union (count of projections in other
%set)
Ahistfull = vertcat(Ahist,Bhist);
Bhistfull = vertcat(Bhist,Ahist);
Afullbasis = vertcat(Aunique,Bproj);
Bfullbasis = vertcat(Bunique,Aproj);
%Form the distance matrix D
D=zeros(length(Ahistfull),length(Bhistfull));
for i=1:length(Ahistfull)
    for j=1:length(Bhistfull)
        if i<=length(Ahist) && j<=length(Bhist) %both off diagonal
            D(i,j)=norm(Afullbasis(i,:)-Bfullbasis(j,:),2);
        elseif i<=length(Ahist) && j>length(Bhist) %one off diagonal
            D(i,j)=norm(Afullbasis(i,:)-Bfullbasis(j,:),2);
        elseif i>length(Ahist) && j<=length(Bhist) %one off diagonal
            D(i,j)=norm(Afullbasis(i,:)-Bfullbasis(j,:),2);
        else %both diagonal
            D(i,j)=0;
        end
    end
end
            
%Then do EMD
%this algorithm has a quirk in that it doesn't count the points that had to
%move the furthest distance. Hence the extra mass penalty is 0 (there should be
%no extra mass), and the flow return time is everything (3).
%Also, D is no longer a metric
[~,F] = emd_hat_mex(Ahistfull,Bhistfull,D,0,3);
%Form the cost matrix (flow and distance)
flowcost = F.*D;
z = sum(flowcost(:)); %This is the EMD
end