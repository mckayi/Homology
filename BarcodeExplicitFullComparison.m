function [y,filter0,filter1,filtertime]=BarcodeExplicitFullComparison(intervals,bcodes)
%import edu.stanford.math.plex4.*; %Don't actually need this
load(intervals);
load(bcodes);
dim0=cell2mat(HDintervals(1));
dim1=cell2mat(HDintervals(2));
y='Failed';
%checkmate = zeros(3);
checkmate = zeros(1,3); %Only need a 1 by three vector
z=0;
%filter=15000;
filter=Inf; %An arbitrarily large value to start out with
filtertime=0;
error=-1;
%These couple of for loops cause a rather long execution time. A quicker
%way to select all those intervals that start & end when you want it to is
%to do something like find((dim1(:,1) > .4) & (dim1(:,2) < .8)) This is all
%the intervals in dim1 that start after .4 and end before .8
for timestep=.04:.08:8
    flag = 0;
    flagZero=0;
   for j=1:size(dim0,1)
       if dim0(j,2)<timestep
           continue
       end
       if strcmp(num2str(dim0(j,2)),'Inf')
           dim0(j,2)=8;
       end
       if dim0(j,1)<timestep
          flagZero=flagZero+1; 
       end
   end
   for j=1:size(dim1,1)
       if dim1(j,2)<timestep
           continue
       end
       if strcmp(num2str(dim1(j,2)),'Inf')
           dim1(j,2)=8;
       end
       if dim1(j,1)<timestep
          flag=flag+1; 
       end
   end
   z=z+1;
   checkmate(z,1)=abs(flag-barcodes(2,1));
   checkmate(z,2)=abs(flagZero-barcodes(1,1));
   checkmate(z,3)=timestep;
   %error=sqrt(checkmate(z,1)^2 + checkmate(z,2)^2);
   error=checkmate(z,1)+checkmate(z,2); %just take the total distance
   if error<filter
       filter=error;
       filter0=checkmate(z,2);
       filter1=checkmate(z,1);
       filtertime=timestep;
   end
   if filter==0
       y='Succeeded';
       break
   end
end
end


    