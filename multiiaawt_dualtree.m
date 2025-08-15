function [sgates,errors]=multiiaawt_dualtree(data_array,numsurrogates,primary,accerror,error_change)

%This is a multivariate time-series variant of the IAAWT algorithm as presented in:

%Keylock, C.J. 2017. Multifractal surrogate-data generation algorithm that preserves pointwise 
%Hölder regularity structure, with initial applications to turbulence, Physical Review E 95, 032123, 
%https://doi.org/10.1103/PhysRevE.95.032123.

%If you use this code or a derivative thereof for your own work, I would be
%grateful if you could acknowledge this source.

%The original algorithm made use of the dual tree complex wavelet transform
%developed by Nick Kingsbury in Cambridge and Nick's implementation thereof.
%This variant uses the dualtree function, which is part of the wavelet 
% toolbox in more recent variants of Matlab. 

%INPUTS

%data_array (length(data_array) by 1 data vector)
%This is a dyadic wavelet transform. If you pass it data that is not of
%length 2^numlevels where numlevels is an integer given by
%log(length(data_array)) / log(2) then it will zero pad up to the next
%dyadic scale.

%numsurrogates (the number of surrogates you wish to generate)

%primary is an integer specifying the column from which the phase
%differences are calculated. If not specified, this is quickly estimated from the
%difference statistics for each column

%accerror (the acceptable error)
%Both the error from the wavelet part and from the restoration of the
%original values must be below this value for the code to terminate

%error_change (acceptable relative error)
%If the change in accerror from one iteration to the next is less than
%accerror / error_change then the code terminates

%OUTPUTS

%sgates (length(data_array) by numsurrogtes array)
%The generated surrogate data

%errors (numsurrogates by 1 vector)
%The value for toterror at the termination of the algorithm for this
%surrogate.

%  DISCLAIMER

%  We make no warranties, explicit or implicit, that this program
%  is free of error, or is consistent with any particular standard 
%  of accuracy, or that it will meet your requirements for any 
%  particular application.

%  It should not be relied on for any purpose where incorrect
%  results could result in loss of property or personal injury.
%  If you do use this program for any such purpose it is at your own
%  risk.  The author disclaims all liability of any kind, either
%  direct or consequential, resulting from your use of this program.

%It is assumed there are more values in the time-series than there are
%time-series. Hence data_array is n by m with n>m
sizer=size(data_array);
if sizer(1)<sizer(2)
    data_array=data_array';
    sizer=size(data_array);   
end
exactlevels=log(length(data_array))/log(2);
numlevels=(floor(log(length(data_array))/log(2)));

if nargin<5
    error_change=100;
end
if nargin<4
    accerror=.001;
end
if nargin<3
    if sizer(2)<=2
        primary=1;
    else
        for loop1=1:sizer(2)
            temp(:,1)=data_array(:,loop1);
            temp(:,2)=circshift(temp(:,1),[-1 0]);
            res(loop1,1)=mean(abs(temp(1:length(temp)-1,1)-temp(1:length(temp)-1,2)))./std(temp(:,1));
        end
        [~,primary]=max(res);
    end
end
primary=round(primary);
if primary<1 || primary>sizer(2)
    primary=1;
end
if nargin<2
    numsurrogates=1;
end


count=1;
if abs(numlevels-exactlevels)>eps
    %We need to zero pad
    temp=zeros(2^(numlevels+1),sizer(2));
    count=2^(numlevels+1)-length(data_array)+1;
    temp(count:length(temp),:)=data_array;
    data_array=temp;
end

not_primary=setdiff(1:sizer(2),primary);

%wavelet decomp of primary
[Yl{primary},Yh] = dualtree(data_array(:,primary),'Level',numlevels,'LevelOneFilter','nearsym5_7','FilterLength',14);
% angs and amplitudes
for loop1=1:numlevels;
    ampYh{loop1}=abs(Yh{loop1});
    phaseYh{loop1}=angle(Yh{loop1});
end
%wavelet decomp of non-primary and wavelet phase diffs stored
for loop2=1:length(not_primary)
    [Yl{not_primary(loop2)},Yh] = dualtree(data_array(:,not_primary(loop2)),'Level',numlevels,'LevelOneFilter','nearsym5_7','FilterLength',14);
    for loop1=1:numlevels;
        ampNP{loop2,loop1}=abs(Yh{loop1});
        phaseDiff{loop2,loop1}=angle(Yh{loop1})-phaseYh{loop1};
    end
end

for loop2=1:sizer(2)
    sortval(:,loop2)=sort(data_array(count:length(data_array),loop2));
    stdval(loop2)=std(sortval(:,loop2));
end
num2sort=length(sortval);


for surrloop=1:numsurrogates
    disp(strcat('Making surrogate number_',num2str(surrloop)))
    
    %make a random dataset and take its imag and phases
    for loop1=1:sizer(2)
        [dummy,shuffind]=sort(rand(num2sort,1));
        shuffind=shuffind-1+count;	
        z(shuffind,loop1) = sortval(:,loop1);
    end
    z(1:count-1,:)=0;
    
    [Zl,Zh] = dualtree(z(:,primary),'Level',numlevels,'LevelOneFilter','nearsym5_7','FilterLength',14);
    %[Zl,Zh] = dtwavexfm(z,numlevels,'near_sym_b','qshift_b');
        
    for loop1=1:numlevels
        newphase{loop1}=angle(Zh{loop1});
    end
    
    
    amperror(1)=100;
    waverror(1)=100;   
    counter=1;
    
    while (amperror(counter) > accerror) && (waverror(counter) > accerror)
        %wavelet construction
        oldz=z;
        for loop1=1:numlevels
            newZh{loop1}=ampYh{loop1}.*exp(1i.*newphase{loop1});
        end
        for loop2=1:length(not_primary)
            for loop1=1:numlevels
                combphase=newphase{loop1}+phaseDiff{loop2,loop1};
                toohigh=find(combphase>pi);
                combphase(toohigh)=-(2*pi-combphase(toohigh));
                toolow=find(combphase<-pi);
                combphase(toolow)=(2*pi+combphase(toolow));
                newNP{loop2}{loop1}=ampNP{loop2,loop1}.*exp(1i.*combphase);
            end
        end
        
        %z=dtwaveifm(Yl,newZh,'near_sym_b','qshift_b');  
        z(:,primary)=idualtree(Yl{primary},newZh,'LevelOneFilter','nearsym5_7','FilterLength',14);
        for loop1=1:length(not_primary)
            z(:,not_primary(loop1))=idualtree(Yl{not_primary(loop1)},newNP{loop1},'LevelOneFilter','nearsym5_7','FilterLength',14);
        end
        wavdiff=(mean(abs(real(z)-real(oldz))))./stdval;
        waverror(counter+1) = mean(wavdiff);
        
        %impose original values
        oldz=z;
        data2sort=z(count:length(data_array),:);
        for loop1=1:sizer(2)
            [dummy,shuffind]=sort(real(data2sort(:,loop1)));
            shuffind=shuffind-1+count;	
            z(shuffind,loop1) = sortval(:,loop1);
        end
        z(1:count-1)=0;
        ampdiff=(mean(abs(real(z)-real(oldz))))./stdval;
        amperror(counter+1) = mean(ampdiff);
        
        %Wavelet step
        [nZl,nZh] = dualtree(z(:,primary),'Level',numlevels,'LevelOneFilter','nearsym5_7','FilterLength',14);  
        %[nZl,nZh] = dtwavexfm(z,numlevels,'near_sym_b','qshift_b');

        %get phases
        for loop1=1:numlevels;
            newphase{loop1}=angle(nZh{loop1});
        end
        
        toterror=amperror(counter+1)+waverror(counter+1);
        oldtoterr=amperror(counter)+waverror(counter);
        if abs((oldtoterr-toterror)/toterror) < (accerror/error_change);
            amperror(counter+1)=-1;
            waverror(counter+1)=-1;
        end
        counter=counter+1;
        clear nZh nZl
    end
    
    clear amperror waverror
    sgates{surrloop}=z(count:length(z),:);
    errors(surrloop,1)=toterror;
end