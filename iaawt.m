function [sgates,errors]=iaawt(data_array,numsurrogates,accerror,error_change)

%This is the IAAWT algorithm as presented in:

%Keylock, C.J. 2017. Multifractal surrogate-data generation algorithm that preserves pointwise 
%Hölder regularity structure, with initial applications to turbulence, Physical Review E 95, 032123, 
%https://doi.org/10.1103/PhysRevE.95.032123.

%If you use this code or a derivative thereof for your own work, I would be
%grateful if you could acknowledge this source.

%The algorithm makes use of the dual tree complex wavelet transform
%developed by Nick Kingsbury in Cambridge. At the moment, we use his
%toolbox. Instructions on how to access it are available at
%http://sigproc.eng.cam.ac.uk/Main/NGK
%You will need to have downloaded this toolbox and for it to be
%in your path for this algorithm to function.

%INPUTS

%data_array (length(data_array) by 1 data vector)
%This is a dyadic wavelet transform. If you pass it data that is not of
%length 2^numlevels where numlevels is an integer given by
%log(length(data_array)) / log(2) then it will zero pad up to the next
%dyadic scale.

%numsurrogates (the number of surrogates you wish to generate)

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

if nargin<4
    error_change=100;
end
if nargin<3
    accerror=.001;
end
if nargin<2
    numsurrogates=1;
end


%randomisation using Kingsbury wavelets
sizer=size(data_array);
if sizer(1)==1
    data_array=data_array';
    sizer=size(data_array);   
end
exactlevels=log(length(data_array))/log(2);
numlevels=(floor(log(length(data_array))/log(2)));

count=1;
if abs(numlevels-exactlevels)>eps
    %We need to zero pad
    temp=zeros(2^(numlevels+1),1);
    count=2^(numlevels+1)-length(data_array)+1;
    temp(count:length(temp),1)=data_array;
    data_array=temp;
end

%wavelet decomp
[Yl,Yh] = dtwavexfm(data_array,numlevels,'near_sym_b','qshift_b');

% real and amplitudes
for loop1=1:numlevels;
    ampYh{loop1}=abs(Yh{loop1});
    phaseYh{loop1}=angle(Yh{loop1});
end

sortval=sort(data_array(count:length(data_array)));
stdval=std(sortval);
num2sort=length(sortval);

for surrloop=1:numsurrogates
    disp(strcat('Making surrogate number_',num2str(surrloop)))
    
    %make a random dataset and take its imag and phases
    [dummy,shuffind]=sort(rand(num2sort,1));
    shuffind=shuffind-1+count;	
    z(shuffind,1) = sortval;
    z(1:count-1)=0;
    [Zl,Zh] = dtwavexfm(z,numlevels,'near_sym_b','qshift_b');
        
    for loop1=1:numlevels
        newphase{loop1}=angle(Zh{loop1});
    end
    
    
    amperror(1)=100;
    waverror(1)=100;   
    counter=1;
    
    while (amperror(counter) > accerror) & (waverror(counter) > accerror)
        %wavelet construction
        oldz=z;
        for loop1=1:numlevels
            newZh{loop1}=ampYh{loop1}.*exp(i.*newphase{loop1});
        end
        
        z=dtwaveifm(Yl,newZh,'near_sym_b','qshift_b');
        wavdiff=mean(mean(abs(real(z)-real(oldz))));
        waverror(counter+1) = wavdiff/stdval;
        
        %impose original values
        oldz=z;
        data2sort=z(count:length(data_array));
        [dummy,shuffind]=sort(real(data2sort));
        shuffind=shuffind-1+count;	
        z(shuffind,1) = sortval;
        z(1:count-1)=0;
        ampdiff=mean(mean(abs(real(z)-real(oldz))));
        amperror(counter+1) = ampdiff/stdval;
        
        %Wavelet step
        [nZl,nZh] = dtwavexfm(z,numlevels,'near_sym_b','qshift_b');
        %get phases and imag
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
    sgates(:,surrloop)=z(count+1:length(z));
    errors(surrloop,1)=toterror;
end