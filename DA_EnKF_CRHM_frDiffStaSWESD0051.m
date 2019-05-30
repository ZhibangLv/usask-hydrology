%%In this code, the interception distribution is varied in different month.

clear
clc
for o=[13] 


Rnum=1;                        %Number of runs
EnNumber=o;
InNum=1369;                        %Number of ensembles
days=2158;
for q=1:Rnum
    
for u=1:EnNumber
copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\OriginalFiles\test1.int',['D:\u of s\InterceptionDA\DAFRsdswe1\test2_',num2str(u),'.int'])
end
    



%define the direction and .prj files will used in this run
cd 'D:\u of s\InterceptionDA\DAFRsdswe1';
mkdir (['Output',num2str(EnNumber),'_',num2str(q)])  %Create a folder for output files
mkdir (['Driving',num2str(EnNumber),'_',num2str(q)])  %Create a folder for forcing files
mkdir (['States',num2str(EnNumber),'_',num2str(q)])  %Create a folder for forcing files
mkdir (['StatesB',num2str(EnNumber),'_',num2str(q)])  %Create a folder for forcing files
fname= 'MarmotCreek_modified_25Jul18_PointTest_Fisera_07-17.prj';                 % the .prj will actually run in the loop                                        
OBSWE_name='OBSWE1DValidation2.obs';           %this file contains measured date SWE. the first and last ones are added on purpose as they are the begin and end of whole run
OP_name='CRHM_output_1.txt';      %name of the simulation output file name
Totalobsfile= 'Marmot_Hourly_ArrayMetData_withT_g_1Oct08-30Sept14_update_12Jul2018.obs';

% determine parameters
Ci=2102;                   %specific heat of Ice [J kg¯1 K¯1 ]
pw=1000;                   %density of water [kg m-3]
Tm=0;                      %m is the melting point of ice (C)
Dmax=249;                  %maximum density [kg m-3]
Wmax=0.0001;               %maximum water content in snow
Tmax=-10;


%import the observed SCA data into matlab
OBSWE= importOBSWE(OBSWE_name);
Obs_time=OBSWE(:,1);

ObmeanRTSD= OBSWE(:,2);
ObsdRTSD= num2cell(cell2mat(OBSWE(:,2))*0.05);
ObmeanSFSD= OBSWE(:,3);
ObsdSFSD= num2cell(cell2mat(OBSWE(:,3))*0.05);
ObmeanNFSD= OBSWE(:,4);
ObsdNFSD= num2cell(cell2mat(OBSWE(:,4))*0.05);
ObmeanNF= OBSWE(:,5);
ObsdNF= num2cell(cell2mat(OBSWE(:,5))*0.05);
ObmeanRT= OBSWE(:,6);
ObsdRT= num2cell(cell2mat(OBSWE(:,6))*0.05);
ObmeanSFT= OBSWE(:,7);
ObsdSFT= num2cell(cell2mat(OBSWE(:,7))*0.05);
ObmeanSFB= OBSWE(:,8);
ObsdSFB= num2cell(cell2mat(OBSWE(:,8))*0.05);
ObmeanLF= OBSWE(:,9);
ObsdLF= num2cell(cell2mat(OBSWE(:,9))*0.05);
Obs_timeN= OBSWE(:,10);

% import the total obs file
Totalobs = importobs1(Totalobsfile,1,52593);
ObsHead = Totalobs (1:9,:); %creat the head file for obs file
ObsT=Totalobs(:,1); % the time of obs 

%define the dimension of MATRICES
EnsembleSWERT=zeros(length(Obs_time)+1,EnNumber);
EnsembleSWENF=zeros(length(Obs_time)+1,EnNumber);
EnsembleSWESFT=zeros(length(Obs_time)+1,EnNumber);
EnsembleSWESFB=zeros(length(Obs_time)+1,EnNumber);
EnsembleSWELF=zeros(length(Obs_time)+1,EnNumber);
EnsembleSDRT=zeros(length(Obs_time)+1,EnNumber);
EnsembleSDNF=zeros(length(Obs_time)+1,EnNumber);
EnsembleSDSFT=zeros(length(Obs_time)+1,EnNumber);
EnsembleSDSFB=zeros(length(Obs_time)+1,EnNumber);
EnsembleSDLF=zeros(length(Obs_time)+1,EnNumber);

EnmeanRT=zeros(length(Obs_time)+1,1);
EnmeanNF=zeros(length(Obs_time)+1,1);
EnmeanSFT=zeros(length(Obs_time)+1,1);
EnmeanSFB=zeros(length(Obs_time)+1,1);
EnmeanLF=zeros(length(Obs_time)+1,1);
EnmeanRTSD=zeros(length(Obs_time)+1,1);
EnmeanNFSD=zeros(length(Obs_time)+1,1);
EnmeanSFTSD=zeros(length(Obs_time)+1,1);
EnmeanSFBSD=zeros(length(Obs_time)+1,1);


ExRT=zeros(length(Obs_time)+1,EnNumber);
ExNF=zeros(length(Obs_time)+1,EnNumber);
ExSFT=zeros(length(Obs_time)+1,EnNumber);
ExSFB=zeros(length(Obs_time)+1,EnNumber);
ExLF=zeros(length(Obs_time)+1,EnNumber);
ExRTSD=zeros(length(Obs_time)+1,EnNumber);
ExNFSD=zeros(length(Obs_time)+1,EnNumber);
ExSFTSD=zeros(length(Obs_time)+1,EnNumber);
ExSFBSD=zeros(length(Obs_time)+1,EnNumber);


EyRT=zeros(length(Obs_time)+1,EnNumber);
EyNF=zeros(length(Obs_time)+1,EnNumber);
EySFT=zeros(length(Obs_time)+1,EnNumber);
EySFB=zeros(length(Obs_time)+1,EnNumber);
EyLF=zeros(length(Obs_time)+1,EnNumber);
EyRTSD=zeros(length(Obs_time)+1,EnNumber);
EyNFSD=zeros(length(Obs_time)+1,EnNumber);
EySFSD=zeros(length(Obs_time)+1,EnNumber);


ObforcastRT=zeros(length(Obs_time)+1,EnNumber);
ObforcastNF=zeros(length(Obs_time)+1,EnNumber);
ObforcastSFT=zeros(length(Obs_time)+1,EnNumber);
ObforcastSFB=zeros(length(Obs_time)+1,EnNumber);
ObforcastLF=zeros(length(Obs_time)+1,EnNumber);
ObforcastRTSD=zeros(length(Obs_time)+1,EnNumber);
ObforcastNFSD=zeros(length(Obs_time)+1,EnNumber);
ObforcastSFSD=zeros(length(Obs_time)+1,EnNumber);



PxxRT=zeros(length(Obs_time)+1,1);
PxxNF=zeros(length(Obs_time)+1,1);
PxxSFT=zeros(length(Obs_time)+1,1);
PxxSFB=zeros(length(Obs_time)+1,1);
PxxLF=zeros(length(Obs_time)+1,1);
PxxRTSD=zeros(length(Obs_time)+1,1);
PxxNFSD=zeros(length(Obs_time)+1,1);
PxxSFTSD=zeros(length(Obs_time)+1,1);
PxxSFBSD=zeros(length(Obs_time)+1,1);


PyyRT=zeros(length(Obs_time)+1,1);
PyyNF=zeros(length(Obs_time)+1,1);
PyySFT=zeros(length(Obs_time)+1,1);
PyySFB=zeros(length(Obs_time)+1,1);
PyyLF=zeros(length(Obs_time)+1,1);
PyyRTSD=zeros(length(Obs_time)+1,1);
PyyNFSD=zeros(length(Obs_time)+1,1);
PyySFSD=zeros(length(Obs_time)+1,1);


KRT=zeros(length(Obs_time)+1,1);
KNF=zeros(length(Obs_time)+1,1);
KSFT=zeros(length(Obs_time)+1,1);
KSFB=zeros(length(Obs_time)+1,1);
KLF=zeros(length(Obs_time)+1,1);
KRTSD=zeros(length(Obs_time)+1,1);
KNFSD=zeros(length(Obs_time)+1,1);
KSFTSD=zeros(length(Obs_time)+1,1);
KSFBSD=zeros(length(Obs_time)+1,1);

EnSWEART=zeros(length(Obs_time)+1,EnNumber);
EnSWEANF=zeros(length(Obs_time)+1,EnNumber);
EnSWEASFT=zeros(length(Obs_time)+1,EnNumber);
EnSWEASFB=zeros(length(Obs_time)+1,EnNumber);
EnSWEALF=zeros(length(Obs_time)+1,EnNumber);
EnSDART=zeros(length(Obs_time)+1,EnNumber);
EnSDANF=zeros(length(Obs_time)+1,EnNumber);
EnSDASFT=zeros(length(Obs_time)+1,EnNumber);
EnSDASFB=zeros(length(Obs_time)+1,EnNumber);
EnSDALF=zeros(length(Obs_time)+1,EnNumber);

EnSWEAmeanRT=zeros(length(Obs_time)+1,1);
EnSWEAmeanNF=zeros(length(Obs_time)+1,1);
EnSWEAmeanSFT=zeros(length(Obs_time)+1,1);
EnSWEAmeanSFB=zeros(length(Obs_time)+1,1);
EnSWEAmeanLF=zeros(length(Obs_time)+1,1);
EnSDAmeanRT=zeros(length(Obs_time)+1,1);
EnSDAmeanNF=zeros(length(Obs_time)+1,1);
EnSDAmeanSFT=zeros(length(Obs_time)+1,1);
EnSDAmeanSFB=zeros(length(Obs_time)+1,1);
EnSDAmeanLF=zeros(length(Obs_time)+1,1);
% the loop to run CRHM and update SWE day by day.

for i=2:length(Obs_time)
    
    %%Make a loop for ensemble simulation
    for j=1:EnNumber
    prjfile = importPRJ(fname);
    prjfile(22,1)=Obs_time(i-1);% Define the simulation start and end time
    prjfile(23,1)=Obs_time(i);
    
    
    prjfile(311,1)=cellstr(['test',num2str(i),'_',num2str(j),'.int']);
    prjfile(315,1)=cellstr(['test',num2str(i+1),'_',num2str(j),'.int']);
    STname=['test',num2str(i+1),'_',num2str(j),'.int']; 
    
    X = sprintf('%s*', ObsT{10:end});
    %     S = CStr2String(OPT, '*');
    V = sscanf(X, '%f*');  % Convert the cells to numbers in output time
    c=9+find(V==cell2mat(Obs_timeN(i-1))); % Find out the location of start Time of the snow fall
    d=9+find(V==cell2mat(Obs_timeN(i)));
    ObstorunN=Totalobs(c:d,:);
    e=d-c+1;
    

    
    ObstorunN=Totalobs(c:d,:);
    %T
    ObstorunN(:,3)=num2cell(5*(1-2*rand(e,1))+cell2mat(ObstorunN(:,3)));
    %rh
    ObstorunN(:,10)=num2cell(max(5,15*(1-2*rand(e,1))+cell2mat(ObstorunN(:,10))));
    %Tg
    ObstorunN(:,17)=num2cell(cell2mat(ObstorunN(:,17))*0.3.*(1-2*rand(e,1))+cell2mat(ObstorunN(:,17)));
    %u
    ObstorunN(:,24)=num2cell(max(0,3*(1-2*rand(e,1))+cell2mat(ObstorunN(:,24))));
%     ObstorunN(:,18)=num2cell(max(0,50*(1-2*rand(e,1))+cell2mat(ObstorunN(:,18))));
    %P
    ObstorunN(:,31)=num2cell(max(0,(cell2mat(ObstorunN(:,31))*0.5+0.01).*(1-2*rand(e,1))+cell2mat(ObstorunN(:,31))));
    %Qsi
    ObstorunN(:,38)=num2cell(max(0,(cell2mat(ObstorunN(:,38))/3.3).*(1-2*rand(e,1))+cell2mat(ObstorunN(:,38))));
    
    Obstorun=[ObsHead;ObstorunN];
 
    prjfile(18,1)={['D:\u of s\InterceptionDA\DAFRsdswe1\Mt',num2str(i),'_',num2str(j),'.obs';]};
    %convert NaN in CRHMobs into blank
       for k = 1:numel(prjfile)
          if isnan(prjfile{k})
             prjfile{k} = '';
          end
       end
            
    %save updated project file to a .prj file
    dlmcell('MarmotCreek_modified_25Jul18_PointTest_Fisera_07-17.prj',prjfile);
    
        %convert NaN in Obstorun into blank
       for k = 1:numel(Obstorun)
          if isnan(Obstorun{k})
             Obstorun{k} = '';
          end
       end
            
    %save updated obs file to a .obs file
    dlmcell(['Mt',num2str(i),'_',num2str(j),'.obs'],Obstorun);

    


    %Run CRHM
    system(['CRHM.exe '  fname ]);    
    %import the model simulation output into matlab 
    OpSWE = importOP(OP_name);
    OPT=OpSWE(:,1);
    OPSNOWNF=OpSWE(:,2);
    OPSNOWRT=OpSWE(:,3);
    OPSNOWSFT=OpSWE(:,4);
    OPSNOWSFB=OpSWE(:,5);
    OPSNOWLF=OpSWE(:,6);
    
    OPSNOWNFSD=OpSWE(:,7);
    OPSNOWRTSD=OpSWE(:,8);
    OPSNOWSFTSD=OpSWE(:,9);
    OPSNOWSFBSD=OpSWE(:,10);
    OPSNOWLFSD=OpSWE(:,11);
    
    OutSWENF=OPSNOWNF(end);
    OutSWERT=OPSNOWRT(end);
    OutSWESFT=OPSNOWSFT(end);
    OutSWESFB=OPSNOWSFB(end);
    OutSWELF=OPSNOWLF(end);
    OutSDNF=OPSNOWNFSD(end);
    OutSDRT=OPSNOWRTSD(end);
    OutSDSFT=OPSNOWSFTSD(end);
    OutSDSFB=OPSNOWSFBSD(end);
    OutSDLF=OPSNOWLFSD(end);
    
    EnsembleSWENF(i,j)=str2double(cell2mat(OutSWENF));
    EnsembleSWERT(i,j)=str2double(cell2mat(OutSWERT));
    EnsembleSWESFT(i,j)=str2double(cell2mat(OutSWESFT));
    EnsembleSWESFB(i,j)=str2double(cell2mat(OutSWESFB));
    EnsembleSWELF(i,j)=str2double(cell2mat(OutSWELF));
    EnsembleSDNF(i,j)=str2double(cell2mat(OutSDNF));
    EnsembleSDRT(i,j)=str2double(cell2mat(OutSDRT));
    EnsembleSDSFT(i,j)=str2double(cell2mat(OutSDSFT));
    EnsembleSDSFB(i,j)=str2double(cell2mat(OutSDSFB));
    EnsembleSDLF(i,j)=str2double(cell2mat(OutSDLF));

     % move the used driving data to Folder Drving
     copyfile(['D:\u of s\InterceptionDA\DAFRsdswe1\Mt',num2str(i),'_',num2str(j),'.obs'],['D:\u of s\InterceptionDA\DAFRsdswe1\Driving',num2str(EnNumber),'_',num2str(q)])
     % copythe output file of this run to folder Output and change its name
     copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\CRHM_output_1.txt',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
     movefile(['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q),'\CRHM_output_1.txt'], ['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q),'\CRHM_output_',num2str(i),'_',num2str(j),'.txt']);
     delete (['D:\u of s\InterceptionDA\DAFRsdswe1\Mt',num2str(i),'_',num2str(j),'.obs']);
     copyfile(['D:\u of s\InterceptionDA\DAFRsdswe1\test',num2str(i+1),'_',num2str(j),'.int'],['D:\u of s\InterceptionDA\DAFRsdswe1\StatesB',num2str(EnNumber),'_',num2str(q)])
    
     copyfile(['D:\u of s\InterceptionDA\DAFRsdswe1\test',num2str(i),'_',num2str(j),'.int'],['D:\u of s\InterceptionDA\DAFRsdswe1\States',num2str(EnNumber),'_',num2str(q)])
     
     if i > 3;
        delete (['D:\u of s\InterceptionDA\DAFRsdswe1\test',num2str(i-1),'_',num2str(j),'.int']);
     end
    
    end
    
    
 %% Calculate Kalman gain and update SWE
 
    %Cslculate ensemble mean for each state vector AND difference from mean
    %of each emsemble
    EnmeanRT(i)=mean(EnsembleSWERT(i,:));
    for j = 1:EnNumber 
    ExRT(i,j)=EnsembleSWERT(i,j)-EnmeanRT(i);
    end
    
    EnmeanRTSD(i)=mean(EnsembleSDRT(i,:));
    for j = 1:EnNumber 
    ExRTSD(i,j)=EnsembleSDRT(i,j)-EnmeanRTSD(i);
    end
    
    EnmeanNF(i)=mean(EnsembleSWENF(i,:));
    for j = 1:EnNumber 
    ExNF(i,j)=EnsembleSWENF(i,j)-EnmeanNF(i);
    end
    
    EnmeanNFSD(i)=mean(EnsembleSDNF(i,:));
    for j = 1:EnNumber 
    ExNFSD(i,j)=EnsembleSDNF(i,j)-EnmeanNFSD(i);
    end
    
    EnmeanSFB(i)=mean(EnsembleSWESFB(i,:));
    for j = 1:EnNumber 
    ExSFB(i,j)=EnsembleSWESFB(i,j)-EnmeanSFB(i);
    end
    
    EnmeanSFBSD(i)=mean(EnsembleSDSFB(i,:));
    for j = 1:EnNumber 
    ExSFBSD(i,j)=EnsembleSDSFB(i,j)-EnmeanSFBSD(i);
    end
    
    EnmeanSFT(i)=mean(EnsembleSWESFT(i,:));
    for j = 1:EnNumber 
    ExSFT(i,j)=EnsembleSWESFT(i,j)-EnmeanSFT(i);
    end
    
%     EnmeanSFTSD(i)=mean(EnsembleSDSFT(i,:));
%     for j = 1:EnNumber 
%     ExSFTSD(i,j)=EnsembleSDSFT(i,j)-EnmeanSFTSD(i);
%     end

    EnmeanLF(i)=mean(EnsembleSWELF(i,:));
    for j = 1:EnNumber 
    ExLF(i,j)=EnsembleSWELF(i,j)-EnmeanLF(i);
    end
   
    %generate forecast observation
    if cell2mat(ObmeanNF(i))<999
      ObforcastNF(i,:)=max(0,cell2mat(ObsdNF(i))*(1-2*rand(EnNumber,1))+cell2mat(ObmeanNF(i)));
      for j = 1:EnNumber 
      EyNF(i,j)=ObforcastNF(i,j)-cell2mat(ObmeanNF(i));
      end
    end
    
    if cell2mat(ObmeanRT(i))<999
       ObforcastRT(i,:)=max(0,cell2mat(ObsdRT(i))*(1-2*rand(EnNumber,1))+cell2mat(ObmeanRT(i)));
       for j = 1:EnNumber 
       EyRT (i,j)=ObforcastRT(i,j)-cell2mat(ObmeanRT(i));
       end
    end
   
    if cell2mat(ObmeanSFT(i))<999
       ObforcastSFT(i,:)=max(0,cell2mat(ObsdSFT(i))*(1-2*rand(EnNumber,1))+cell2mat(ObmeanSFT(i)));
       for j = 1:EnNumber 
       EySFT(i,j)=ObforcastSFT(i,j)-cell2mat(ObmeanSFT(i));
       end
    end
    
    if cell2mat(ObmeanSFB(i))<999
        ObforcastSFB(i,:)=max(0,cell2mat(ObsdSFB(i))*(1-2*rand(EnNumber,1))+cell2mat(ObmeanSFB(i)));
        for j = 1:EnNumber 
        EySFB(i,j)=ObforcastSFB(i,j)-cell2mat(ObmeanSFB(i));
        end
    end
    
    
    if cell2mat(ObmeanLF(i))<999
        ObforcastLF(i,:)=max(0,cell2mat(ObsdLF(i))*(1-2*rand(EnNumber,1))+cell2mat(ObmeanLF(i)));
        for j = 1:EnNumber 
        EyLF(i,j)=ObforcastLF(i,j)-cell2mat(ObmeanLF(i));
        end
    end
    
    if cell2mat(ObmeanNFSD(i))<999
      ObforcastNFSD(i,:)=max(0,cell2mat(ObsdNFSD(i))*(1-2*rand(EnNumber,1))+cell2mat(ObmeanNFSD(i)));
      for j = 1:EnNumber 
      EyNFSD(i,j)=ObforcastNFSD(i,j)-cell2mat(ObmeanNFSD(i));
      end
    end

    if cell2mat(ObmeanRTSD(i))<999
       ObforcastRTSD(i,:)=max(0,cell2mat(ObsdRTSD(i))*(1-2*rand(EnNumber,1))+cell2mat(ObmeanRTSD(i)));
       for j = 1:EnNumber 
       EyRTSD(i,j)=ObforcastRTSD(i,j)-cell2mat(ObmeanRTSD(i));
       end
    end
   
    if cell2mat(ObmeanSFSD(i))<999
       ObforcastSFSD(i,:)=max(0,cell2mat(ObsdSFSD(i))*(1-2*rand(EnNumber,1))+cell2mat(ObmeanSFSD(i)));
       for j = 1:EnNumber 
       EySFSD(i,j)=ObforcastSFSD(i,j)-cell2mat(ObmeanSFSD(i));
       end
    end
    

    
     
    %Calculate covariance of model(Pxx) and observation(Pyy)
    PxxRT(i)=ExRT(i,:)*ExRT(i,:)'/(EnNumber-1);
    PxxNF(i)=ExNF(i,:)*ExNF(i,:)'/(EnNumber-1);
    PxxSFT(i)=ExSFT(i,:)*ExSFT(i,:)'/(EnNumber-1);
    PxxSFB(i)=ExSFB(i,:)*ExSFB(i,:)'/(EnNumber-1);
    PxxLF(i)=ExLF(i,:)*ExLF(i,:)'/(EnNumber-1);
    PxxRTSD(i)=ExRTSD(i,:)*ExRTSD(i,:)'/(EnNumber-1);
    PxxNFSD(i)=ExNFSD(i,:)*ExNFSD(i,:)'/(EnNumber-1);
    PxxSFBSD(i)=ExSFBSD(i,:)*ExSFBSD(i,:)'/(EnNumber-1);
   
    
    PyyRT(i)=EyRT(i,:)*EyRT(i,:)'/(EnNumber-1);
    PyyNF(i)=EyNF(i,:)*EyNF(i,:)'/(EnNumber-1);
    PyySFT(i)=EySFT(i,:)*EySFT(i,:)'/(EnNumber-1);
    PyySFB(i)=EySFB(i,:)*EySFB(i,:)'/(EnNumber-1);
    PyyLF(i)=EyLF(i,:)*EyLF(i,:)'/(EnNumber-1);
    PyyRTSD(i)=EyRTSD(i,:)*EyRTSD(i,:)'/(EnNumber-1);
    PyyNFSD(i)=EyNFSD(i,:)*EyNFSD(i,:)'/(EnNumber-1);
    PyySFSD(i)=EySFSD(i,:)*EySFSD(i,:)'/(EnNumber-1);    
    
    %For those model ensemble all equals zero, use everage Pxx for that
    %date
    if  PxxRT(i)==0;
        PxxRT(i)=mean(PxxRT(2:i));
    end
    
    if  PxxNF(i)==0;
        PxxNF(i)=mean(PxxNF(2:i));
    end
    
    if  PxxSFT(i)==0;
        PxxSFT(i)=mean(PxxSFT(2:i));
    end
    
    if  PxxSFB(i)==0;
        PxxSFB(i)=mean(PxxSFB(2:i));
    end
    
    if  PxxLF(i)==0;
        PxxLF(i)=mean(PxxLF(2:i));
    end
    
    if  PxxRTSD(i)==0;
        PxxRTSD(i)=mean(PxxRTSD(2:i));
    end
    
    if  PxxNFSD(i)==0;
        PxxNFSD(i)=mean(PxxNFSD(2:i));
    end
    
    if  PxxSFBSD(i)==0;
        PxxSFBSD(i)=mean(PxxSFBSD(2:i));
    end
    

    
    
        %For those observation equals zero, use everage Pxx for that
    %date
    if  PyyRT(i)==0;
        PyyRT(i)=mean(PyyRT(2:i));
    end
    
    if  PyyNF(i)==0;
        PyyNF(i)=mean(PyyNF(2:i));
    end
    
    if  PxxSFT(i)==0;
        PyySFT(i)=mean(PyySFT(2:i));
    end
    
    if  PyySFB(i)==0;
        PyySFB(i)=mean(PyySFB(2:i));
    end
    
    if  PyyLF(i)==0;
        PyyLF(i)=mean(PyyLF(2:i));
    end
    
    if  PyyRTSD(i)==0;
        PyyRTSD(i)=mean(PyyRTSD(2:i));
    end
    
    if  PyyNFSD(i)==0;
        PyyNFSD(i)=mean(PyyNFSD(2:i));
    end
    
    if  PyySFSD(i)==0;
        PyySFSD(i)=mean(PyySFSD(2:i));
    end
    
    
    %Calculate Kalman gain for each variables    
    KRT(i)=PxxRT(i)*inv(PxxRT(i)+PyyRT(i));
    KNF(i)=PxxNF(i)*inv(PxxNF(i)+PyyNF(i));
    KSFT(i)=PxxSFT(i)*inv(PxxSFT(i)+PyySFT(i));
    KSFB(i)=PxxSFB(i)*inv(PxxSFB(i)+PyySFB(i));
    KLF(i)=PxxLF(i)*inv(PxxLF(i)+PyyLF(i));
    KRTSD(i)=PxxRTSD(i)*inv(PxxRTSD(i)+PyyRTSD(i));
    KNFSD(i)=PxxNFSD(i)*inv(PxxNFSD(i)+PyyNFSD(i));
    KSFBSD(i)=PxxSFBSD(i)*inv(PxxSFBSD(i)+PyySFSD(i));
    
    %Update state and calculate the updated ensemble mean
    if cell2mat(ObmeanRT(i))<999
    EnSWEART(i,:)=EnsembleSWERT(i,:)+KRT(i)*(ObforcastRT(i,:)-EnsembleSWERT(i,:));
    EnSWEAmeanRT(i)=mean(EnSWEART(i,:));
    end
    
    if cell2mat(ObmeanNF(i))<999
    EnSWEANF(i,:)=EnsembleSWENF(i,:)+KNF(i)*(ObforcastNF(i,:)-EnsembleSWENF(i,:));
    EnSWEAmeanNF(i)=mean(EnSWEANF(i,:));
    end
    
    if cell2mat(ObmeanSFT(i))<999
    EnSWEASFT(i,:)=EnsembleSWESFT(i,:)+KSFT(i)*(ObforcastSFT(i,:)-EnsembleSWESFT(i,:));
    EnSWEAmeanSFT(i)=mean(EnSWEASFT(i,:));
    end
    
    if cell2mat(ObmeanSFB(i))<999
    EnSWEASFB(i,:)=EnsembleSWESFB(i,:)+KSFB(i)*(ObforcastSFB(i,:)-EnsembleSWESFB(i,:));
    EnSWEAmeanSFB(i)=mean(EnSWEASFB(i,:));
    end
    
    if cell2mat(ObmeanLF(i))<999
    EnSWEALF(i,:)=EnsembleSWELF(i,:)+KLF(i)*(ObforcastLF(i,:)-EnsembleSWELF(i,:));  
    EnSWEAmeanLF(i)=mean(EnSWEALF(i,:));
    end
    
    if cell2mat(ObmeanRTSD(i))<999
    EnSDART(i,:)=EnsembleSDRT(i,:)+KRTSD(i)*(ObforcastRTSD(i,:)-EnsembleSDRT(i,:));
    EnSDAmeanRT(i)=mean(EnSDART(i,:));
    end
    
    if cell2mat(ObmeanNFSD(i))<999
    EnSDANF(i,:)=EnsembleSDNF(i,:)+KNFSD(i)*(ObforcastNFSD(i,:)-EnsembleSDNF(i,:));
    EnSDAmeanNF(i)=mean(EnSDANF(i,:));
    end
    
    if cell2mat(ObmeanSFSD(i))<999
    EnSDASFB(i,:)=EnsembleSDSFB(i,:)+KSFBSD(i)*(ObforcastSFSD(i,:)-EnsembleSDSFB(i,:));
    EnSDAmeanSFB(i)=mean(EnSDASFB(i,:));
    end
    
    if cell2mat(ObmeanSFSD(i))<999
    EnSDASFT(i,:)=EnsembleSDSFT(i,:)+KSFTSD(i)*(ObforcastSFSD(i,:)-EnsembleSDSFB(i,:))/1.13;
    EnSDAmeanSFB(i)=mean(EnSDASFT(i,:));
    end
    
    if cell2mat(ObmeanSFSD(i))<999
    EnSDALF(i,:)=EnsembleSDLF(i,:)+KSFTSD(i)*(ObforcastSFSD(i,:)-EnsembleSDSFB(i,:))*1.03;
    EnSDAmeanLF(i)=mean(EnSDALF(i,:));
    end
    
    display(i)

    
    for j=1:EnNumber
        
    %import the state file in to matlab
    STfile= importST(['test',num2str(i+1),'_',num2str(j),'.int']);
    
    
    for h=1:5;
    %T_s  snow pack Temperature
     if cell2mat(STfile(85,h))>0
       STfile(88,h)=num2cell(str2double(cell2mat(STfile(88,h))));
     else
       STfile(88,h)=num2cell(Tmax);
     end  
     
    %T_s_0  snow pack active layer Temperature
     if cell2mat(STfile(118,h))>0
       STfile(91,h)=num2cell(str2double(cell2mat(STfile(91,h))));
     else
       STfile(91,h)=num2cell(Tmax);
     end      
     
     %T_s_1  snow pack low layer Temperature
     if cell2mat(STfile(121,h))>0
       STfile(94,h)=num2cell(str2double(cell2mat(STfile(94,h))));
     else
       STfile(94,h)=num2cell(Tmax);
     end
    end 

    
    %Update SWE
    if cell2mat(ObmeanNF(i))<999
        STfile(85,1)=num2cell(max(0,EnSWEANF(i,j)));
    else
        STfile(85,1)=num2cell(str2double(cell2mat(STfile(85,1))));
    end
    
    if cell2mat(ObmeanRT(i))<999
        STfile(85,2)=num2cell(max(0,EnSWEART(i,j)));
    else
        STfile(85,2)=num2cell(str2double(cell2mat(STfile(85,2))));
    end
    
    if cell2mat(ObmeanSFT(i))<999
        STfile(85,3)=num2cell(max(0,EnSWEASFT(i,j)));
    else
        STfile(85,3)=num2cell(str2double(cell2mat(STfile(85,3))));
    end
    
    if cell2mat(ObmeanSFB(i))<999
        STfile(85,4)=num2cell(max(0,EnSWEASFB(i,j)));
    else
        STfile(85,4)=num2cell(str2double(cell2mat(STfile(85,4))));
    end
    
    if cell2mat(ObmeanLF(i))<999
        STfile(85,5)=num2cell(max(0,EnSWEALF(i,j)));
    else
        STfile(85,5)=num2cell(str2double(cell2mat(STfile(85,5))));
    end
    
    % update Snow depth
    if cell2mat(ObmeanNFSD(i))<999
        STfile(130,1)=num2cell(max(0,EnSDANF(i,j)));
    else
        STfile(130,1)=num2cell(str2double(cell2mat(STfile(130,1))));
    end
    
    if cell2mat(ObmeanRTSD(i))<999
        STfile(130,2)=num2cell(max(0,EnSDART(i,j)));
    else
        STfile(130,2)=num2cell(str2double(cell2mat(STfile(130,2))));
    end
    
    if cell2mat(ObmeanSFSD(i))<999
        STfile(130,3)=num2cell(max(0,EnSDASFT(i,j)));
    else
        STfile(130,3)=num2cell(str2double(cell2mat(STfile(130,3))));
    end
    
    if cell2mat(ObmeanSFSD(i))<999
        STfile(130,4)=num2cell(max(0,EnSDASFB(i,j)));
    else
        STfile(130,4)=num2cell(str2double(cell2mat(STfile(130,4))));
    end
        
    if cell2mat(ObmeanSFSD(i))<999
        STfile(130,5)=num2cell(max(0,EnSDALF(i,j)));
    else
        STfile(130,5)=num2cell(str2double(cell2mat(STfile(130,5))));
    end
    
    %change cell to double-Snow Cover
    SCNF=str2double(cell2mat(STfile(127,1)));
    SCRT=str2double(cell2mat(STfile(127,2)));
    SCSFT=str2double(cell2mat(STfile(127,3)));
    SCSFB=str2double(cell2mat(STfile(127,4)));
    SCLF=str2double(cell2mat(STfile(127,5)));
    
          %%Update other states for HRU NF
          
          if cell2mat(ObmeanNF(i))<999 & cell2mat(ObmeanNFSD(i))<999
          
                 if cell2mat(STfile(130,1))>0;
                    
                    %snow cover (1/0)
                    STfile(127,1)=num2cell('1');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,1))>0.1
                    STfile(115,1)=num2cell('2'); 
                    else
                    STfile(115,1)=num2cell('1');
                    end
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if cell2mat(STfile(85,1))>0
                    STfile(124,1)=num2cell(cell2mat(STfile(85,1))/cell2mat(STfile(130,1)));
                    else
                    STfile(124,1)=num2cell(Dmax);
                    end
                                       
                    
                    %SWE  (mm)
                    STfile(85,1)=num2cell(cell2mat(STfile(124,1))*cell2mat(STfile(130,1)));
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,1)=num2cell(min(cell2mat(STfile(130,1)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,1)=num2cell(max(cell2mat(STfile(130,1))-cell2mat(STfile(133,1)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,1)=num2cell(cell2mat(STfile(133,1))*cell2mat(STfile(124,1)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,1)=num2cell(cell2mat(STfile(136,1))*cell2mat(STfile(124,1)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,1)=num2cell(Ci*pw*(cell2mat(STfile(88,1))-Tm)*cell2mat(STfile(85,1))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,1)=num2cell(Ci*pw*(cell2mat(STfile(91,1))-Tm)*cell2mat(STfile(118,1))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,1)=num2cell(Ci*pw*(cell2mat(STfile(94,1))-Tm)*cell2mat(STfile(121,1))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,1)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,1)=num2cell(Wmax);
                    
                 elseif (cell2mat(STfile(130,1))==0) & (cell2mat(STfile(85,1))>0);
                   
                    %snow cover (1/0)
                    STfile(127,1)=num2cell('1');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,1)))>0
                    STfile(124,1)=num2cell(str2double(cell2mat(STfile(124,1))));
                    else
                    STfile(124,1)=num2cell(Dmax);
                    end
                    
                    
                    %SD  (mm)
                    STfile(130,1)=num2cell(cell2mat(STfile(85,1))/cell2mat(STfile(124,1)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,1))>0.1
                    STfile(115,1)=num2cell('2'); 
                    else
                    STfile(115,1)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,1)=num2cell(min(cell2mat(STfile(130,1)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,1)=num2cell(max(cell2mat(STfile(130,1))-cell2mat(STfile(133,1)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,1)=num2cell(cell2mat(STfile(133,1))*cell2mat(STfile(124,1)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,1)=num2cell(cell2mat(STfile(136,1))*cell2mat(STfile(124,1)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,1)=num2cell(Ci*pw*(cell2mat(STfile(88,1))-Tm)*cell2mat(STfile(85,1))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,1)=num2cell(Ci*pw*(cell2mat(STfile(91,1))-Tm)*cell2mat(STfile(118,1))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,1)=num2cell(Ci*pw*(cell2mat(STfile(94,1))-Tm)*cell2mat(STfile(121,1))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,1)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,1)=num2cell(Wmax);

                    
                 elseif   (cell2mat(STfile(130,1))==0)& (cell2mat(STfile(85,1))==0);

                     STfile(127,1)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');

 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,1)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,1)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,1)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,1)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,1)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,1)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,1)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,1)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,1)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,1)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,1)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,1)=num2cell('0');

                end
              
           elseif cell2mat(ObmeanNF(i))==999 & cell2mat(ObmeanNFSD(i))<999
          
                 if cell2mat(STfile(130,1))>0;
                    
                     %snow cover (1/0)
                    STfile(127,1)=num2cell('1');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,1))>0.1
                    STfile(115,1)=num2cell('2'); 
                    else
                    STfile(115,1)=num2cell('1');
                    end


                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,1)))>0
                    STfile(124,1)=num2cell(cell2mat(STfile(85,1))/cell2mat(STfile(130,1)));
                    else
                    STfile(124,1)=num2cell(Dmax);
                    end
                    
                    %SWE  (mm)
                    STfile(85,1)=num2cell(cell2mat(STfile(124,1))*cell2mat(STfile(130,1)));
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,1)=num2cell(min(cell2mat(STfile(130,1)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,1)=num2cell(max(cell2mat(STfile(130,1))-cell2mat(STfile(133,1)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,1)=num2cell(cell2mat(STfile(133,1))*cell2mat(STfile(124,1)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,1)=num2cell(cell2mat(STfile(136,1))*cell2mat(STfile(124,1)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,1)=num2cell(Ci*pw*(cell2mat(STfile(88,1))-Tm)*cell2mat(STfile(85,1))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,1)=num2cell(Ci*pw*(cell2mat(STfile(91,1))-Tm)*cell2mat(STfile(118,1))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,1)=num2cell(Ci*pw*(cell2mat(STfile(94,1))-Tm)*cell2mat(STfile(121,1))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,1)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,1)=num2cell(Wmax);

                 elseif (cell2mat(STfile(130,1))==0) & (cell2mat(STfile(85,1))>0);

                     
                    %snow cover (1/0)
                    STfile(127,1)=num2cell('1');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,1)))>0
                    STfile(124,1)=num2cell(str2double(cell2mat(STfile(124,1))));
                    else
                    STfile(124,1)=num2cell(Dmax);
                    end
                    
                    
                    %SD  (mm)
                    STfile(130,1)=num2cell(cell2mat(STfile(85,1))/cell2mat(STfile(124,1)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,1))>0.1
                    STfile(115,1)=num2cell('2'); 
                    else
                    STfile(115,1)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,1)=num2cell(min(cell2mat(STfile(130,1)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,1)=num2cell(max(cell2mat(STfile(130,1))-cell2mat(STfile(133,1)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,1)=num2cell(cell2mat(STfile(133,1))*cell2mat(STfile(124,1)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,1)=num2cell(cell2mat(STfile(136,1))*cell2mat(STfile(124,1)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,1)=num2cell(Ci*pw*(cell2mat(STfile(88,1))-Tm)*cell2mat(STfile(85,1))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,1)=num2cell(Ci*pw*(cell2mat(STfile(91,1))-Tm)*cell2mat(STfile(118,1))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,1)=num2cell(Ci*pw*(cell2mat(STfile(94,1))-Tm)*cell2mat(STfile(121,1))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,1)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,1)=num2cell(Wmax);

                    
                 elseif   (cell2mat(STfile(130,1))==0)& (cell2mat(STfile(85,1))==0);

                    STfile(127,1)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');

                    %SWE  (mm)
                    STfile(85,1)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,1)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,1)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,1)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,1)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,1)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,1)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,1)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,1)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,1)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,1)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,1)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,1)=num2cell('0');

                 end             
              
          elseif cell2mat(ObmeanNF(i))<999 & cell2mat(ObmeanNFSD(i))==999
          
                    if cell2mat(STfile(85,1))>0;

                    %snow cover (1/0)
                    STfile(127,1)=num2cell('1');

                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,1)))>0
                    STfile(124,1)=num2cell(str2double(cell2mat(STfile(124,1))));
                    else
                    STfile(124,1)=num2cell(Dmax);
                    end
                    
                    %SD (m)
                    STfile(130,1)=num2cell(cell2mat(STfile(85,1))/cell2mat(STfile(124,1)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,1))>0.1
                    STfile(115,1)=num2cell('2'); 
                    else
                    STfile(115,1)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,1)=num2cell(min(cell2mat(STfile(130,1)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,1)=num2cell(max(cell2mat(STfile(130,1))-cell2mat(STfile(133,1)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,1)=num2cell(cell2mat(STfile(133,1))*cell2mat(STfile(124,1)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,1)=num2cell(cell2mat(STfile(136,1))*cell2mat(STfile(124,1)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,1)=num2cell(Ci*pw*(cell2mat(STfile(88,1))-Tm)*cell2mat(STfile(85,1))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,1)=num2cell(Ci*pw*(cell2mat(STfile(91,1))-Tm)*cell2mat(STfile(118,1))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,1)=num2cell(Ci*pw*(cell2mat(STfile(94,1))-Tm)*cell2mat(STfile(121,1))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,1)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,1)=num2cell(Wmax);

                    
                    elseif   (cell2mat(STfile(85,1))==0);
                    STfile(127,1)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,1)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,1)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,1)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,1)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,1)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,1)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,1)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,1)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,1)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,1)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,1)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,1)=num2cell('0');

                    end              
              
          end      
 
              
              
           %%Update other states for RT
          if cell2mat(ObmeanRT(i))<999 & cell2mat(ObmeanRTSD(i))<999;
                    
                 if cell2mat(STfile(130,2))>0;

                    %snow cover (1/0)
                    STfile(127,2)=num2cell('1');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,2))>0.1
                    STfile(115,2)=num2cell('2'); 
                    else
                    STfile(115,2)=num2cell('1');
                    end
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if cell2mat(STfile(85,2))>0
                    STfile(124,2)=num2cell(cell2mat(STfile(85,2))/cell2mat(STfile(130,2)));
                    else
                    STfile(124,2)=num2cell(Dmax);
                    end
                    
                    %SWE  (mm)
                    STfile(85,2)=num2cell(cell2mat(STfile(124,2))*cell2mat(STfile(130,2)));
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,2)=num2cell(min(cell2mat(STfile(130,2)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,2)=num2cell(max(cell2mat(STfile(130,2))-cell2mat(STfile(133,2)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,2)=num2cell(cell2mat(STfile(133,2))*cell2mat(STfile(124,2)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,2)=num2cell(cell2mat(STfile(136,2))*cell2mat(STfile(124,2)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,2)=num2cell(Ci*pw*(cell2mat(STfile(88,2))-Tm)*cell2mat(STfile(85,2))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,2)=num2cell(Ci*pw*(cell2mat(STfile(91,2))-Tm)*cell2mat(STfile(118,2))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,2)=num2cell(Ci*pw*(cell2mat(STfile(94,2))-Tm)*cell2mat(STfile(121,2))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,2)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,2)=num2cell(Wmax);
                    
                    
                 elseif (cell2mat(STfile(130,2))==0) & (cell2mat(STfile(85,2))>0);

                    %snow cover (1/0)
                    STfile(127,2)=num2cell('1');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,2)))>0
                    STfile(124,2)=num2cell(str2double(cell2mat(STfile(124,2))));
                    else
                    STfile(124,2)=num2cell(Dmax);
                    end
                    
                    
                    %SD  (mm)
                    STfile(130,2)=num2cell(cell2mat(STfile(85,2))/cell2mat(STfile(124,2)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,2))>0.1
                    STfile(115,2)=num2cell('2'); 
                    else
                    STfile(115,2)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,2)=num2cell(min(cell2mat(STfile(130,2)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,2)=num2cell(max(cell2mat(STfile(130,2))-cell2mat(STfile(133,2)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,2)=num2cell(cell2mat(STfile(133,2))*cell2mat(STfile(124,2)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,2)=num2cell(cell2mat(STfile(136,2))*cell2mat(STfile(124,2)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,2)=num2cell(Ci*pw*(cell2mat(STfile(88,2))-Tm)*cell2mat(STfile(85,2))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,2)=num2cell(Ci*pw*(cell2mat(STfile(91,2))-Tm)*cell2mat(STfile(118,2))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,2)=num2cell(Ci*pw*(cell2mat(STfile(94,2))-Tm)*cell2mat(STfile(121,2))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,2)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,2)=num2cell(Wmax);


                    
                 elseif   (cell2mat(STfile(130,2))==0)&(cell2mat(STfile(85,2))==0);
                    STfile(127,2)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,2)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,2)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,2)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,2)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,2)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,2)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,2)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,2)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,2)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,2)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,2)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,2)=num2cell('0');

                 end
          
          elseif cell2mat(ObmeanRT(i))==999 & cell2mat(ObmeanRTSD(i))<999
                    
                 if cell2mat(STfile(130,2))>0;

                    %snow cover (1/0)
                    STfile(127,2)=num2cell('1');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,2))>0.1
                    STfile(115,2)=num2cell('2'); 
                    else
                    STfile(115,2)=num2cell('1');
                    end


                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,2)))>0
                    STfile(124,2)=num2cell(cell2mat(STfile(85,2))/cell2mat(STfile(130,2)));
                    else
                    STfile(124,2)=num2cell(Dmax);
                    end
                    
                    %SWE  (mm)
                    STfile(85,2)=num2cell(cell2mat(STfile(124,2))*cell2mat(STfile(130,2)));
                    
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,2)=num2cell(min(cell2mat(STfile(130,2)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,2)=num2cell(max(cell2mat(STfile(130,2))-cell2mat(STfile(133,2)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,2)=num2cell(cell2mat(STfile(133,2))*cell2mat(STfile(124,2)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,2)=num2cell(cell2mat(STfile(136,2))*cell2mat(STfile(124,2)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,2)=num2cell(Ci*pw*(cell2mat(STfile(88,2))-Tm)*cell2mat(STfile(85,2))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,2)=num2cell(Ci*pw*(cell2mat(STfile(91,2))-Tm)*cell2mat(STfile(118,2))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,2)=num2cell(Ci*pw*(cell2mat(STfile(94,2))-Tm)*cell2mat(STfile(121,2))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,2)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,2)=num2cell(Wmax);

                 elseif (cell2mat(STfile(130,2))==0) & (cell2mat(STfile(85,2))>0);

                     
                    %snow cover (1/0)
                    STfile(127,2)=num2cell('1');   
                     
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,2)))>0
                    STfile(124,2)=num2cell(str2double(cell2mat(STfile(124,2))));
                    else
                    STfile(124,2)=num2cell(Dmax);
                    end
                    
                    
                    %SD  (mm)
                    STfile(130,2)=num2cell(cell2mat(STfile(85,2))/cell2mat(STfile(124,2)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,2))>0.1
                    STfile(115,2)=num2cell('2'); 
                    else
                    STfile(115,2)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,2)=num2cell(min(cell2mat(STfile(130,2)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,2)=num2cell(max(cell2mat(STfile(130,2))-cell2mat(STfile(133,2)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,2)=num2cell(cell2mat(STfile(133,2))*cell2mat(STfile(124,2)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,2)=num2cell(cell2mat(STfile(136,2))*cell2mat(STfile(124,2)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,2)=num2cell(Ci*pw*(cell2mat(STfile(88,2))-Tm)*cell2mat(STfile(85,2))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,2)=num2cell(Ci*pw*(cell2mat(STfile(91,2))-Tm)*cell2mat(STfile(118,2))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,2)=num2cell(Ci*pw*(cell2mat(STfile(94,2))-Tm)*cell2mat(STfile(121,2))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,2)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,2)=num2cell(Wmax);
                    
                    
                 elseif   (cell2mat(STfile(130,2))==0)&(cell2mat(STfile(85,2))==0);
                    STfile(127,2)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
                    
                    %SWE  (mm)
                    STfile(85,2)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,2)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,2)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,2)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,2)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,2)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,2)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,2)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,2)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,2)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,2)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,2)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,2)=num2cell('0');

                 end

          elseif cell2mat(ObmeanRT(i))<999 & cell2mat(ObmeanRTSD(i))==999
                    
                 if cell2mat(STfile(85,2))>0;

                    %snow cover (1/0)
                    STfile(127,2)=num2cell('1');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,2)))>0
                    STfile(124,2)=num2cell(str2double(cell2mat(STfile(124,2))));
                    else
                    STfile(124,2)=num2cell(Dmax);
                    end
                    
                    %SD(m)
                    STfile(130,2)=num2cell(cell2mat(STfile(85,2))/cell2mat(STfile(124,2)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,2))>0.1
                    STfile(115,2)=num2cell('2'); 
                    else
                    STfile(115,2)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,2)=num2cell(min(cell2mat(STfile(130,2)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,2)=num2cell(max(cell2mat(STfile(130,2))-cell2mat(STfile(133,2)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,2)=num2cell(cell2mat(STfile(133,2))*cell2mat(STfile(124,2)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,2)=num2cell(cell2mat(STfile(136,2))*cell2mat(STfile(124,2)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,2)=num2cell(Ci*pw*(cell2mat(STfile(88,2))-Tm)*cell2mat(STfile(85,2))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,2)=num2cell(Ci*pw*(cell2mat(STfile(91,2))-Tm)*cell2mat(STfile(118,2))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,2)=num2cell(Ci*pw*(cell2mat(STfile(94,2))-Tm)*cell2mat(STfile(121,2))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,2)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,2)=num2cell(Wmax);

                    
                 elseif   (cell2mat(STfile(85,2))==0);
                    STfile(127,2)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,2)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,2)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,2)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,2)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,2)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,2)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,2)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,2)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,2)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,2)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,2)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,2)=num2cell('0');

                 end
          end         
                    
                    
       
            %%Update other states for SFT 
           
            if cell2mat(ObmeanSFT(i))<999 & cell2mat(ObmeanSFSD(i))<999
                    
                 if cell2mat(STfile(130,3))>0;

                    %snow cover (1/0)
                    STfile(127,3)=num2cell('1');                        
                        
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,3))>0.1
                    STfile(115,3)=num2cell('2'); 
                    else
                    STfile(115,3)=num2cell('1');
                    end
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if cell2mat(STfile(85,3))>0
                    STfile(124,3)=num2cell(cell2mat(STfile(85,3))/cell2mat(STfile(130,3)));
                    else
                    STfile(124,3)=num2cell(Dmax);
                    end
                    
                    %SWE  (mm)
                    STfile(85,3)=num2cell(cell2mat(STfile(124,3))*cell2mat(STfile(130,3)));
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,3)=num2cell(min(cell2mat(STfile(130,3)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,3)=num2cell(max(cell2mat(STfile(130,3))-cell2mat(STfile(133,3)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,3)=num2cell(cell2mat(STfile(133,3))*cell2mat(STfile(124,3)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,3)=num2cell(cell2mat(STfile(136,3))*cell2mat(STfile(124,3)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,3)=num2cell(Ci*pw*(cell2mat(STfile(88,3))-Tm)*cell2mat(STfile(85,3))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,3)=num2cell(Ci*pw*(cell2mat(STfile(91,3))-Tm)*cell2mat(STfile(118,3))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,3)=num2cell(Ci*pw*(cell2mat(STfile(94,3))-Tm)*cell2mat(STfile(121,3))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,3)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,3)=num2cell(Wmax);
                    
                    
                    
                 elseif (cell2mat(STfile(130,3))==0) & (cell2mat(STfile(85,3))>0);

                    %snow cover (1/0)
                    STfile(127,3)=num2cell('1');                        
                        
                        
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,3)))>0
                    STfile(124,3)=num2cell(str2double(cell2mat(STfile(124,3))));
                    else
                    STfile(124,3)=num2cell(Dmax);
                    end
                    
                    
                    %Snow depth  (m)
                    STfile(130,3)=num2cell(cell2mat(STfile(85,3))/cell2mat(STfile(124,3)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,3))>0.1
                    STfile(115,3)=num2cell('2'); 
                    else
                    STfile(115,3)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,3)=num2cell(min(cell2mat(STfile(130,3)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,3)=num2cell(max(cell2mat(STfile(130,3))-cell2mat(STfile(133,3)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,3)=num2cell(cell2mat(STfile(133,3))*cell2mat(STfile(124,3)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,3)=num2cell(cell2mat(STfile(136,2))*cell2mat(STfile(124,3)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,3)=num2cell(Ci*pw*(cell2mat(STfile(88,3))-Tm)*cell2mat(STfile(85,3))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,3)=num2cell(Ci*pw*(cell2mat(STfile(91,3))-Tm)*cell2mat(STfile(118,3))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,3)=num2cell(Ci*pw*(cell2mat(STfile(94,3))-Tm)*cell2mat(STfile(121,3))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,3)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,3)=num2cell(Wmax);
                    
                 elseif   (cell2mat(STfile(130,3))==0)&(cell2mat(STfile(85,3))==0);
                    STfile(127,3)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,3)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,3)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,3)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,3)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,3)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,3)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,3)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,3)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,3)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,3)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,3)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,3)=num2cell('0');

                 end
         
            elseif cell2mat(ObmeanSFT(i))==999 & cell2mat(ObmeanSFSD(i))<999
                    
                                        
                 if cell2mat(STfile(130,3))>0;

                    %snow cover (1/0)
                    STfile(127,3)=num2cell('1');                     
                     
                     
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,3))>0.1
                    STfile(115,3)=num2cell('2'); 
                    else
                    STfile(115,3)=num2cell('1');
                    end


                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,3)))>0
                    STfile(124,3)=num2cell(cell2mat(STfile(85,3))/cell2mat(STfile(130,3)));
                    else
                    STfile(124,3)=num2cell(Dmax);
                    end
                    
                    %SWE  (mm)
                    STfile(85,3)=num2cell(cell2mat(STfile(124,3))*cell2mat(STfile(130,3)));
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,3)=num2cell(min(cell2mat(STfile(130,3)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,3)=num2cell(max(cell2mat(STfile(130,3))-cell2mat(STfile(133,3)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,3)=num2cell(cell2mat(STfile(133,3))*cell2mat(STfile(124,3)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,3)=num2cell(cell2mat(STfile(136,3))*cell2mat(STfile(124,3)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,3)=num2cell(Ci*pw*(cell2mat(STfile(88,3))-Tm)*cell2mat(STfile(85,3))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,3)=num2cell(Ci*pw*(cell2mat(STfile(91,3))-Tm)*cell2mat(STfile(118,3))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,3)=num2cell(Ci*pw*(cell2mat(STfile(94,3))-Tm)*cell2mat(STfile(121,3))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,3)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,3)=num2cell(Wmax);

                elseif (cell2mat(STfile(130,3))==0) & (cell2mat(STfile(85,3))>0);

                    %snow cover (1/0)
                    STfile(127,3)=num2cell('1');                    
                    
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,3)))>0
                    STfile(124,3)=num2cell(str2double(cell2mat(STfile(124,3))));
                    else
                    STfile(124,3)=num2cell(Dmax);
                    end
                    
                    
                    %Snow depth  (m)
                    STfile(130,3)=num2cell(cell2mat(STfile(85,3))/cell2mat(STfile(124,3)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,3))>0.1
                    STfile(115,3)=num2cell('2'); 
                    else
                    STfile(115,3)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,3)=num2cell(min(cell2mat(STfile(130,3)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,3)=num2cell(max(cell2mat(STfile(130,3))-cell2mat(STfile(133,3)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,3)=num2cell(cell2mat(STfile(133,3))*cell2mat(STfile(124,3)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,3)=num2cell(cell2mat(STfile(136,3))*cell2mat(STfile(124,3)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,3)=num2cell(Ci*pw*(cell2mat(STfile(88,3))-Tm)*cell2mat(STfile(85,3))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,3)=num2cell(Ci*pw*(cell2mat(STfile(91,3))-Tm)*cell2mat(STfile(118,3))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,3)=num2cell(Ci*pw*(cell2mat(STfile(94,3))-Tm)*cell2mat(STfile(121,3))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,3)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,3)=num2cell(Wmax);
                    
                 elseif   (cell2mat(STfile(130,3))==0)&(cell2mat(STfile(85,3))==0);
                    STfile(127,3)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %SWE  mm
                    STfile(85,3)=num2cell('0');
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,3)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,3)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,3)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,3)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,3)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,3)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,3)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,3)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,3)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,3)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,3)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,3)=num2cell('0');

                 end
          
                    
            elseif cell2mat(ObmeanSFT(i))<999 & cell2mat(ObmeanSFSD(i))==999
                    
                                        
                    if cell2mat(STfile(85,3))>0;

                    %snow cover (1/0)
                    STfile(127,3)=num2cell('1');                        

                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,3)))>0
                    STfile(124,3)=num2cell(str2double(cell2mat(STfile(124,3))));
                    else
                    STfile(124,3)=num2cell(Dmax);
                    end
                    
                    %z_s active layer depth (m)
                    STfile(130,3)=num2cell(cell2mat(STfile(85,3))/cell2mat(STfile(124,3)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,3))>0.1
                    STfile(115,3)=num2cell('2'); 
                    else
                    STfile(115,3)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,3)=num2cell(min(cell2mat(STfile(130,3)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,3)=num2cell(max(cell2mat(STfile(130,3))-cell2mat(STfile(133,3)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,3)=num2cell(cell2mat(STfile(133,3))*cell2mat(STfile(124,3)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,3)=num2cell(cell2mat(STfile(136,3))*cell2mat(STfile(124,3)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,3)=num2cell(Ci*pw*(cell2mat(STfile(88,3))-Tm)*cell2mat(STfile(85,3))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,3)=num2cell(Ci*pw*(cell2mat(STfile(91,3))-Tm)*cell2mat(STfile(118,3))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,3)=num2cell(Ci*pw*(cell2mat(STfile(94,3))-Tm)*cell2mat(STfile(121,3))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,3)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,3)=num2cell(Wmax);

                    
                    elseif   (cell2mat(STfile(85,3))==0);
                    STfile(127,3)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,3)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,3)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,3)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,3)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,3)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,3)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,3)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,3)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,3)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,3)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,3)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,3)=num2cell('0');

                    end                    
            end
         
             %%Update other states for SFB 
           
            if cell2mat(ObmeanSFB(i))<999 & cell2mat(ObmeanSFSD(i))<999
                    
                  if cell2mat(STfile(130,4))>0;

                    %snow cover (1/0)
                    STfile(127,4)=num2cell('1');                        
                        
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,4))>0.1
                    STfile(115,4)=num2cell('2'); 
                    else
                    STfile(115,4)=num2cell('1');
                    end
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if cell2mat(STfile(85,4))>0
                    STfile(124,4)=num2cell(cell2mat(STfile(85,4))/cell2mat(STfile(130,4)));
                    else
                    STfile(124,4)=num2cell(Dmax);
                    end
                    
                    %SWE  (mm)
                    STfile(85,4)=num2cell(cell2mat(STfile(124,4))*cell2mat(STfile(130,4)));
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,4)=num2cell(min(cell2mat(STfile(130,4)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,4)=num2cell(max(cell2mat(STfile(130,4))-cell2mat(STfile(133,4)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,4)=num2cell(cell2mat(STfile(133,3))*cell2mat(STfile(124,4)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,4)=num2cell(cell2mat(STfile(136,3))*cell2mat(STfile(124,4)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,4)=num2cell(Ci*pw*(cell2mat(STfile(88,4))-Tm)*cell2mat(STfile(85,4))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,4)=num2cell(Ci*pw*(cell2mat(STfile(91,4))-Tm)*cell2mat(STfile(118,4))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,4)=num2cell(Ci*pw*(cell2mat(STfile(94,4))-Tm)*cell2mat(STfile(121,4))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,4)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,4)=num2cell(Wmax);
                    
                    
                    
                 elseif (cell2mat(STfile(130,4))==0) & (cell2mat(STfile(85,4))>0);

                    %snow cover (1/0)
                    STfile(127,4)=num2cell('1');                           
                        
                        
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,4)))>0
                    STfile(124,4)=num2cell(str2double(cell2mat(STfile(124,4))));
                    else
                    STfile(124,4)=num2cell(Dmax);
                    end
                    
                    
                    %Snow depth  (m)
                    STfile(130,4)=num2cell(cell2mat(STfile(85,4))/cell2mat(STfile(124,4)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,4))>0.1
                    STfile(115,4)=num2cell('2'); 
                    else
                    STfile(115,4)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,4)=num2cell(min(cell2mat(STfile(130,4)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,4)=num2cell(max(cell2mat(STfile(130,4))-cell2mat(STfile(133,4)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,4)=num2cell(cell2mat(STfile(133,4))*cell2mat(STfile(124,4)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,4)=num2cell(cell2mat(STfile(136,4))*cell2mat(STfile(124,4)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,4)=num2cell(Ci*pw*(cell2mat(STfile(88,4))-Tm)*cell2mat(STfile(85,4))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,4)=num2cell(Ci*pw*(cell2mat(STfile(91,4))-Tm)*cell2mat(STfile(118,4))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,4)=num2cell(Ci*pw*(cell2mat(STfile(94,4))-Tm)*cell2mat(STfile(121,4))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,4)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,4)=num2cell(Wmax);
                    
                  elseif   (cell2mat(STfile(130,4))==0)&(cell2mat(STfile(85,4))==0);
                    STfile(127,4)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,4)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,4)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,4)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,4)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,4)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,4)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,4)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,4)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,4)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,4)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,4)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,4)=num2cell('0');

                  end
         
            elseif cell2mat(ObmeanSFB(i))==999 & cell2mat(ObmeanSFSD(i))<999
                    
                                        
                 if cell2mat(STfile(130,4))>0;
                    
                    %snow cover (1/0)
                    STfile(127,4)=num2cell('1');   
                     
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,4))>0.1
                    STfile(115,4)=num2cell('2'); 
                    else
                    STfile(115,4)=num2cell('1');
                    end


                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,4)))>0
                    STfile(124,4)=num2cell(cell2mat(STfile(85,4))/cell2mat(STfile(130,4)));
                    else
                    STfile(124,4)=num2cell(Dmax);
                    end
                    
                    %SWE  (mm)
                    STfile(85,4)=num2cell(cell2mat(STfile(124,4))*cell2mat(STfile(130,4)));
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,4)=num2cell(min(cell2mat(STfile(130,4)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,4)=num2cell(max(cell2mat(STfile(130,4))-cell2mat(STfile(133,4)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,4)=num2cell(cell2mat(STfile(133,3))*cell2mat(STfile(124,4)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,4)=num2cell(cell2mat(STfile(136,3))*cell2mat(STfile(124,4)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,4)=num2cell(Ci*pw*(cell2mat(STfile(88,4))-Tm)*cell2mat(STfile(85,4))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,4)=num2cell(Ci*pw*(cell2mat(STfile(91,4))-Tm)*cell2mat(STfile(118,4))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,4)=num2cell(Ci*pw*(cell2mat(STfile(94,4))-Tm)*cell2mat(STfile(121,4))/1000);
                    
                    %h20  water  (0-1)
                    STfile(106,4)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,4)=num2cell(Wmax);

                elseif (cell2mat(STfile(130,4))==0) & (cell2mat(STfile(85,4))>0);

                    
                    %snow cover (1/0)
                    STfile(127,4)=num2cell('1');                       
                    
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,4)))>0
                    STfile(124,4)=num2cell(str2double(cell2mat(STfile(124,4))));
                    else
                    STfile(124,4)=num2cell(Dmax);
                    end
                    
                    
                    %Snow depth  (m)
                    STfile(130,4)=num2cell(cell2mat(STfile(85,4))/cell2mat(STfile(124,4)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,4))>0.1
                    STfile(115,4)=num2cell('2'); 
                    else
                    STfile(115,4)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,4)=num2cell(min(cell2mat(STfile(130,4)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,4)=num2cell(max(cell2mat(STfile(130,4))-cell2mat(STfile(133,4)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,4)=num2cell(cell2mat(STfile(133,4))*cell2mat(STfile(124,4)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,4)=num2cell(cell2mat(STfile(136,4))*cell2mat(STfile(124,4)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,4)=num2cell(Ci*pw*(cell2mat(STfile(88,4))-Tm)*cell2mat(STfile(85,4))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,4)=num2cell(Ci*pw*(cell2mat(STfile(91,4))-Tm)*cell2mat(STfile(118,4))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,4)=num2cell(Ci*pw*(cell2mat(STfile(94,4))-Tm)*cell2mat(STfile(121,4))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,4)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,4)=num2cell(Wmax);
                    
                 elseif   (cell2mat(STfile(130,4))==0)&(cell2mat(STfile(85,4))==0);
                    STfile(127,4)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %SWE  mm
                    STfile(85,4)=num2cell('0');
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,4)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,4)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,4)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,4)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,4)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,4)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,4)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,4)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,4)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,4)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,4)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,4)=num2cell('0');

                 end
          
                    
            elseif cell2mat(ObmeanSFB(i))<999 & cell2mat(ObmeanSFSD(i))==999
                    
                                        
                    if cell2mat(STfile(85,4))>0;

                    %snow cover (1/0)
                    STfile(127,4)=num2cell('1');   

                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,4)))>0
                    STfile(124,4)=num2cell(str2double(cell2mat(STfile(124,4))));
                    else
                    STfile(124,4)=num2cell(Dmax);
                    end
                    
                    %z_s active layer depth (m)
                    STfile(130,4)=num2cell(cell2mat(STfile(85,4))/cell2mat(STfile(124,4)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,4))>0.1
                    STfile(115,4)=num2cell('2'); 
                    else
                    STfile(115,4)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,4)=num2cell(min(cell2mat(STfile(130,4)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,4)=num2cell(max(cell2mat(STfile(130,4))-cell2mat(STfile(133,4)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,4)=num2cell(cell2mat(STfile(133,4))*cell2mat(STfile(124,4)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,4)=num2cell(cell2mat(STfile(136,4))*cell2mat(STfile(124,4)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,4)=num2cell(Ci*pw*(cell2mat(STfile(88,4))-Tm)*cell2mat(STfile(85,4))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,4)=num2cell(Ci*pw*(cell2mat(STfile(91,4))-Tm)*cell2mat(STfile(118,4))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,4)=num2cell(Ci*pw*(cell2mat(STfile(94,4))-Tm)*cell2mat(STfile(121,4))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,4)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,4)=num2cell(Wmax);

                    
                    elseif   (cell2mat(STfile(85,4))==0);
                    STfile(127,4)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,4)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,4)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,4)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,4)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,4)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,4)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,4)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,4)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,4)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,4)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,4)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,4)=num2cell('0');

                    end                    
            end         
          

             %%Update other states for LF 
           
            if cell2mat(ObmeanLF(i))<999 & cell2mat(ObmeanSFSD(i))<999
                    
                  if cell2mat(STfile(130,5))>0;

                    %snow cover (1/0)
                    STfile(127,5)=num2cell('1');                           
                        
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,5))>0.1
                    STfile(115,5)=num2cell('2'); 
                    else
                    STfile(115,5)=num2cell('1');
                    end
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if cell2mat(STfile(85,5))>0
                    STfile(124,5)=num2cell(cell2mat(STfile(85,5))/cell2mat(STfile(130,5)));
                    else
                    STfile(124,5)=num2cell(Dmax);
                    end
                    
                    %SWE  (mm)
                    STfile(85,5)=num2cell(cell2mat(STfile(124,5))*cell2mat(STfile(130,5)));
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,5)=num2cell(min(cell2mat(STfile(130,5)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,5)=num2cell(max(cell2mat(STfile(130,5))-cell2mat(STfile(133,5)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,5)=num2cell(cell2mat(STfile(133,5))*cell2mat(STfile(124,5)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,5)=num2cell(cell2mat(STfile(136,5))*cell2mat(STfile(124,5)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,5)=num2cell(Ci*pw*(cell2mat(STfile(88,5))-Tm)*cell2mat(STfile(85,5))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,5)=num2cell(Ci*pw*(cell2mat(STfile(91,5))-Tm)*cell2mat(STfile(118,5))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,5)=num2cell(Ci*pw*(cell2mat(STfile(94,5))-Tm)*cell2mat(STfile(121,5))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,5)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,5)=num2cell(Wmax);
                    
                    
                    
                  elseif (cell2mat(STfile(130,5))==0) & (cell2mat(STfile(85,5))>0);

                    %snow cover (1/0)
                    STfile(127,5)=num2cell('1');                       
                        
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,5)))>0
                    STfile(124,5)=num2cell(str2double(cell2mat(STfile(124,5))));
                    else
                    STfile(124,5)=num2cell(Dmax);
                    end
                    
                    
                    %Snow depth  (m)
                    STfile(130,5)=num2cell(cell2mat(STfile(85,5))/cell2mat(STfile(124,5)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,5))>0.1
                    STfile(115,5)=num2cell('2'); 
                    else
                    STfile(115,5)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,5)=num2cell(min(cell2mat(STfile(130,5)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,5)=num2cell(max(cell2mat(STfile(130,5))-cell2mat(STfile(133,5)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,5)=num2cell(cell2mat(STfile(133,5))*cell2mat(STfile(124,5)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,5)=num2cell(cell2mat(STfile(136,5))*cell2mat(STfile(124,5)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,5)=num2cell(Ci*pw*(cell2mat(STfile(88,5))-Tm)*cell2mat(STfile(85,5))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,5)=num2cell(Ci*pw*(cell2mat(STfile(91,5))-Tm)*cell2mat(STfile(118,5))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,5)=num2cell(Ci*pw*(cell2mat(STfile(94,5))-Tm)*cell2mat(STfile(121,5))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,5)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,5)=num2cell(Wmax);
                    
                  elseif   (cell2mat(STfile(130,5))==0)&(cell2mat(STfile(85,5))==0);
                    STfile(127,5)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,5)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,5)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,5)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,5)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,5)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,5)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,5)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,5)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,5)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,5)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,5)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,5)=num2cell('0');

                  end
         
            elseif cell2mat(ObmeanLF(i))==999 & cell2mat(ObmeanSFSD(i))<999
                    
                                        
                 if cell2mat(STfile(130,5))>0;

                    %snow cover (1/0)
                    STfile(127,5)=num2cell('1');                      
                     
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,5))>0.1
                    STfile(115,5)=num2cell('2'); 
                    else
                    STfile(115,5)=num2cell('1');
                    end


                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,5)))>0
                    STfile(124,5)=num2cell(cell2mat(STfile(85,5))/cell2mat(STfile(130,5)));
                    else
                    STfile(124,5)=num2cell(Dmax);
                    end
                    
                    %SWE  (mm)
                    STfile(85,5)=num2cell(cell2mat(STfile(124,5))*cell2mat(STfile(130,5)));
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,5)=num2cell(min(cell2mat(STfile(130,5)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,5)=num2cell(max(cell2mat(STfile(130,5))-cell2mat(STfile(133,5)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,5)=num2cell(cell2mat(STfile(133,5))*cell2mat(STfile(124,5)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,5)=num2cell(cell2mat(STfile(136,5))*cell2mat(STfile(124,5)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,5)=num2cell(Ci*pw*(cell2mat(STfile(88,5))-Tm)*cell2mat(STfile(85,5))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,5)=num2cell(Ci*pw*(cell2mat(STfile(91,5))-Tm)*cell2mat(STfile(118,5))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,5)=num2cell(Ci*pw*(cell2mat(STfile(94,5))-Tm)*cell2mat(STfile(121,5))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,5)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,5)=num2cell(Wmax);

                elseif (cell2mat(STfile(130,5))==0) & (cell2mat(STfile(85,5))>0);

                    %snow cover (1/0)
                    STfile(127,5)=num2cell('1');                     
                    
                    
                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,5)))>0
                    STfile(124,5)=num2cell(str2double(cell2mat(STfile(124,5))));
                    else
                    STfile(124,5)=num2cell(Dmax);
                    end
                    
                    
                    %Snow depth  (m)
                    STfile(130,5)=num2cell(cell2mat(STfile(85,5))/cell2mat(STfile(124,5)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,5))>0.1
                    STfile(115,5)=num2cell('2'); 
                    else
                    STfile(115,5)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,5)=num2cell(min(cell2mat(STfile(130,5)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,5)=num2cell(max(cell2mat(STfile(130,5))-cell2mat(STfile(133,5)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,5)=num2cell(cell2mat(STfile(133,5))*cell2mat(STfile(124,5)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,5)=num2cell(cell2mat(STfile(136,5))*cell2mat(STfile(124,5)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,5)=num2cell(Ci*pw*(cell2mat(STfile(88,5))-Tm)*cell2mat(STfile(85,5))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,5)=num2cell(Ci*pw*(cell2mat(STfile(91,5))-Tm)*cell2mat(STfile(118,5))/1000); 
                    
                    %cc_s_1  lower layer cold content  (J/m^2)
                    STfile(103,5)=num2cell(Ci*pw*(cell2mat(STfile(94,5))-Tm)*cell2mat(STfile(121,5))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,5)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,5)=num2cell(Wmax);
                    
                 elseif   (cell2mat(STfile(130,5))==0)&(cell2mat(STfile(85,5))==0);
                    STfile(127,5)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %SWE  mm
                    STfile(85,5)=num2cell('0');
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,5)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,5)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,5)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,5)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,5)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,5)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,5)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,5)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,5)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,5)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,5)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,5)=num2cell('0');

                 end
          
                    
            elseif cell2mat(ObmeanLF(i))<999 & cell2mat(ObmeanSFSD(i))==999
                    
                                        
                    if cell2mat(STfile(85,5))>0;

                    %snow cover (1/0)
                    STfile(127,5)=num2cell('1');                         

                    %rh0  average snowcover density   (kg/m^3)
                    if str2double(cell2mat(STfile(124,5)))>0
                    STfile(124,5)=num2cell(str2double(cell2mat(STfile(124,5))));
                    else
                    STfile(124,5)=num2cell(Dmax);
                    end
                    
                    %z_s active layer depth (m)
                    STfile(130,5)=num2cell(cell2mat(STfile(85,5))/cell2mat(STfile(124,5)));
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    if cell2mat(STfile(130,5))>0.1
                    STfile(115,5)=num2cell('2'); 
                    else
                    STfile(115,5)=num2cell('1');
                    end
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,5)=num2cell(min(cell2mat(STfile(130,5)),0.1));
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,5)=num2cell(max(cell2mat(STfile(130,5))-cell2mat(STfile(133,5)), 0));                
                    
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,5)=num2cell(cell2mat(STfile(133,5))*cell2mat(STfile(124,5)));
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,5)=num2cell(cell2mat(STfile(136,5))*cell2mat(STfile(124,5)));
                    
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,5)=num2cell(Ci*pw*(cell2mat(STfile(88,5))-Tm)*cell2mat(STfile(85,5))/1000);  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,5)=num2cell(Ci*pw*(cell2mat(STfile(91,5))-Tm)*cell2mat(STfile(118,5))/1000); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,5)=num2cell(Ci*pw*(cell2mat(STfile(94,5))-Tm)*cell2mat(STfile(121,5))/1000);

                    
                    %h20  water  (0-1)
                    STfile(106,5)=num2cell(Wmax); 
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,5)=num2cell(Wmax);

                    
                    elseif   (cell2mat(STfile(85,5))==0);
                    STfile(127,5)=num2cell('0');
%                     STfile(76,2)=num2cell('0');
%                     STfile(76,3)=num2cell('0');
 
                    %cc_s  snowcover's cold content   (J/m^2)
                    STfile(97,5)=num2cell('0');  %% specific heat of ice, the layer temperature and the specific mass
                    
                    %cc_s_0  active layer cold content  (J/m^2)
                    STfile(100,5)=num2cell('0'); 
                    
                    %cc_s_0  lower layer cold content  (J/m^2)
                    STfile(103,5)=num2cell('0');
                    
                    %h20  water  (0-1)
                    STfile(106,5)=num2cell('0');  
                    
                    %h20_sat  % of liquid H2O saturation (relative water, 0-1). 
                    STfile(109,5)=num2cell('0');
                    
                    %layer_count  number of layers in snowcover (0/1/2)
                    STfile(115,5)=num2cell('0');   
                    
                    %m_s_0  active layer specific mass (kg/m^2)
                    STfile(118,5)=num2cell('0');
                    
                    %m_s_1   lower layer specific mass (kg/m^2)
                    STfile(121,5)=num2cell('0');
                    
                    %rh0  average snowcover density   (kg/m^3)
                    STfile(124,5)=num2cell('0');
                    
                    % z_s total snowcover thickness (m)
                    STfile(130,5)=num2cell('0');
                    
                    %z_s_0 active layer depth (m)
                    STfile(133,5)=num2cell('0');
                   
                    %z_s_1 lower layer depth (m)
                    STfile(136,5)=num2cell('0');

                    end                    
            end         
          
            
            
            
            
    for k = 1:numel(STfile)
    if isnan(STfile{k})
           STfile{k} = '';
    end
    end
            
    %save updated state file to a .obs file
    dlmcell(['test',num2str(i+1),'_',num2str(j),'.int'],STfile);
     
    end


    

end 
    
    %save updated SWE to a .obs file
    dlmcell('UpdateNF.obs',num2cell(EnSWEAmeanNF));
    dlmcell('UpdateRT.obs',num2cell(EnSWEAmeanRT));
    dlmcell('UpdateSFT.obs',num2cell(EnSWEAmeanSFT));
    dlmcell('UpdateSFB.obs',num2cell(EnSWEAmeanSFB));
    dlmcell('UpdateLF.obs',num2cell(EnSWEAmeanLF));
    dlmcell('UpdateNFSD.obs',num2cell(EnSDAmeanNF));
    dlmcell('UpdateRTSD.obs',num2cell(EnSDAmeanRT));
    dlmcell('UpdateSFTSD.obs',num2cell(EnSDAmeanSFT));
    dlmcell('KNF.obs',num2cell(KNF));
    dlmcell('KRT.obs',num2cell(KRT));
    dlmcell('KSFT.obs',num2cell(KSFT));
    dlmcell('KSFB.obs',num2cell(KSFB));
    dlmcell('KLF.obs',num2cell(KLF));
    dlmcell('KNFSD.obs',num2cell(KNFSD));
    dlmcell('KRTSD.obs',num2cell(KRTSD));
    dlmcell('KSFTSD.obs',num2cell(KSFTSD));
    
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\UpdateNF.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\UpdateRT.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\UpdateSFT.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\UpdateSFB.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\UpdateLF.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\UpdateNFSD.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\UpdateRTSD.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\UpdateSFTSD.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])

    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\KNF.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\KRT.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\KSFT.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\KSFB.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\KLF.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\KNFSD.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\KRTSD.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\KSFTSD.obs',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    
    
    % Import the output files from Data assimilation to matlab 
    %and combine them together into a new file
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\dlmcell.m',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\importOP.m',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)])
    copyfile('D:\u of s\InterceptionDA\DAFRsdswe1\importOP1.m',['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)]) 
    
    cd (['D:\u of s\InterceptionDA\DAFRsdswe1\Output',num2str(EnNumber),'_',num2str(q)]);
    Outputtotal=[];
    DailyOutput=[];
    DailyOutput1=[];
    Output=[];
    Length=[];
    NFtotal=[];
    RTtotal=[];
    SFTtotal=[];
    SFBtotal=[];
    LFtotal=[];
    NFtotalSD=[];
    RTtotalSD=[];
    SFTtotalSD=[];
    SFBtotalSD=[];
    LFtotalSD=[];
    NFtotalDEN=[];
    RTtotalDEN=[];
    SFTtotalDEN=[];
    SFBtotalDEN=[];
    LFtotalDEN=[];
    NFDay=[];
    RTDay=[];
    SFTDay=[];
    SFBDay=[];
    LFDay=[];
    NFSDDay=[];
    RTSDDay=[];
    SFTSDDay=[];
    SFBSDDay=[];
    LFSDDay=[];
    NFDEDay=[];
    RTDEDay=[];
    SFTDEDay=[];
    SFBDEDay=[];
    LFDEDay=[];
    Date=[];
    
for j=1:EnNumber
  for i=2:InNum
      OP_name=['CRHM_output_',num2str(i),'_',num2str(j),'.txt'];
      OpSWE = importOP(OP_name);
      
      Output=[Output; OpSWE];
    
  end
  
  dlmcell(['Output',num2str(j),'.obs'],Output);
  Output=[];
end


for j=1:EnNumber;
    
      OP_name=['Output',num2str(j),'.obs'];
      OpSWE1 = importOP1(OP_name);
      Outputtotal=[Outputtotal, OpSWE1];
    
end

Time=Outputtotal(:,1);

for i=1:EnNumber;
NFtotal=[NFtotal,Outputtotal(:,16*(i-1)+2)];
RTtotal=[RTtotal,Outputtotal(:,16*(i-1)+3)];
SFTtotal=[SFTtotal,Outputtotal(:,16*(i-1)+4)];
SFBtotal=[SFBtotal,Outputtotal(:,16*(i-1)+5)];
LFtotal=[LFtotal,Outputtotal(:,16*(i-1)+6)];
NFtotalSD=[NFtotalSD,Outputtotal(:,16*(i-1)+7)];
RTtotalSD=[RTtotalSD,Outputtotal(:,16*(i-1)+8)];
SFTtotalSD=[SFTtotalSD,Outputtotal(:,16*(i-1)+9)];
SFBtotalSD=[SFBtotalSD,Outputtotal(:,16*(i-1)+10)];
LFtotalSD=[LFtotalSD,Outputtotal(:,16*(i-1)+11)];
NFtotalDEN=[NFtotalDEN,Outputtotal(:,16*(i-1)+12)];
RTtotalDEN=[RTtotalDEN,Outputtotal(:,16*(i-1)+13)];
SFTtotalDEN=[SFTtotalDEN,Outputtotal(:,16*(i-1)+14)];
SFBtotalDEN=[SFBtotalDEN,Outputtotal(:,16*(i-1)+15)];
LFtotalDEN=[LFtotalDEN,Outputtotal(:,16*(i-1)+16)];
end



NFtotalN=cellfun(@(x)str2double(x), NFtotal);
RTtotalN=cellfun(@(x)str2double(x), RTtotal);
SFTtotalN=cellfun(@(x)str2double(x), SFTtotal);
SFBtotalN=cellfun(@(x)str2double(x), SFBtotal);
LFtotalN=cellfun(@(x)str2double(x), LFtotal);
NFtotalSDN=cellfun(@(x)str2double(x), NFtotalSD);
RTtotalSDN=cellfun(@(x)str2double(x), RTtotalSD);
SFTtotalSDN=cellfun(@(x)str2double(x), SFTtotalSD);
SFBtotalSDN=cellfun(@(x)str2double(x), SFBtotalSD);
LFtotalSDN=cellfun(@(x)str2double(x), LFtotalSD);
NFtotalDENN=cellfun(@(x)str2double(x), NFtotalDEN);
RTtotalDENN=cellfun(@(x)str2double(x), RTtotalDEN);
SFTtotalDENN=cellfun(@(x)str2double(x), SFTtotalDEN);
SFBtotalDENN=cellfun(@(x)str2double(x), SFBtotalDEN);
LFtotalDENN=cellfun(@(x)str2double(x), LFtotalDEN);

NFmean=mean(NFtotalN,2);
RTmean=mean(RTtotalN,2);
SFTmean=mean(SFTtotalN,2);
SFBmean=mean(SFBtotalN,2);
LFmean=mean(LFtotalN,2);
NFmeanSD=mean(NFtotalSDN,2);
RTmeanSD=mean(RTtotalSDN,2);
SFTmeanSD=mean(SFTtotalSDN,2);
SFBmeanSD=mean(SFBtotalSDN,2);
LFmeanSD=mean(LFtotalSDN,2);
NFmeanDEN=mean(NFtotalDENN,2);
RTmeanDEN=mean(RTtotalDENN,2);
SFTmeanDEN=mean(SFTtotalDENN,2);
SFBmeanDEN=mean(SFBtotalDENN,2);
LFmeanDEN=mean(LFtotalDENN,2);


for d=1:days;
    Date(d,1)=cellfun(@(x)str2double(x), Time(d*24));
    NFDay(d,1)=mean(NFmean((d-1)*24+1:d*24));
    RTDay(d,1)=mean(RTmean((d-1)*24+1:d*24));
    SFTDay(d,1)=mean(SFTmean((d-1)*24+1:d*24));
    SFBDay(d,1)=mean(SFBmean((d-1)*24+1:d*24));
    LFDay(d,1)=mean(LFmean((d-1)*24+1:d*24));
    NFSDDay(d,1)=mean(NFmeanSD((d-1)*24+1:d*24));
    RTSDDay(d,1)=mean(RTmeanSD((d-1)*24+1:d*24));
    SFTSDDay(d,1)=mean(SFTmeanSD((d-1)*24+1:d*24));
    SFBSDDay(d,1)=mean(SFBmeanSD((d-1)*24+1:d*24));
    LFSDDay(d,1)=mean(LFmeanSD((d-1)*24+1:d*24));
    NFDEDay(d,1)=mean(NFmeanDEN((d-1)*24+1:d*24));
    RTDEDay(d,1)=mean(RTmeanDEN((d-1)*24+1:d*24));
    SFTDEDay(d,1)=mean(SFTmeanDEN((d-1)*24+1:d*24));
    SFBDEDay(d,1)=mean(SFBmeanDEN((d-1)*24+1:d*24));
    LFDEDay(d,1)=mean(LFmeanDEN((d-1)*24+1:d*24));
    
end


for d=1:days;
    NFDay1(d,1)=NFmean((d-1)*24+1);
    RTDay1(d,1)=RTmean((d-1)*24+1);
    SFTDay1(d,1)=SFTmean((d-1)*24+1);
    SFBDay1(d,1)=SFBmean((d-1)*24+1);
    LFDay1(d,1)=LFmean((d-1)*24+1);
    NFDaySD1(d,1)=NFmeanSD((d-1)*24+1);
    RTDaySD1(d,1)=RTmeanSD((d-1)*24+1);
    SFTDaySD1(d,1)=SFTmeanSD((d-1)*24+1);
    SFBDaySD1(d,1)=SFBmeanSD((d-1)*24+1);
    LFDaySD1(d,1)=LFmeanSD((d-1)*24+1);
    NFDayDE1(d,1)=NFmeanDEN((d-1)*24+1);
    RTDayDE1(d,1)=RTmeanDEN((d-1)*24+1);
    SFTDayDE1(d,1)=SFTmeanDEN((d-1)*24+1);
    SFBDayDE1(d,1)=SFBmeanDEN((d-1)*24+1);
    LFDayDE1(d,1)=LFmeanDEN((d-1)*24+1);
end

NFmeanC=num2cell(NFmean);
RTmeanC=num2cell(RTmean);
SFTmeanC=num2cell(SFTmean);
SFBmeanC=num2cell(SFBmean);
LFmeanC=num2cell(LFmean);
NFmeanSDC=num2cell(NFmeanSD);
RTmeanSDC=num2cell(RTmeanSD);
SFTmeanSDC=num2cell(SFTmeanSD);
SFBmeanSDC=num2cell(SFBmeanSD);
LFmeanSDC=num2cell(LFmeanSD);
NFmeanDENC=num2cell(NFmeanDEN);
RTmeanDENC=num2cell(RTmeanDEN);
SFTmeanDENC=num2cell(SFTmeanDEN);
SFBmeanDENC=num2cell(SFBmeanDEN);
LFmeanDENC=num2cell(LFmeanDEN);

DateC=num2cell(Date);
NFDayC=num2cell(NFDay);
RTDayC=num2cell(RTDay);
SFTDayC=num2cell(SFTDay);
SFBDayC=num2cell(SFBDay);
LFDayC=num2cell(LFDay);
NFDaySDC=num2cell(NFSDDay);
RTDaySDC=num2cell(RTSDDay);
SFTDaySDC=num2cell(SFTSDDay);
SFBDaySDC=num2cell(SFBSDDay);
LFDaySDC=num2cell(LFSDDay);
NFDayDENC=num2cell(NFDEDay);
RTDayDENC=num2cell(RTDEDay);
SFTDayDENC=num2cell(SFTDEDay);
SFBDayDENC=num2cell(SFBDEDay);
LFDayDENC=num2cell(LFDEDay);

NFDayC1=num2cell(NFDay1);
RTDayC1=num2cell(RTDay1);
SFTDayC1=num2cell(SFTDay1);
SFBDayC1=num2cell(SFBDay1);
LFDayC1=num2cell(LFDay1);
NFDaySDC1=num2cell(NFDaySD1);
RTDaySDC1=num2cell(RTDaySD1);
SFTDaySDC1=num2cell(SFTDaySD1);
SFBDaySDC1=num2cell(SFBDaySD1);
LFDaySDC1=num2cell(LFDaySD1);
NFDayDEC1=num2cell(NFDayDE1);
RTDayDEC1=num2cell(RTDayDE1);
SFTDayDEC1=num2cell(SFTDayDE1);
SFBDayDEC1=num2cell(SFBDayDE1);
LFDayDEC1=num2cell(LFDayDE1);

Outputtotal=[Outputtotal, NFmeanC];
Outputtotal=[Outputtotal, RTmeanC];
Outputtotal=[Outputtotal, SFTmeanC];
Outputtotal=[Outputtotal, SFBmeanC];
Outputtotal=[Outputtotal, LFmeanC];
Outputtotal=[Outputtotal, NFmeanSDC];
Outputtotal=[Outputtotal, RTmeanSDC];
Outputtotal=[Outputtotal, SFTmeanSDC];
Outputtotal=[Outputtotal, SFBmeanSDC];
Outputtotal=[Outputtotal, LFmeanSDC];
Outputtotal=[Outputtotal, NFmeanDENC];
Outputtotal=[Outputtotal, RTmeanDENC];
Outputtotal=[Outputtotal, SFTmeanDENC];
Outputtotal=[Outputtotal, SFBmeanDENC];
Outputtotal=[Outputtotal, LFmeanDENC];

DailyOutput=[DailyOutput, DateC];
DailyOutput=[DailyOutput, NFDayC];
DailyOutput=[DailyOutput, RTDayC];
DailyOutput=[DailyOutput, SFTDayC];
DailyOutput=[DailyOutput, SFBDayC];
DailyOutput=[DailyOutput, LFDayC];
DailyOutput=[DailyOutput, NFDaySDC];
DailyOutput=[DailyOutput, RTDaySDC];
DailyOutput=[DailyOutput, SFTDaySDC];
DailyOutput=[DailyOutput, SFBDaySDC];
DailyOutput=[DailyOutput, LFDaySDC];
DailyOutput=[DailyOutput, NFDayDENC];
DailyOutput=[DailyOutput, RTDayDENC];
DailyOutput=[DailyOutput, SFTDayDENC];
DailyOutput=[DailyOutput, SFBDayDENC];
DailyOutput=[DailyOutput, LFDayDENC];

DailyOutput1=[DailyOutput1,DateC];
DailyOutput1=[DailyOutput1,NFDayC1];
DailyOutput1=[DailyOutput1,RTDayC1];
DailyOutput1=[DailyOutput1,SFTDayC1];
DailyOutput1=[DailyOutput1,SFBDayC1];
DailyOutput1=[DailyOutput1,LFDayC1];
DailyOutput1=[DailyOutput1,NFDaySDC1];
DailyOutput1=[DailyOutput1,RTDaySDC1];
DailyOutput1=[DailyOutput1,SFTDaySDC1];
DailyOutput1=[DailyOutput1,SFBDaySDC1];
DailyOutput1=[DailyOutput1,LFDaySDC1];
DailyOutput1=[DailyOutput1,NFDayDEC1];
DailyOutput1=[DailyOutput1,RTDayDEC1];
DailyOutput1=[DailyOutput1,SFTDayDEC1];
DailyOutput1=[DailyOutput1,SFBDayDEC1];
DailyOutput1=[DailyOutput1,LFDayDEC1];

dlmcell('DailyOutput.obs',DailyOutput);
dlmcell('OutputDaily1.obs',DailyOutput1);
dlmcell('Outputtotal.obs',Outputtotal);
  
   
end

end
    disp('It is done now, Zhibang')
