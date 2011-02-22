clear all
%% 

%% 2-D model to study landslide cosmo mixing requirements/sampling
%% procedures.  For reference see Yanites, B.J, G.E. Tucker, and R.S. Anderson (2008,
%% in press) Numerical and analytical models of cosmogenic radionuclide
%% dynamics in landslide-dominated drainage basins, JGR-Earth Surface.
    % -Turn Matlab 'Cell Mode' to 'On'
    % -Change parameters in "MODEL INPUT".
    % -Use in conjunction with 'SedMix.m' to allow for fluvial storage
    % -Will save a .mat file with model outputs with format " 'ba_' BA '_backgE_' Background_E '_return_' Ri '_alpha_' alpha '.mat'  " where variable names are replaced with the actual variable values, can change on line 154
    % -Code calls powerlaw_distribution.m and uniform_distribution.m which were both written by D. Nathan Bradley
    % -need .m file 'powerlaw_distribution.m' and 'uniform_distribution.m' in directory to run
    %Output: 
        %'massout' records the annual mass coming off of hillslopes.
        %'tatoms' recores the annual number of atoms coming off of hillslopes in a given year.
        %'lsmass' records the annual mass contributed by landslides
        %'lsatoms' recoreds the annual atoms contributed by landslides
        %'Atrack' tracks landslide area
        %'Etotal' total mean erosion rate over entire model run
    
    %lines 286-299 uncomment to get distribution of landslide material in
    %variable "nfreqsls"

    %% MODEL INPUT
%Landslide scaling exponent
alpha=2.1;

%Background erosion rate in mm/yr
Background_E=.001
   
%Landslide recurrence interval for a 1km^2 basin in years
Ri=500

%Drainage area in km^2
BA=1

%scaling factor for lenght of landslide to depth, unitless.  Hovius et al.,
%1997 found et=0.05;
et=.05;

%total time period the model runs in years
model_time=200000;

%max landslide area m^2
Amax=1000000;

%minimum landslide area m^2
Amin=10000; 

%Cosmogenic production/decay, etc. parameters
att=160; %attenuantion of CRN production curve
rho=2;  %density of substrate/bedrock g/cm^3
P=6; %production rate atoms/g/yr
lambda=(4.62*(10^-7)); %CRN decay rate 1/yr
zstar=att/rho;

%model cell size in m^2
cellsize=10;

%if tracking histogram of landslide vs. background erosion, this variable
%controls the depth interval for averaging the concentration and then adding to the histogram.
hint=1; %%depth interval in cm for adding to histogram of concentrations

%Time step to begin capturing surface concentration snapshots
begincapture=50000;
%Number of snapshots to record
n_surfconc=50;


%% Simulated landscape set up: 2-D matricies
mBA=BA*1000*1000;  %convert BA to meters^2
cell=mBA/cellsize; %#of cells of size cellsize
n=round(sqrt(cell));%compute number of cells along x for x by x array
scape=zeros(n,n); %%create cosmo array to track surface concs
timesince=zeros(n,n);  %%array to track time since last landslide for each cell

%% model time set up
dt=1;
tarray=0:dt:model_time;


%% Landslide model set up
 %landslide frequency, value used for Poisson model
retrn=1./Ri;
% initial landslide volume
Vtotal=0;
% initial landslide count
landslidecount=0;


%% The following pulls random variables before time loop to save computation time
splcount=0;
splcount2=0;
splcount3=0;
splcount4=0;
plawd1=powerlaw_distribution(ceil(1+(length(tarray)*BA)/4),-alpha, Amin, Amax);
chance1=rand(1,ceil(1+(length(tarray)*BA)/4));
plawd2=powerlaw_distribution(ceil(1+(length(tarray)*BA)/4),-alpha, Amin, Amax);
chance2=rand(1,ceil(1+(length(tarray)*BA)/4));
plawd3=powerlaw_distribution(ceil(1+(length(tarray)*BA)/4),-alpha, Amin, Amax);
chance3=rand(1,ceil(1+(length(tarray)*BA)/4));
plawd4=powerlaw_distribution(ceil(1+(length(tarray)*BA)/4),-alpha, Amin, Amax);
chance4=rand(1,ceil(1+(length(tarray)*BA)/4));

%% array initiation
massout=zeros(1,length(tarray));
tatoms=zeros(1,length(tarray));
lsmass=zeros(1,length(tarray));
backmass=zeros(1,length(tarray));
conctout=zeros(1,length(tarray));
ErRate=zeros(1,length(tarray));


%% build up initial concentration array
t=5000000;
No=0;
E=.05;
Ntx=(No*(exp(-lambda*t)))+((P/(lambda+(rho*E/att)))*...
    (1-(exp(-(lambda+(rho*E/att))*t)))*(exp(-rho*(0)/att)));
scape=scape+Ntx;
% randomly distribute disturbance depths to scape
zint=powerlaw_distribution(n^2,-alpha, et*(Amin^.5), et*(Amax^.5));
    for scapefix=1:n
      gfix(scapefix,1:n)=zint(((scapefix*n)-(n-1)):(scapefix*n));
    end
    scape=scape.*exp(-gfix./zstar);
    initialscape=scape;

%% freq distribution set up
binsize=100;
bins=0:binsize:8e6;
binInd=bins/binsize;
nfreq=zeros(1,length(bins));
nfreqls=zeros(1,length(bins));

%% cap dist of surface concs
captime=model_time-begincapture;
capint=captime/n_surfconc;
cap_surfconc=begincapture-1;
surfdist=zeros(n_surfconc,n^2);
% surfdist2=zeros(n_surfconc,round(n*(n/3)));
% surfdist3=zeros(n_surfconc,round(n*(n/3)));
surf_freq=zeros(n_surfconc,length(bins));
surfcount=1;
%track_time=zeros(cell,n_surfconc);

%% convert variable units or calculate other variables if needed.  Model works in cm, g, yr, atoms
%Background_E to cm/yr
backgE=Background_E*.1;
report=model_time/100;
%% begin time loop
filename=['ba_' num2str(BA) '_backgE_' num2str(Background_E) '_return_' num2str(1/retrn) '_alpha_' num2str(alpha) '.mat']
for ttime=1:length(tarray);
   
    % track model progress
    if ttime>report
        percent_done=100*(1-(model_time-ttime)/model_time)
        report=report+(model_time/100);
    end
    
    time=ttime*dt;
    spatial_interval=0;
    % update surface concentrations for background erosion
    scape=(scape.*(exp(-lambda*dt)))+((P/(lambda+(rho*0/att)))*...
        (1-(exp(-(lambda+(rho*0/att))*dt)))*(exp(-rho*(0)/att)));
    concbackeout=-(scape.*(att/(backgE*rho))).*((exp(-(backgE*rho/att)))-1);
    backEatoms=(-(scape.*(att/(backgE*rho))).*((exp(-(backgE*rho/att)))-1)).*rho*backgE*(cellsize*10000); %%Amin=min cell size; convert to cm^2
    tatoms(1,ttime)=tatoms(1,ttime)+(sum(sum(backEatoms)));
    massout(1,ttime)=massout(1,ttime)+(backgE*rho*mBA*10000); %calculate mass added from background E- convert mBA to cm^2
     scape=scape.*exp(-backgE/zstar);%
     backmass(1,ttime)=backmass(1,ttime)+(backgE*rho*mBA*10000);
     backatms(1,ttime)=(sum(sum(backEatoms)));
    %Capture eroded concentration distributions and mass
%       if ttime>begincapture
%         binconc=ceil(concbackeout./binsize);
%         binconccol=binconc(:);
%         nv=hist(binconccol,binInd);%% each count= mass of rho*backgE
%         nfreq=nfreq+(nv.*(rho*backgE*(cellsize)*10000)); 
%         nv=[];
%       end
% %       
%      capture surface concentration distributions
       if ttime>begincapture && ttime>cap_surfconc;
          surfdist(surfcount,:)=scape(:);
%           surfdist2(surfcount,:)=scape((n/3)+1:2*n/3);
%           surfdist3(surfcount,:)=scape((2*n/3)+1:end)
            surf_freq(surfcount,:)=hist(scape(:),binInd);
            ff=timesince(:);
%            track_time(:,surfcount)=(time-ff);
            cap_surfconc=cap_surfconc+capint;
            surfcount=surfcount+1
       end
    
    
    
   %% Space loop for landslide decision making (i.e. Stark and Hovius data for # l.s. in 1 km^2)
   
    for spatial_interval=1:BA;

      %%finds previously pulled random variable for timestep and spatial
      %%interval.  Necessary when model_time and/or drainage basin size is
      %%large.  Variables are for conditions of l.s. generation below
        if splcount<((length(tarray))*BA)/4
            splcount=splcount+1;
            r=plawd1(splcount);
            g=chance1(1,splcount);
        elseif splcount==((length(tarray))*BA)/4
            r=plawd1(splcount);
            g=chance1(1,splcount);
            splcount=splcount+1;
        elseif splcount>((length(tarray))*BA)/4 && splcount2<((length(tarray))*BA)/4
            splcount2=splcount2+1;
            r=plawd2(splcount2);
            g=chance2(1,splcount2);
        elseif splcount2==((length(tarray))*BA)/4
            r=plawd2(splcount2);
            g=chance2(1,splcount2);
            splcount2=splcount2+1;
        elseif splcount2>((length(tarray))*BA)/4 && splcount3<((length(tarray))*BA)/4
            splcount3=splcount3+1;
            r=plawd3(splcount3);
            g=chance3(1,splcount3);
        elseif splcount3==((length(tarray))*BA)/4
            r=plawd3(splcount3);
            g=chance3(1,splcount3);
            splcount3=splcount3+1;
        elseif splcount3>((length(tarray))*BA)/4 %& splcount4<((length(tarray))*BA)/4
            splcount4=splcount4+1;
            r=plawd4(splcount4);
            g=chance4(1,splcount4);
        end
        
        %% Generate landslide if conditions met
        if g<retrn && r>Amin% 
            landslidecount=landslidecount+1;
            SA=r;
            Atrack(1,landslidecount)=SA; %track landslide area
            H=et*(SA^.5); %landslide depth
            volume=SA*H;
            Vtotal=volume+Vtotal; %%keeps track of total volume landslided away
            Hc=H*100; %%converts depth to cm
            depth(1,landslidecount)=H;
            
            %pick landslide location
            xloc=(round(rand*(n-1)))+1;
            yloc=(round(rand*(n-1)))+1;
            
            %ls size
            disturbed_cells=SA/(cellsize);
            sdc=sqrt(disturbed_cells);
            
            %%fit size to array and location
            ydir=(round((sdc)));
            xdir=(round((sdc)));
            exx=[];
            exy=[];
            if (ydir-1+yloc)>n
                exy=yloc-(n-ydir);
                ydir=n-yloc;
            end
            if (xdir-1+xloc)>n
                exx=xloc-(n-xdir);
                xdir=n-xloc;
            end
            locxarray=[xloc:1:(xloc+xdir-1) 1:1:exx];
            locyarray=[yloc:1:(yloc+ydir-1) 1:1:exy];
            
            %build landslide depth array
            depth_p=[0:hint:floor(Hc) Hc];
            
            % get surface concentrations at landslide location 
            lsscape(1:length(locxarray),1:length(locyarray))=scape(locxarray,locyarray);
            
            %measure atoms and mass eroded by landslide 
                for r1=1:length(locxarray)
                    for r2=1:length(locyarray)
                        locconc=lsscape(r1,r2);
                        meanconc=-(locconc*zstar/Hc)*((exp(-Hc/zstar)-1));
                        lsatoms=meanconc*rho*Hc*(cellsize)*10000;
                        tatoms(1,ttime)=tatoms(1,ttime)+lsatoms;
                        massout(1,ttime)=massout(1,ttime)+(rho*Hc*(cellsize)*10000);
                        lsmass(1,ttime)=lsmass(1,ttime)+(rho*Hc*(cellsize)*10000);
                        
                        % Capture concentration distribution of landslide
%                        eroded material.
%                           if ttime>begincapture
%                             dzconcs=locconc.*(exp(-depth_p(1:end-1)./zstar));
%                             binlsd=ceil(dzconcs./binsize);
%                             nv=hist(binlsd,binInd);
%                             nfreqls=nfreqls+(nv.*(rho*hint*(cellsize)*10000));
%                             nv=[];
%                             boconc=locconc.*(exp(-depth_p(end)./zstar));
%                             binlsdl=ceil(boconc./binsize);
%                             nvl=hist(binlsdl,binInd);
%                             nfreqls=nfreqls+(nvl.*((depth_p(end)-depth_p(end-1))*rho*(cellsize)*10000));
%                             nvl=[];
%                           end
                          
                                          
                                     
                    end
                end
                
            %% fix landscape surface concs to new surface concs based on
            %%l.s. depth
            timesince(locxarray,locyarray)=time;
            scape(locxarray,locyarray)=scape(locxarray,locyarray).*exp(-Hc/zstar);%adjust SCAPE concentration to new conc based on depht of landslide
            lsscape=[];
            locxarray=[];
            locyarray=[];
        end
    end
end


 tmass=sum(massout); %calculates total mass eroded
 Etotal=tmass/(rho*mBA*10000*time) %reports total long term erosion rate
 save(filename) %saves .mat file


    