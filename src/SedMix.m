%% SedMix.m written by Brian J. Yanites see Yanites, B.J, G.E. Tucker, and R.S. Anderson (2008, in press) 
%%Numerical and analytical models of cosmogenic radionuclide dynamics in landslide-dominated drainage basins, JGR-Earth Surface.

    %Run this model right after running "CosmoLand" or after loading .mat file
    %created by "CosmoLand"
    
    %Adjust "Model Input" to match your basin parameters
    
    %Output
       % 'stormass' tracks the mass held in fluvial storage
       % 'storatoms' tracks the number of atoms held in fluvial storage
       % 'ptime_ero' percent time the system produces an erosion rate
                    % estimate within the threshold value
       % 'ptime_eb' percent time the system produces an erosion rate
                    % estimate below the threshold value
       % 'ptime_ea' percent time the system produces an erosion rate
                    % estimate above the threshold value
       
       
       
%% Model Input
Ch=.0002; %coefficient to relate sed depth to drainage area (includes converstion from km to m)
Cw=.010; %coefficient to relate sed width to drainage area (includes converstion from km to m)
%Ca=5*(10^3)%8624%  actual regression gives exponent of 4, but that makes mean cross section larger than max for small drainage areas.  -2.25;   %coefficient to relate probability of drainage area-includes input of km2 to m2--derived from Peikang Hsi Accum data--PDF coeffiecient
Acrit=10^5; %Critical area needed to maintain a channel in m^2
be=1.4; %exponent relating probability to basin area from Rod-Iturbe pg 134-136
As=0;
thresh=.5  %% threshold to calculate the 'Reliability of the sample' i.e. a 0.5 threshold is a CRN concentration within %50 of the long term average.
%%

%Drainage area
Dd=1/(Acrit^.5); 

%counter if running multiple storage scenarios
gb=0
gb=gb+1;

        %Capture variables from "____" model run.
        run_area(1,gb)=BA;
         mod_time(1,gb)=model_time;
        run_alpha(1,gb)=alpha;
        run_RI(1,gb)=1/retrn;
        run_backE(1,gb)=backgE;
        run_Etotal(1,gb)=Etotal;
        run_concvar(1,gb)=var(tatoms./massout);
        run_meanlsA(1,gb)=mean(Atrack);
 
        %calculate coeffiecients 
        Ca=(1-be)./(((run_area(1,gb).*1000000).^(1-be))-(As.^(1-be))); %PDF of Area coeffiecint  
        k=((Ch*Cw*Ca)); 
        
       % Calculate sediment resevoir variables
        maxxa(1,gb)=Ch.*Cw.*(run_area(1,gb).*1000000); %maximum cross sectional area
        mxa(1,gb)=((k)./(2-be)).*(((run_area(1,gb).*1000000).^(2-be))-(Acrit.^(2-be)));%% Mean cross sectional area  
        V=Dd.*(run_area(1,gb).*1000000).*mxa(1,gb); %fluvial storage volume
        vol(1,gb)=V; %track volume if running multiple mixing scenarios
        erom=run_Etotal(1,gb)*.01; %total erosion rate in m/yr
        Qsed=(erom*run_area(1,gb)*1000000*(rho*(100^3))); %sed flux in m^3/yr
        K(1,gb)=run_Etotal(1,gb).*(run_area(1,gb)*1000*1000*100*100)/(vol(1,gb)*(100^3));  % mixing constant
       
        if K(1,gb)>1
            K(1,gb)=1;
        end
        
        
 %% Set up storage arrays
        
        storatoms=zeros(1,length(tarray));
        stormass=storatoms;
        carryoveratoms=storatoms;
        carryovermass=storatoms;
        atomsout=storatoms;
        msout=storatoms;
        
%% Run storage algorithm         
        for ii=1:length(tarray)
            if ii==1 %no 'carryover' atoms or mass from previous year (first timestep)
                storatoms(1,ii)=tatoms(1,ii);
                stormass(1,ii)=massout(1,ii);
            elseif ii>1
                storatoms(1,ii)=carryoveratoms(1,ii-1)+tatoms(1,ii);
                stormass(1,ii)=carryovermass(1,ii-1)+massout(1,ii);
            end
           Qsed=stormass(1,ii)*K(1,gb); 
           atomsout(1,ii)=Qsed*(storatoms(1,ii)/stormass(1,ii));
            msout(1,ii)=Qsed;
%         
            %update carryover values for next iterations (i.e. storage
            %starting point for next timestep)
            carryoveratoms(1,ii)=storatoms(1,ii)-atomsout(1,ii);
            carryovermass(1,ii)=stormass(1,ii)-msout(1,ii);
        end
        
        %caclculates the CRN concentration for each year of the storage volume
        conqs_conc=storatoms./stormass;
               
                
        
        
     %% Calculate storage parameters   
        mt(1,gb)=vol(1,gb)/((Etotal/100)*mBA);%mixing time/residence time
        varnomix(1,gb)=var(tatoms(begincapture:end)./massout(begincapture:end));%variance with no mixing
        varsmix(1,gb)=var(conqs_conc(begincapture:end));%variance with mixing
        stdnomix(1,gb)=sqrt(varnomix(1,gb));%standard deviation with no mixing
        stdmix(1,gb)=sqrt(varsmix(1,gb)); %standard deviation with mixing
        wconc(1,gb)=sum(tatoms(begincapture:end))/sum(massout(begincapture:end));%total fluxed weighted concentration
        meanc_nw(1,gb)=mean(tatoms(begincapture:end)./massout(begincapture:end));
        meanbc_nw(1,gb)=mean(conqs_conc(begincapture:end));
        maxstormass(1,gb)=max(stormass); %maximum mass in storage
        minstormass(1,gb)=min(stormass); %minimum mass in storage time series
        meanstormass(1,gb)=mean(stormass); %mean mass in storage time series
        meanstormass2(1,gb)=mean(stormass(begincapture:end));
        econc=P*zstar/Etotal; %%expected concentration for volumetric erosion rate
        conqs_e=P.*zstar./conqs_conc(begincapture:end); %%calc erosion interp from bucket conc
        
        
      %% calculate system reliability  
        ero_rat=conqs_e./Etotal;%%ratio of measured vs. expected erosion
        e_meetc=find(ero_rat>=(1-thresh) & ero_rat<=(1+thresh)); %finds timesteps where erosion rate is within threshold
        ptime_ero(1,gb)=100*length(e_meetc)./length(conqs_conc(begincapture:end));%finds percent time the system is within threshold
        
        nm_er=P*zstar./(tatoms(begincapture:end)./massout(begincapture:end));
        nm_rat=nm_er./Etotal;
        nm_meetc=find(nm_rat>=(1-thresh) & nm_rat<=(1+thresh));
        no_ptime_ero(1,gb)=100*length(nm_meetc)./length(conqs_conc(begincapture:end));
        
        
        %find timesteps erosion rate is below and above
        e_below=find(ero_rat<(1-thresh));
        e_above=find(ero_rat>(1+thresh));
        
        %find percent time erosion rate is below and above
        ptime_eb(1,gb)=100.*length(e_below)./length(conqs_conc(begincapture:end));
        ptime_ea(1,gb)=100.*length(e_above)./length(conqs_conc(begincapture:end));
        
        
        
        %same as above but for % time concentration is within threshold
        conc_rat=conqs_conc./econc; %%ration of bucket concentration to expected
        a_meetc=find(conc_rat>=(1-thresh) & conc_rat<=(1+thresh)); %Finds where ratio is with +/- threshold of unity
        ptime_et(1,gb)=100*length(a_meetc)./length(conqs_conc(begincapture:end)); %%counts number found in line above and divides by number of years in sample