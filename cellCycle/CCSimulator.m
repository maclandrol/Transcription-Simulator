function [liste, handle]= CCSimulator(g1, g2, m, arnnmber, protnmber, time, ncycle);
%Emmanuel Noutahi- Juillet 2013 Matlab R2013a
%Stochastic Biology Process Simulation with the direct method of Gillepsie
%(Gillespie,J.Phys.Chem.,1977).
%Promoter A --> Promoter I
%Promoter I --> Promoter A
%Promoter A --> ARNm (transcription)
%ARNm --> 0 (ARN degradation)
%ARNm --> Protein (translation)
%Protein --> 0 (Protein degradation)
%This function works with burst mode and continuous mode
clc;

assert (numel(g1)==numel(g2) && numel(g2)==numel(m), 'Error, paramètres manquants');
if ((g1(5)~=0 || g1(6) ~=0) && (g2(5)~=0 || (g2(6)~=0)) && (m(5)~=0 || m(6)~=0))
    
    [liste, handle] = burst(g1, g2, m, arnnmber, protnmber, time, ncycle);  
else
    [liste, handle] = cont(g1, g2, m, arnnmber, protnmber, time, ncycle);
end

end

function [liste, handle]= cont(g1, g2, m, arnnmber, protnmber, time, ncycle)

t_int=1
simul.time={};
simul.prot={};
simul.arn={};
simul.prom={};
simul.tau={};

if((g1(2)~=0 || g1(4) ~=0) && (g2(2)~=0 || (g2(4)~=0)) && (m(2)~=0 || m(4)~=0))
    
    for n=1:ncycle
        
        %% Phase 1 : G1/S
        %%%  Thattai - van Oudenaarden model (single gene model)
        %%% reaction prom -> ARN + prom (k1)
        %%% reaction ARN -> 0 (k2)
          
        t=0;
        arnnumber=arnnmber; protnumber=protnmber;
        t_end= g1(end);  % End Time
        trans= g1(1);
        trad=g1(2);
        darn=g1(3);
        dprot=g1(4);
        prom=1; % Only one promoters, and it's always actif
        %%% array to store results
        i=1; %indice of array
        t_array(i)=t;
        tau_array(i)=0;
        prom_array(i)=prom;
        arn_array(i)=arnnumber;
        prot_array(i)=protnumber; %protein array and initial value
        tic;
        
        while t< t_end
            
            %%% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            a= [trans*prom arnnumber*darn arnnumber*trad protnumber*dprot];
            a0= sum(a(:));
            
            %%% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
            while(r1==0 || r2==0), r1=rand; r2=rand; end
            tau= ((1/a0)* (log(1/r1)));
            
            j=1; mu=0; amu=0;
            
            while amu < r2*a0,
                mu = mu + 1;
                amu = amu + a(j);
                j = j + 1;
            end
            
            %%% Step 3 : Using the tau and r1 values obtained in step 2, increase
            % t by tau, and adjust the molecular population levels to reflect
            % the occurrence of one R reaction
            t = t+tau;
            if mu == 1 %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                arnnumber = arnnumber + 1;
            elseif mu == 2 %mRNA decay
                arnnumber = arnnumber - 1;
                if(arnnumber<=0)
                    arnnumber=0;
                end
            elseif mu == 3 % translation : protein production
                protnumber=protnumber+1;
            elseif mu == 4 %protein decay
                protnumber= protnumber - 1;
                if(protnumber<=0)
                    protnumber=0;
                end
            end
            % update  molecule number, species and time for later output
            if t >= i*t_int
                i=i+1;
                t_array(i)=t;
                prom_array(i)=prom;
                arn_array(i)=arnnumber;
                prot_array(i)= protnumber;
                tau_array(i)=tau;
            end
        end
        
        exc_time= toc;
        fprintf('G1/S simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        G1mean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        G1noise_mRNA_theoretical = sqrt(1/G1mean_mRNA_theoretical); %steady-state noise mRNA
        G1mean_protein_theoretical = (trans*trad)/(darn*dprot); %steady-state mean protein
        G1noise_protein_theoretical = sqrt((1/G1mean_protein_theoretical)+ (1/(G1mean_mRNA_theoretical*(1+(darn/dprot))))); %steady-state noise protein
        
        fprintf('\n-----------------G1/S Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Protein)= %.2f\nnoise(Protein)= %.2f\n',  G1mean_mRNA_theoretical , G1noise_mRNA_theoretical,G1mean_protein_theoretical,G1noise_protein_theoretical);
        
        %%% simulation statistics
        G1mean_mRNA_simulation = mean(arn_array);
        G1stdev_mRNA_simulation = std(arn_array);
        G1noise_mRNA_simulation = G1stdev_mRNA_simulation/G1mean_mRNA_simulation;
        G1mean_protein_simulation = mean(prot_array);
        G1stdev_protein_simulation = std(prot_array);
        G1noise_protein_simulation = G1stdev_protein_simulation/G1mean_protein_simulation;
        
        fprintf('\n-----------------G1/S Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Prot)= %.2f\nnoise(Prot)= %.2f\n',  G1mean_mRNA_simulation , G1noise_mRNA_simulation,G1mean_protein_simulation, G1noise_protein_simulation);
        
        
        %% G2 phase
        
        totG1 = i;
        t_end= g1(end)+g2(end);  % End Time
        trans= g2(1);
        trad=g2(2);
        darn=g2(3);
        dprot=g2(4);
        prom=1; %Two promoters at the G2 phase
        %%% array to store results
        tic;
        while t< t_end
            
            %%% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            a= [trans*prom arnnumber*darn arnnumber*trad protnumber*dprot];
            a0= sum(a(:));
            
            %%% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
            while(r1==0 || r2==0), r1=rand; r2=rand; end
            tau= ((1/a0)* (log(1/r1)));
            
            j=1; mu=0; amu=0;
            
            while amu < r2*a0,
                mu = mu + 1;
                amu = amu + a(j);
                j = j + 1;
            end
            
            %%% Step 3 : Using the tau and r1 values obtained in step 2, increase
            % t by tau, and adjust the molecular population levels to reflect
            % the occurrence of one R reaction
            t = t+tau;
            if mu == 1 %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                arnnumber = arnnumber + 1;
            elseif mu == 2 %mRNA decay
                arnnumber = arnnumber - 1;
                if(arnnumber<=0)
                    arnnumber=0;
                end
            elseif mu == 3 % translation : protein production
                protnumber=protnumber+1;
            elseif mu == 4 %protein decay
                protnumber= protnumber - 1;
                if(protnumber<=0)
                    protnumber=0;
                end
            end
            % update  molecule number, species and time for later output
            if t >= i*t_int
                i=i+1;
                t_array(i)=t;
                prom_array(i)=prom;
                arn_array(i)=arnnumber;
                prot_array(i)= protnumber;
                tau_array(i)=tau;
            end
        end
        
        exc_time= toc;
        fprintf('\nG2 simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        G2mean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        G2noise_mRNA_theoretical = sqrt(1/G2mean_mRNA_theoretical); %steady-state noise mRNA
        G2mean_protein_theoretical = (trans*trad)/(darn*dprot); %steady-state mean protein
        G2noise_protein_theoretical = sqrt((1/G2mean_protein_theoretical)+ (1/(G2mean_mRNA_theoretical*(1+(darn/dprot))))); %steady-state noise protein
        
        fprintf('\n-----------------G2 Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Protein)= %.2f\nnoise(Protein)= %.2f\n',  G2mean_mRNA_theoretical , G2noise_mRNA_theoretical,G2mean_protein_theoretical,G2noise_protein_theoretical);
        
        %%% simulation statistics
        G2mean_mRNA_simulation = mean(arn_array);
        G2stdev_mRNA_simulation = std(arn_array);
        G2noise_mRNA_simulation = G2stdev_mRNA_simulation/G2mean_mRNA_simulation;
        G2mean_protein_simulation = mean(prot_array);
        G2stdev_protein_simulation = std(prot_array);
        G2noise_protein_simulation = G2stdev_protein_simulation/G2mean_protein_simulation;
        
        fprintf('\n-----------------G2 Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Prot)= %.2f\nnoise(Prot)= %.2f\n',  G2mean_mRNA_simulation , G2noise_mRNA_simulation,G2mean_protein_simulation, G2noise_protein_simulation);
        
        
        %% M phase
        
        totG2 = i;
        t_end= time;  % End Time
        trans= m(1);
        arnnumber= round(rand *arnnumber);
        protnumber = round(rand*protnumber);
        trad=m(2);
        darn=m(3);
        dprot=m(4);
        prom=1; %Two promoters at the G2 phase
        %%% array to store results
        tic;
        while t< t_end
            
            %%% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            a= [trans*prom arnnumber*darn arnnumber*trad protnumber*dprot];
            a0= sum(a(:));
            
            %%% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
            while(r1==0 || r2==0), r1=rand; r2=rand; end
            tau= ((1/a0)* (log(1/r1)));
            
            j=1; mu=0; amu=0;
            
            while amu < r2*a0,
                mu = mu + 1;
                amu = amu + a(j);
                j = j + 1;
            end
            
            %%% Step 3 : Using the tau and r1 values obtained in step 2, increase
            % t by tau, and adjust the molecular population levels to reflect
            % the occurrence of one R reaction
            t = t+tau;
            if mu == 1 %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                arnnumber = arnnumber + 1;
            elseif mu == 2 %mRNA decay
                arnnumber = arnnumber - 1;
                if(arnnumber<=0)
                    arnnumber=0;
                end
            elseif mu == 3 % translation : protein production
                protnumber=protnumber+1;
            elseif mu == 4 %protein decay
                protnumber= protnumber - 1;
                if(protnumber<=0)
                    protnumber=0;
                end
            end
            % update  molecule number, species and time for later output
            if t >= i*t_int
                i=i+1;
                t_array(i)=t;
                prom_array(i)=prom;
                arn_array(i)=arnnumber;
                prot_array(i)= protnumber;
                tau_array(i)=tau;
            end
        end
        
        exc_time= toc;
        fprintf('\nM simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        Mmean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        Mnoise_mRNA_theoretical = sqrt(1/Mmean_mRNA_theoretical); %steady-state noise mRNA
        Mmean_protein_theoretical = (trans*trad)/(darn*dprot); %steady-state mean protein
        Mnoise_protein_theoretical = sqrt((1/Mmean_protein_theoretical)+ (1/(Mmean_mRNA_theoretical*(1+(darn/dprot))))); %steady-state noise protein
        
        fprintf('\n-----------------M Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Protein)= %.2f\nnoise(Protein)= %.2f\n',  Mmean_mRNA_theoretical , Mnoise_mRNA_theoretical,Mmean_protein_theoretical,Mnoise_protein_theoretical);
        
        %%% simulation statistics
        Mmean_mRNA_simulation = mean(arn_array);
        Mstdev_mRNA_simulation = std(arn_array);
        Mnoise_mRNA_simulation = Mstdev_mRNA_simulation/Mmean_mRNA_simulation;
        Mmean_protein_simulation = mean(prot_array);
        Mstdev_protein_simulation = std(prot_array);
        Mnoise_protein_simulation = Mstdev_protein_simulation/Mmean_protein_simulation;
        
        fprintf('\n-----------------M Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Prot)= %.2f\nnoise(Prot)= %.2f\n',  Mmean_mRNA_simulation , Mnoise_mRNA_simulation,Mmean_protein_simulation, Mnoise_protein_simulation);
        
        simul.time{n}=t_array;
        simul.prot{n} = prot_array;
        simul.arn{n}= arn_array;
        simul.prom{n}=prom_array;
        simul.tau{n}=tau_array;
        
    end
    
     taillemin=  min(cellfun('length',simul.time));
    
    for i=1:ncycle
        simul.time{i} = simul.time{i}(1:taillemin);
        simul.arn{i} = simul.arn{i}(1:taillemin);
        simul.tau{i} = simul.tau{i}(1:taillemin);
        simul.prot{i} = simul.prot{i}(1:taillemin);
        
    end
    t_array= mean(reshape(cell2mat(simul.time), taillemin, ncycle),2);
    arn_array= round(mean(reshape(cell2mat(simul.arn), taillemin, ncycle),2));
    tau_array= mean(reshape(cell2mat(simul.tau), taillemin, ncycle),2);
    prot_array= round(mean(reshape(cell2mat(simul.prot),taillemin, ncycle),2));

    
  

    %% Plotting result in a fancy way
    hold on;
    
    %%% ARNm output
    handle{1}= subplot(2,2,1);
    plot(t_array,arn_array); %mRNA time series total
    xlabel('time (s)');
    ylabel('mRNA no.');
    
    handle{2}= subplot(2,2,2);
    x= min(arn_array):max(arn_array);
    liste.arn= [0:max(arn_array); hist(arn_array, max(arn_array)+1)]';
    bar(x, hist(arn_array, max(arn_array)+1-min(arn_array))./ sum(hist(arn_array)), 'hist');%mRNA histogram
    legend(sprintf(' mean = %g \n stdev = %g',mean(arn_array),std(arn_array)))
    set(gca,'XTick',min(arn_array):round((max(arn_array)-min(arn_array))/2):max(arn_array));
    %%% Fitting a poisson distribution
    prob_distribution = fitdist(arn_array(:),'Poisson');
    hold on;
    probdf = pdf(prob_distribution,x);
    plot(x,probdf,'LineWidth',2, 'color', 'r');
    xlabel('mRNA no.');
    ylabel('x100 (%)');
    fprintf('\n-----------------Fitting poisson distribution-----------------\n');
    [H, P, STATS] = chi2gof(arn_array,'cdf', @(z)poisscdf(z,mean(arn_array)), 'nparams', 1, 'nbins', max(arn_array)+1-min(arn_array))
    
    %%% Protein output
    handle{3}= subplot(2,2,3);
    plot(t_array,prot_array);
    xlabel('time (s)');
    ylabel('protein no.')
    handle{4}=subplot(2,2,4);
    hist(prot_array, max(prot_array)-min(prot_array)+1); %protein histogram
    xlabel('protein no.');
    ylabel('occurrence');
    legend(sprintf(' mean = %g \n stdev = %g',mean(prot_array),std(prot_array)))
    liste.prot= [0:max(prot_array); hist(prot_array, max(prot_array)+1)]';
    save('data.mat', 't_array', 'prom_array', 'arn_array', 'prot_array', 'tau_array', 'simul', 'totG1', 'totG2');
    
else
    for n=1:ncycle
        
        %% Phase 1 : G1/S
        %%%  Thattai - van Oudenaarden model (single gene model)
        %%% reaction prom -> ARN + prom (k1)
        %%% reaction ARN -> 0 (k2)
        arnnumber=arnnmber; protnumber=protnmber;

        t=0;
        t_end= g1(end);  % End Time
        trans= g1(1);
        darn=g1(3);
        prom=1; % Only one promoters, and it's always actif
        %%% array to store results
        i=1; %indice of array
        t_array(i)=t;
        tau_array(i)=0;
        prom_array(i)=prom;
        arn_array(i)=arnnumber;
        tic;
        
        while t< t_end
            
            %%% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            a= [trans*prom arnnumber*darn];
            a0= sum(a(:));
            
            %%% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            
            r1=rand; r2=rand;
            while(r1==0 || r2==0), r1=rand; r2=rand; end
            tau= ((1/a0)* (log(1/r1)));
            
            j=1; mu=0; amu=0;
            
            while amu < r2*a0,
                mu = mu + 1;
                amu = amu + a(j);
                j = j + 1;
            end
            
            %% Step 3 : Using the tau and r1 values obtained in step 2, increase
            % t by tau, and adjust the molecular population levels to reflect
            % the occurrence of one R reaction
            t = t+tau;
            if mu == 1 %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                arnnumber = arnnumber + 1;
            elseif mu == 2 %mRNA decay
                arnnumber = arnnumber - 1;
                if(arnnumber<=0)
                    arnnumber=0;
                end
            end
            % update  molecule number, species and time for later
            % output each 1s
            if t >= i*t_int
                i=i+1;
                t_array(i)=t;
                prom_array(i)=prom;
                arn_array(i)=arnnumber;
                tau_array(i)=tau;
            end
        end
        
        exc_time= toc;
        fprintf('G1/S simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        G1mean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        G1noise_mRNA_theoretical = sqrt(1/G1mean_mRNA_theoretical); %steady-state noise mRNA
        fprintf('\n-----------------G1/S Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  G1mean_mRNA_theoretical , G1noise_mRNA_theoretical);
        
        %%% simulation statistics
        G1mean_mRNA_simulation = mean(arn_array);
        G1stdev_mRNA_simulation = std(arn_array);
        G1noise_mRNA_simulation = G1stdev_mRNA_simulation/G1mean_mRNA_simulation;
        fprintf('\n-----------------G1/S Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  G1mean_mRNA_simulation , G1noise_mRNA_simulation);
        
        
        %% G2 phase
        
        totG1 = i;
        t_end= g1(end)+g2(end);  % End Time
        trans= g2(1);
        darn=g2(3);
        prom=1; %Two promoters at the G2 phase
        %%% array to store results
        tic;
        while t< t_end
            %% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            a= [trans*prom arnnumber*darn];
            a0= sum(a(:));
            
            %% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
            while(r1==0 || r2==0), r1=rand; r2=rand; end
            tau= ((1/a0)* (log(1/r1)));
            
            j=1; mu=0; amu=0;
            
            while amu < r2*a0,
                mu = mu + 1;
                amu = amu + a(j);
                j = j + 1;
            end
            
            %% Step 3 : Using the tau and r1 values obtained in step 2, increase
            % t by tau, and adjust the molecular population levels to reflect
            % the occurrence of one R reaction
            t = t+tau;
            if mu == 1 %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                arnnumber = arnnumber + 1;
            elseif mu == 2 %mRNA decay
                arnnumber = arnnumber - 1;
                if(arnnumber<=0)
                    arnnumber=0;
                end
            end
            % update  molecule number, species and time for later output
            if t >= i*t_int
                i=i+1;
                t_array(i)=t;
                prom_array(i)=prom;
                arn_array(i)=arnnumber;
                tau_array(i)=tau;
            end
        end
        
        exc_time= toc;
        fprintf('\nG2 simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        G2mean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        G2noise_mRNA_theoretical = sqrt(1/G2mean_mRNA_theoretical); %steady-state noise mRNA
        fprintf('\n-----------------G2 Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  G2mean_mRNA_theoretical , G2noise_mRNA_theoretical);
        
        %%% simulation statistics
        G2mean_mRNA_simulation = mean(arn_array);
        G2stdev_mRNA_simulation = std(arn_array);
        G2noise_mRNA_simulation = G2stdev_mRNA_simulation/G2mean_mRNA_simulation;
        fprintf('\n-----------------G2 Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  G2mean_mRNA_simulation , G2noise_mRNA_simulation);
        
        
        %% M phase
        
        totG2 = i;
        t_end= time;  % End Time
        trans= m(1);
        arnnumber= round(rand *arnnumber);
        darn=m(3);
        prom=1; %Two promoters at the G2 phase
        %%% array to store results
        tic;
        while t< t_end
            
            %% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            a= [trans*prom arnnumber*darn];
            a0= sum(a(:));
            
            %% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
            while(r1==0 || r2==0), r1=rand; r2=rand; end
            tau= ((1/a0)* (log(1/r1)));
            
            j=1; mu=0; amu=0;
            
            while amu < r2*a0,
                mu = mu + 1;
                amu = amu + a(j);
                j = j + 1;
            end
            
            %% Step 3 : Using the tau and r1 values obtained in step 2, increase
            % t by tau, and adjust the molecular population levels to reflect
            % the occurrence of one R reaction
            t = t+tau;
            if mu == 1 %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                arnnumber = arnnumber + 1;
            elseif mu == 2 %mRNA decay
                arnnumber = arnnumber - 1;
                if(arnnumber<=0)
                    arnnumber=0;
                end
            end
            % update  molecule number, species and time for later output
            if t >= i*t_int
                i=i+1;
                t_array(i)=t;
                prom_array(i)=prom;
                arn_array(i)=arnnumber;
                tau_array(i)=tau;
            end
        end
        
        exc_time= toc;
        fprintf('\nM simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        Mmean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        Mnoise_mRNA_theoretical = sqrt(1/Mmean_mRNA_theoretical); %steady-state noise mRNA
        fprintf('\n-----------------M Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  Mmean_mRNA_theoretical , Mnoise_mRNA_theoretical);
        
        %%% simulation statistics
        Mmean_mRNA_simulation = mean(arn_array);
        Mstdev_mRNA_simulation = std(arn_array);
        Mnoise_mRNA_simulation = Mstdev_mRNA_simulation/Mmean_mRNA_simulation;
        fprintf('\n-----------------M Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  Mmean_mRNA_simulation , Mnoise_mRNA_simulation);
        
        simul.time{n}=t_array;
        simul.arn{n}= arn_array;
        simul.prom{n}=prom_array;
        simul.tau{n}=tau_array;
    end
    
    %% Make it to the right size in order to calculate min
    taillemin=  min(cellfun('length',simul.time));
    
    for i=1:ncycle
        simul.time{i} = simul.time{i}(1:taillemin);
        simul.arn{i} = simul.arn{i}(1:taillemin);
        simul.tau{i} = simul.tau{i}(1:taillemin);
        
    end
    t_array= mean(reshape(cell2mat(simul.time), taillemin, ncycle),2);
    arn_array= round(mean(reshape(cell2mat(simul.arn), taillemin, ncycle),2));
    tau_array= mean(reshape(cell2mat(simul.tau), taillemin, ncycle),2);
    
    %% Plotting result in a fancy way
    hold on;
    
    %%% ARNm output
    handle{1}= subplot(2,1,1);
    plot(t_array,arn_array); %mRNA time series total
    xlabel('time (s)');
    ylabel('mRNA no.');
    
    handle{2}= subplot(2,1,2);
    x= min(arn_array):max(arn_array);
    liste.arn= [0:max(arn_array); hist(arn_array, max(arn_array)+1)]';
    bar(x, hist(arn_array, max(arn_array)+1-min(arn_array))./ sum(hist(arn_array)), 'hist');%mRNA histogram
    legend(sprintf(' mean = %g \n stdev = %g',mean(arn_array),std(arn_array)))
    set(gca,'XTick',min(arn_array):round(log2(max(arn_array))):max(arn_array));
    %%% Fitting a poisson distribution
    prob_distribution = fitdist(arn_array(:),'Poisson');
    hold on;
    probdf = pdf(prob_distribution,x);
    plot(x,probdf,'LineWidth',2, 'color', 'r');
    xlabel('mRNA no.');
    ylabel('x100 (%)');
    fprintf('\n-----------------Fitting poisson distribution-----------------\n');
    [H, P, STATS] = chi2gof(arn_array,'cdf', @(z)poisscdf(z,mean(arn_array)), 'nparams', 1, 'nbins', max(arn_array)+1-min(arn_array))
    
    save('data.mat', 't_array', 'prom_array', 'arn_array', 'tau_array', 'simul', 'totG1', 'totG2');
    
end
end


function [liste, handle]= burst(g1, g2, m, arnnumber, protnumber, time, ncycle)

t_int=1;
simul.time={};
simul.prot={};
simul.arn={};
simul.promnumber={};
simul.promstate={};
simul.tau={};

if((g1(2)~=0 || g1(4) ~=0) && (g2(2)~=0 || (g2(4)~=0)) && (m(2)~=0 || m(4)~=0))
    
    for n=1:ncycle
        
        %% Phase 1 : G1/S
        %%%  Thattai - van Oudenaarden model (single gene model)
        %%% reaction prom -> ARN + prom (k1)
        %%% reaction ARN -> 0 (k2)
        arnnumber=arnnmber; protnumber=protnmber;

        t=0;
        t_end= g1(end);  % End Time
        trans= g1(1);
        trad=g1(2);
        darn=g1(3);
        dprot=g1(4);
        ig=g1(5);
        ag=g1(6);
        prom.number=1; % Only one promoters, and it's always actif
        prom.state=1;
        %%% array to store results
        i=1; %indice of array
        t_array(i)=t;
        tau_array(i)=0;
        prom_array(i)=prom.number;
        promstat_array(i) = prom.state;
        arn_array(i)=arnnumber;
        prot_array(i)=protnumber; %protein array and initial value
        tic;
        
        while t< t_end
            
            %%% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            if(prom.state==1)
                a= [prom.number*ig trans*prom.number arnnumber*darn arnnumber*trad protnumber*dprot];
            elseif(prom.state==0)
                a= [prom.number*ag arnnumber*darn arnnumber*trad protnumber*dprot];
            end
            a0= sum(a(:));
            
            %%% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
        while(r1==0 || r2==0), r1=rand; r2=rand; end
        tau= ((1/a0)* (log(1/r1)));
        
        j=1; mu=0; amu=0;
        
        while amu < r2*a0,
            mu = mu + 1;
            amu = amu + a(j);
            j = j + 1;
        end
        
        %% Step 3 : Using the tau and r1 values obtained in step 2, increase
        % t by tau, and adjust the molecular population levels to reflect
        % the occurrence of one R reaction
        t = t+tau;
        switch prom.state
            case 0
                if mu==1 %prom activation
                    prom.state=1;
                elseif mu == 2 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                elseif mu==3 %traduction
                    protnumber=protnumber+1;
                elseif mu == 4 % prot decay
                    protnumber = protnumber - 1;
                    if(protnumber<=0)
                        protnumber=0;
                    end
                end
                
            case 1
                
                if mu ==1   % prom inactivation
                    prom.state=0;
                elseif mu==2 && prom.state==1  %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                    arnnumber = arnnumber + 1;
                elseif mu == 3 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                elseif mu==4 %traduction
                    protnumber = protnumber+1;
                elseif mu==5
                    protnumber = protnumber - 1;
                    if(protnumber<=0)
                        protnumber=0;
                    end
                end
        end
        % update  molecule number, species and time for later output
        if t >= i*t_int
            i=i+1;
            t_array(i)=t;
            prom_array(i)=prom.number;
            promstat_array(i)=prom.state;
            arn_array(i)=arnnumber;
            prot_array(i)= protnumber;
            tau_array(i)=tau;
        end
        end
        
        exc_time= toc;
        fprintf('G1/S simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        G1mean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        G1noise_mRNA_theoretical = sqrt(1/G1mean_mRNA_theoretical); %steady-state noise mRNA
        G1mean_protein_theoretical = (trans*trad)/(darn*dprot); %steady-state mean protein
        G1noise_protein_theoretical = sqrt((1/G1mean_protein_theoretical)+ (1/(G1mean_mRNA_theoretical*(1+(darn/dprot))))); %steady-state noise protein
        
        fprintf('\n-----------------G1/S Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Protein)= %.2f\nnoise(Protein)= %.2f\n',  G1mean_mRNA_theoretical , G1noise_mRNA_theoretical,G1mean_protein_theoretical,G1noise_protein_theoretical);
        
        %%% simulation statistics
        G1mean_mRNA_simulation = mean(arn_array);
        G1stdev_mRNA_simulation = std(arn_array);
        G1noise_mRNA_simulation = G1stdev_mRNA_simulation/G1mean_mRNA_simulation;
        G1mean_protein_simulation = mean(prot_array);
        G1stdev_protein_simulation = std(prot_array);
        G1noise_protein_simulation = G1stdev_protein_simulation/G1mean_protein_simulation;
        
        fprintf('\n-----------------G1/S Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Prot)= %.2f\nnoise(Prot)= %.2f\n',  G1mean_mRNA_simulation , G1noise_mRNA_simulation,G1mean_protein_simulation, G1noise_protein_simulation);
        
        
        %% G2 phase
        
        totG1 = i;
        t_end= g1(end)+g2(end);  % End Time
        trans= g2(1);
        trad=g2(2);
        darn=g2(3);
        dprot=g2(4);
        ig=g2(5);
        ag=g2(6);
        prom.number=2; % Only one promoters, and it's always actif       
        tic;
        while t< t_end
            
            %%% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            if(prom.state==1)
                a= [prom.number*ig trans*prom.number arnnumber*darn arnnumber*trad protnumber*dprot];
            elseif(prom.state==0)
                a= [prom.number*ag arnnumber*darn arnnumber*trad protnumber*dprot];
            end
            a0= sum(a(:));
            
            %%% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
        while(r1==0 || r2==0), r1=rand; r2=rand; end
        tau= ((1/a0)* (log(1/r1)));
        
        j=1; mu=0; amu=0;
        
        while amu < r2*a0,
            mu = mu + 1;
            amu = amu + a(j);
            j = j + 1;
        end
        
        %% Step 3 : Using the tau and r1 values obtained in step 2, increase
        % t by tau, and adjust the molecular population levels to reflect
        % the occurrence of one R reaction
        t = t+tau;
        switch prom.state
            case 0
                if mu==1 %prom activation
                    prom.state=1;
                elseif mu == 2 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                elseif mu==3 %traduction
                    protnumber=protnumber+1;
                elseif mu == 4 % prot decay
                    protnumber = protnumber - 1;
                    if(protnumber<=0)
                        protnumber=0;
                    end
                end
                
            case 1
                
                if mu ==1   % prom inactivation
                    prom.state=0;
                elseif mu==2 && prom.state==1  %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                    arnnumber = arnnumber + 1;
                elseif mu == 3 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                elseif mu==4 %traduction
                    protnumber = protnumber+1;
                elseif mu==5
                    protnumber = protnumber - 1;
                    if(protnumber<=0)
                        protnumber=0;
                    end
                end
        end
        % update  molecule number, species and time for later output
        if t >= i*t_int
            i=i+1;
            t_array(i)=t;
            prom_array(i)=prom.number;
            promstat_array(i)=prom.state;
            arn_array(i)=arnnumber;
            prot_array(i)= protnumber;
            tau_array(i)=tau;
        end
        end
        
        exc_time= toc;
        fprintf('\nG2 simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        G2mean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        G2noise_mRNA_theoretical = sqrt(1/G2mean_mRNA_theoretical); %steady-state noise mRNA
        G2mean_protein_theoretical = (trans*trad)/(darn*dprot); %steady-state mean protein
        G2noise_protein_theoretical = sqrt((1/G2mean_protein_theoretical)+ (1/(G2mean_mRNA_theoretical*(1+(darn/dprot))))); %steady-state noise protein
        
        fprintf('\n-----------------G2 Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Protein)= %.2f\nnoise(Protein)= %.2f\n',  G2mean_mRNA_theoretical , G2noise_mRNA_theoretical,G2mean_protein_theoretical,G2noise_protein_theoretical);
        
        %%% simulation statistics
        G2mean_mRNA_simulation = mean(arn_array);
        G2stdev_mRNA_simulation = std(arn_array);
        G2noise_mRNA_simulation = G2stdev_mRNA_simulation/G2mean_mRNA_simulation;
        G2mean_protein_simulation = mean(prot_array);
        G2stdev_protein_simulation = std(prot_array);
        G2noise_protein_simulation = G2stdev_protein_simulation/G2mean_protein_simulation;
        
        fprintf('\n-----------------G2 Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Prot)= %.2f\nnoise(Prot)= %.2f\n',  G2mean_mRNA_simulation , G2noise_mRNA_simulation,G2mean_protein_simulation, G2noise_protein_simulation);
        
        
        %% M phase
        
        totG2 = i;
        t_end= time;  % End Time
        trans= m(1);
        arnnumber= round(rand *arnnumber);
        protnumber = round(rand*protnumber);
        trad=m(2);
        darn=m(3);
        dprot=m(4);
        ig=m(5);
        ag=m(6);
        prom.number=1; %Two promoters at the G2 phase
        %%% array to store results
        tic;
        while t< t_end
            
            %%% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            if(prom.state==1)
                a= [prom.number*ig trans*prom.number arnnumber*darn arnnumber*trad protnumber*dprot];
            elseif(prom.state==0)
                a= [prom.number*ag arnnumber*darn arnnumber*trad protnumber*dprot];
            end
            a0= sum(a(:));
            
            %%% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
        while(r1==0 || r2==0), r1=rand; r2=rand; end
        tau= ((1/a0)* (log(1/r1)));
        
        j=1; mu=0; amu=0;
        
        while amu < r2*a0,
            mu = mu + 1;
            amu = amu + a(j);
            j = j + 1;
        end
        
        %% Step 3 : Using the tau and r1 values obtained in step 2, increase
        % t by tau, and adjust the molecular population levels to reflect
        % the occurrence of one R reaction
        t = t+tau;
        switch prom.state
            case 0
                if mu==1 %prom activation
                    prom.state=1;
                elseif mu == 2 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                elseif mu==3 %traduction
                    protnumber=protnumber+1;
                elseif mu == 4 % prot decay
                    protnumber = protnumber - 1;
                    if(protnumber<=0)
                        protnumber=0;
                    end
                end
                
            case 1
                
                if mu ==1   % prom inactivation
                    prom.state=0;
                elseif mu==2 && prom.state==1  %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                    arnnumber = arnnumber + 1;
                elseif mu == 3 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                elseif mu==4 %traduction
                    protnumber = protnumber+1;
                elseif mu==5
                    protnumber = protnumber - 1;
                    if(protnumber<=0)
                        protnumber=0;
                    end
                end
        end
        % update  molecule number, species and time for later output
        
        if t >= i*t_int
            i=i+1;
            t_array(i)=t;
            prom_array(i)=prom.number;
            promstat_array(i)=prom.state;
            arn_array(i)=arnnumber;
            prot_array(i)= protnumber;
            tau_array(i)=tau;
        end
        end
        
        exc_time= toc;
        fprintf('\nM simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        Mmean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        Mnoise_mRNA_theoretical = sqrt(1/Mmean_mRNA_theoretical); %steady-state noise mRNA
        Mmean_protein_theoretical = (trans*trad)/(darn*dprot); %steady-state mean protein
        Mnoise_protein_theoretical = sqrt((1/Mmean_protein_theoretical)+ (1/(Mmean_mRNA_theoretical*(1+(darn/dprot))))); %steady-state noise protein
        
        fprintf('\n-----------------M Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Protein)= %.2f\nnoise(Protein)= %.2f\n',  Mmean_mRNA_theoretical , Mnoise_mRNA_theoretical,Mmean_protein_theoretical,Mnoise_protein_theoretical);
        
        %%% simulation statistics
        Mmean_mRNA_simulation = mean(arn_array);
        Mstdev_mRNA_simulation = std(arn_array);
        Mnoise_mRNA_simulation = Mstdev_mRNA_simulation/Mmean_mRNA_simulation;
        Mmean_protein_simulation = mean(prot_array);
        Mstdev_protein_simulation = std(prot_array);
        Mnoise_protein_simulation = Mstdev_protein_simulation/Mmean_protein_simulation;
        
        fprintf('\n-----------------M Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Prot)= %.2f\nnoise(Prot)= %.2f\n',  Mmean_mRNA_simulation , Mnoise_mRNA_simulation,Mmean_protein_simulation, Mnoise_protein_simulation);
        
        simul.time{n}=t_array;
        simul.prot{n} = prot_array;
        simul.arn{n}= arn_array;
        simul.promnumber{n}=prom_array;
        simul.promstate{n}=promstat_array;
        simul.tau{n}=tau_array;
    end
    
    %% Plotting result in a fancy way
    cla;
    hold on;
    handle{1}= subplot(3,2,1);
    plot(t_array,arn_array); %mRNA time series
    set(gca, 'Ylim', [0 max(arn_array)+1]);
    xlabel('time (s)');
    ylabel('mRNA no.');
    handle{2}= subplot(3,2,2);
    %         length(arn_array), max(arn_array)
    %         for j=0:max(arn_array)
    %             sum(arn_array==j)/ length(arn_array)
    %         end
    %         hist(arn_array, max(arn_array)); %mRNA histogram
    
    liste.arn= [0:max(arn_array); hist(arn_array, max(arn_array)+1)]';
    hist(arn_array, max(arn_array)-min(arn_array)+1);
    set(gca, 'XLim', [min(arn_array) max(arn_array)+1]);
    set(gca,'XTick',min(arn_array):round((max(arn_array)-min(arn_array))/2):max(arn_array));
    legend(sprintf(' mean = %g \n stdev = %g',mean(arn_array),std(arn_array)))
    
    xlabel('mRNA no.');
    ylabel('count');
    
    handle{3}= subplot(3,2,[5 6]);
    t=1;
    for i=1:numel(promstat_array)-1
        timearray(t)=t_array(i);
        timearray(t+1)=t_array(i);
        promoarray(t)=promstat_array(i);
        if(prom_array(i)==promstat_array(i+1))
            promoarray(t+1)=promstat_array(i);
        else
            promoarray(t+1)=promstat_array(i+1);
        end
        t=t+2;
    end
    plot(timearray, promoarray, '-r*');
    
    ylabel('Promoter state.');
    xlabel('time (s)');
    set(gca, 'YLim', [-1 2]);
    set(gca, 'YTick', [0 1]);
    set(gca,'YTickLabel',{'off', 'on'});
    
%     fprintf('\n-----------------Fitting poisson distribution-----------------\n');
%     [H P STATS] = chi2gof(arn_array,'cdf', @(z)poisscdf(z,mean(arn_array)), 'nparams', 1, 'nbins', max(arn_array)+1-min(arn_array))
        
    %%% Protein output
    handle{3}= subplot(3,2,3);
    plot(t_array,prot_array);
    xlabel('time (s)');
    ylabel('protein no.');
    handle{4}=subplot(3,2,4);
    hist(prot_array, max(prot_array)-min(prot_array)+1); %protein histogram
    xlabel('protein no.');
    ylabel('occurrence');
    legend(sprintf(' mean = %g \n stdev = %g',mean(prot_array),std(prot_array)))
    liste.prot= [0:max(prot_array); hist(prot_array, max(prot_array)+1)]'; 
    save('data.mat', 't_array', 'prom_array', 'arn_array', 'tau_array', 'simul', 'totG1', 'totG2');
    
else
    for n=1:ncycle
        
        %% Phase 1 : G1/S
        %%%  Thattai - van Oudenaarden model (single gene model)
        %%% reaction prom -> ARN + prom (k1)
        %%% reaction ARN -> 0 (k2)
        arnnumber=arnnmber; protnumber=protnmber;
        t=0;
        t_end= g1(end);  % End Time
        trans= g1(1);
        darn=g1(3);
        ig=g1(5);
        ag=g1(6);
        prom.number=1; % Only one promoters, and it's always actif
        prom.state=1;
        %%% array to store results
        i=1; %indice of array
        t_array(i)=t;
        tau_array(i)=0;
        prom_array(i)=prom.number;
        promstat_array=prom.state;
        arn_array(i)=arnnumber;
        tic;
        
        while t< t_end
            
            %%% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            if(prom.state==1)
                a= [prom.number*ig trans*prom.number arnnumber*darn];
            elseif(prom.state==0)
                a= [prom.number*ag arnnumber*darn];
            end
            a0= sum(a(:));
            
            %%% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
        while(r1==0 || r2==0), r1=rand; r2=rand; end
        tau= ((1/a0)* (log(1/r1)));
        
        j=1; mu=0; amu=0;
        
        while amu < r2*a0,
            mu = mu + 1;
            amu = amu + a(j);
            j = j + 1;
        end
        
        %% Step 3 : Using the tau and r1 values obtained in step 2, increase
        % t by tau, and adjust the molecular population levels to reflect
        % the occurrence of one R reaction
        t = t+tau;
        switch prom.state
            case 0
                if mu==1 %prom activation
                    prom.state=1;
                elseif mu == 2 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                end
                
            case 1
                
                if mu ==1   % prom inactivation
                    prom.state=0;
                elseif mu==2 && prom.state==1  %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                    arnnumber = arnnumber + 1;
                elseif mu == 3 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                end
        end
            % update  molecule number, species and time for later output
         if t >= i*t_int
            i=i+1;
            t_array(i)=t;
            prom_array(i)=prom.number;
            promstat_array(i)=prom.state;
            arn_array(i)=arnnumber;
            tau_array(i)=tau;
            
         end
        end
        
        exc_time= toc;
        fprintf('G1/S simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        G1mean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        G1noise_mRNA_theoretical = sqrt(1/G1mean_mRNA_theoretical); %steady-state noise mRNA
        fprintf('\n-----------------G1/S Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  G1mean_mRNA_theoretical , G1noise_mRNA_theoretical);
        
        %%% simulation statistics
        G1mean_mRNA_simulation = mean(arn_array);
        G1stdev_mRNA_simulation = std(arn_array);
        G1noise_mRNA_simulation = G1stdev_mRNA_simulation/G1mean_mRNA_simulation;
        fprintf('\n-----------------G1/S Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  G1mean_mRNA_simulation , G1noise_mRNA_simulation);
        
        
        %% G2 phase
        
        totG1 = i;
        t_end= g1(end)+g2(end);  % End Time
        trans= g2(1);
        darn=g2(3);
        ig=g2(5);
        ag=g2(6);
        prom.number=2; %Two promoters at the G2 phase
        %%% array to store results
        tic;
        while t< t_end
            %% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            if(prom.state==1)
                a= [prom.number*ig trans*prom.number arnnumber*darn];
            elseif(prom.state==0)
                a= [prom.number*ag arnnumber*darn];
            end
            a0= sum(a(:));
            %% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
            while(r1==0 || r2==0), r1=rand; r2=rand; end
            tau= ((1/a0)* (log(1/r1)));
            
            j=1; mu=0; amu=0;
            
            while amu < r2*a0,
                mu = mu + 1;
                amu = amu + a(j);
                j = j + 1;
            end
            
            %% Step 3 : Using the tau and r1 values obtained in step 2, increase
            % t by tau, and adjust the molecular population levels to reflect
            % the occurrence of one R reaction
            t = t+tau;
              t = t+tau;
        switch prom.state
            case 0
                if mu==1 %prom activation
                    prom.state=1;
                elseif mu == 2 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                end
                
            case 1
                
                if mu ==1   % prom inactivation
                    prom.state=0;
                elseif mu==2 && prom.state==1  %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                    arnnumber = arnnumber + 1;
                elseif mu == 3 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                end
        end
            % update  molecule number, species and time for later output
          if t >= i*t_int
            i=i+1;
            t_array(i)=t;
            prom_array(i)=prom.number;
            promstat_array(i)=prom.state;
            arn_array(i)=arnnumber;
            tau_array(i)=tau;
          end 
        end
        
        exc_time= toc;
        fprintf('\nG2 simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        G2mean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        G2noise_mRNA_theoretical = sqrt(1/G2mean_mRNA_theoretical); %steady-state noise mRNA
        fprintf('\n-----------------G2 Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  G2mean_mRNA_theoretical , G2noise_mRNA_theoretical);
        
        %%% simulation statistics
        G2mean_mRNA_simulation = mean(arn_array);
        G2stdev_mRNA_simulation = std(arn_array);
        G2noise_mRNA_simulation = G2stdev_mRNA_simulation/G2mean_mRNA_simulation;
        fprintf('\n-----------------G2 Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  G2mean_mRNA_simulation , G2noise_mRNA_simulation);
        
        
        %% M phase
        
        totG2 = i;
        t_end= time;  % End Time
        trans= m(1);
        arnnumber= round(rand *arnnumber);
        darn=m(3);
        ig=m(5);
        ag=m(6);
        prom.number=1; %Two promoters at the G2 phase
        %%% array to store results
        tic;
        while t< t_end
            
            %% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
            % am = hm*cm for the current molecular population numbers.
            if(prom.state==1)
                a= [prom.number*ig trans*prom.number arnnumber*darn];
            elseif(prom.state==0)
                a= [prom.number*ag arnnumber*darn];
            end
            a0= sum(a(:));
            
            %% Step 2 : Generate two random numbers r1 and r2 using the
            % unit-interval uniform random number generator, and calculate tau
            % and mu
            
            r1=rand; r2=rand;
            while(r1==0 || r2==0), r1=rand; r2=rand; end
            tau= ((1/a0)* (log(1/r1)));
            
            j=1; mu=0; amu=0;
            
            while amu < r2*a0,
                mu = mu + 1;
                amu = amu + a(j);
                j = j + 1;
            end
            
            %% Step 3 : Using the tau and r1 values obtained in step 2, increase
            % t by tau, and adjust the molecular population levels to reflect
            % the occurrence of one R reaction
            t = t+tau;
            switch prom.state
            case 0
                if mu==1 %prom activation
                    prom.state=1;
                elseif mu == 2 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                end
                
            case 1
                
                if mu ==1   % prom inactivation
                    prom.state=0;
                elseif mu==2 && prom.state==1  %transcription, this can be changed if we consider simultaneous multi trans for the same prom
                    arnnumber = arnnumber + 1;
                elseif mu == 3 %mRNA decay
                    arnnumber = arnnumber - 1;
                    if(arnnumber<=0)
                        arnnumber=0;
                    end
                end
        end
            % update  molecule number, species and time for later output
           if t >= i*t_int
            i=i+1;
            t_array(i)=t;
            prom_array(i)=prom.number;
            promstat_array(i)=prom.state;
            arn_array(i)=arnnumber;
            tau_array(i)=tau;
            
           end
        end
        
        exc_time= toc;
        fprintf('\nM simulation time: %f seconds\n', exc_time);
        
        %%% analytical results for the model
        Mmean_mRNA_theoretical = trans/darn; %steady-state mean mRNA
        Mnoise_mRNA_theoretical = sqrt(1/Mmean_mRNA_theoretical); %steady-state noise mRNA
        fprintf('\n-----------------M Analytical result-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  Mmean_mRNA_theoretical , Mnoise_mRNA_theoretical);
        
        %%% simulation statistics
        Mmean_mRNA_simulation = mean(arn_array);
        Mstdev_mRNA_simulation = std(arn_array);
        Mnoise_mRNA_simulation = Mstdev_mRNA_simulation/Mmean_mRNA_simulation;
        fprintf('\n-----------------M Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  Mmean_mRNA_simulation , Mnoise_mRNA_simulation);
        
        simul.time{n}=t_array;
        simul.arn{n}= arn_array;
        simul.promnumber{n}=prom_array;
        simul.promstate{n}=promstat_array;
        simul.tau{n}=tau_array;
    end
     %% Plotting result in a fancy way
    cla;
    hold on;
    handle{1}= subplot(2,2,1);
    plot(t_array,arn_array); %mRNA time series
    set(gca, 'Ylim', [0 max(arn_array)+1]);
    xlabel('time (s)');
    ylabel('mRNA no.');
    handle{2}= subplot(2,2,2);
    %         length(arn_array), max(arn_array)
    %         for j=0:max(arn_array)
    %             sum(arn_array==j)/ length(arn_array)
    %         end
    %         hist(arn_array, max(arn_array)); %mRNA histogram
    
    liste.arn= [0:max(arn_array); hist(arn_array, max(arn_array)+1)]';
    hist(arn_array, max(arn_array)-min(arn_array)+1);
    set(gca, 'XLim', [min(arn_array) max(arn_array)+1]);
    set(gca,'XTick',min(arn_array):round((max(arn_array)-min(arn_array))/2):max(arn_array));
    legend(sprintf(' mean = %g \n stdev = %g',mean(arn_array),std(arn_array)))
    
    xlabel('mRNA no.');
    ylabel('count');
    
    handle{3}= subplot(2,2,[3 4]);
    t=1;
    for i=1:numel(promstat_array)-1
        timearray(t)=t_array(i);
        timearray(t+1)=t_array(i);
        promoarray(t)=promstat_array(i);
        if(prom_array(i)==promstat_array(i+1))
            promoarray(t+1)=promstat_array(i);
        else
            promoarray(t+1)=promstat_array(i+1);
        end
        t=t+2;
    end
    plot(timearray, promoarray, '-r*');
    
    ylabel('Promoter state.');
    xlabel('time (s)');
    set(gca, 'YLim', [-1 2]);
    set(gca, 'YTick', [0 1]);
    set(gca,'YTickLabel',{'off', 'on'});
    hold off;
%     
%     fprintf('\n-----------------Fitting poisson distribution-----------------\n');
%     [H P STATS] = chi2gof(arn_array,'cdf', @(z)poisscdf(z,mean(arn_array)), 'nparams', 1, 'nbins', max(arn_array)+1-min(arn_array))
    save('data.mat', 't_array', 'prom_array', 'arn_array', 'tau_array', 'simul', 'totG1', 'totG2');
    
end
end