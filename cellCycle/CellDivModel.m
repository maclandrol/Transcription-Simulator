function [liste, handle]=CellDivModel(mode, time, trans, darn ,arnnumber,  trad, dprot, protnumber,ag , ig, celldiv)
%Emmanuel Noutahi- Juillet 2013 Matlab R2013a
%Stochastic Biology Process Simulation with the direct method of Gillepsie
%(Gillespie,J.Phys.Chem.,1977).
%Cell division model
%Promoter A --> Promoter I
%Promoter I --> Promoter A
%Promoter A --> ARNm
%ARNm --> 0
%ARNm --> Protein
%Protein --> 0
%Cell --> Cell/2 (ARNm --> ARNm/2 && Protein --> Protein/2)

clc;
if(isequal(mode, 'Burst'))
    [liste, handle]= divburst(time, trans, darn ,arnnumber,  ag , trad, dprot, protnumber, ig , celldiv);
elseif(isequal(mode,'Continuous'))
    [liste, handle] = divcont(time, trans, darn, arnnumber, trad, dprot, protnumber, celldiv);
    
end
end

%%%%%%%%%%Function divcont

function [liste, handle] = divcont(time, trans, darn, arnnumber, trad, dprot, protnumber, celldiv)

%% Setting Time

t=0;
t_end= time;  % End Time
t_int= 1;  %Time interval to save a point, You can change it
k=1; %wait bar update counter
tot_wait=10^2; %waitbar percent total

if(trad~=0 || dprot ~=0)
    %% Initialization (case with protein)
    %%%  Thattai - van Oudenaarden model (single gene model)
    %%% reaction prom -> ARN + prom (k1)
    %%% reaction ARN -> 0 (k2)
    prom=1; % Promoter state = always on for cont mode
    cell=1;
    %%% array to store results
    i=1; %indice of array
    t_array(i)=t;
    prom_array(i)=prom;
    arn_array(i)=arnnumber;
    prot_array(i)=protnumber; %protein array and initial value
    tic;
    w= waitbar(0, 'updating molecule number...');
    
    while t< t_end
        
        %% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
        % am = hm*cm for the current molecular population numbers.
        a= [trans*prom arnnumber*darn arnnumber*trad protnumber*dprot cell*celldiv];
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
        elseif mu == 3 % translation : protein production
            protnumber=protnumber+1;
        elseif mu == 4 %protein decay
            protnumber= protnumber - 1;
            if(protnumber<=0)
                protnumber=0;
            end
        elseif mu==5 %cell division
            protnumber=round(protnumber*rand);
            arnnumber=round(arnnumber*rand);
        end
        % update  molecule number, species and time for later output
        if t >= i*t_int
            i=i+1;
            t_array(i)=i;
            prom_array(i)=prom;
            arn_array(i)=arnnumber;
            prot_array(i)= protnumber;
        end
        
        % update waitbar
        if t>=k*tot_wait*t_int
            k=k+1;
            waitbar(t/t_end);
        end
    end
    close(w);
    exc_time= toc;
    fprintf('Simulation total time: %f seconds\n', exc_time);
    
    %% simulation statistics
    mean_mRNA_simulation = mean(arn_array);
    stdev_mRNA_simulation = std(arn_array);
    noise_mRNA_simulation = stdev_mRNA_simulation/mean_mRNA_simulation;
    mean_protein_simulation = mean(prot_array);
    stdev_protein_simulation = std(prot_array);
    noise_protein_simulation = stdev_protein_simulation/mean_protein_simulation;
    
    fprintf('\n-----------------Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Prot)= %.2f\nnoise(Prot)= %.2f\n',  mean_mRNA_simulation , noise_mRNA_simulation,mean_protein_simulation, noise_protein_simulation);
    
    %% Plotting result in a fancy way
    hold on;
    
    %%% ARNm output
    handle{1}= subplot(2,2,1);
    plot(t_array,arn_array); %mRNA time series
    xlabel('time (s)');
    ylabel('mRNA no.');
    
    handle{2}= subplot(2,2,2);
    x= min(arn_array):max(arn_array);
    liste.arn= [0:max(arn_array); hist(arn_array, max(arn_array)+1)]';
    bar(x, hist(arn_array, max(arn_array)+1-min(arn_array))./ sum(hist(arn_array)), 'hist');%mRNA histogram
    legend(sprintf(' mean = %g \n stdev = %g',mean_mRNA_simulation,stdev_mRNA_simulation))
    set(gca,'XTick',min(arn_array):round((max(arn_array)-min(arn_array))/2):max(arn_array));
    %%% Fitting a poisson distribution
    prob_distribution = fitdist(arn_array(:),'Poisson');
    hold on;
    probdf = pdf(prob_distribution,x);
    plot(x,probdf,'LineWidth',2, 'color', 'r');
    xlabel('mRNA no.');
    ylabel('percent');
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
    legend(sprintf(' mean = %g \n stdev = %g',mean_protein_simulation,stdev_protein_simulation))
    liste.prot= [0:max(prot_array); hist(prot_array, max(prot_array)+1)]';
    save('data.mat', 't_array', 'prom_array', 'arn_array', 'prot_array');
    
else
    %% Initialization
    %%%  Thattai - van Oudenaarden model (single gene model)
    %%% reaction prom -> ARN + prom (k1)
    %%% reaction ARN -> 0 (k2)
    prom=1; % Promoter state = always on for cont mode
    cell=1;
    %%% array to store results
    i=1; %indice of arrayi=1; %indice of array
    t_array(i)=t;
    prom_array(i)=prom;
    arn_array(i)=arnnumber;
    tic;
    w= waitbar(0, 'updating molecule number...');
    
    while t< t_end
        
        %% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
        % am = hm*cm for the current molecular population numbers.
        a= [trans*prom arnnumber*darn cell*celldiv];
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
            
        elseif mu==3
            arnnumber=round(arnnumber*rand);%*rand);
        end
        % update  molecule number, species and time for later output
        if t >= i*t_int
            i=i+1;
            t_array(i)=i;
            prom_array(i)=prom;
            arn_array(i)=arnnumber;
        end
        
        % update waitbar
        if t>=k*tot_wait*t_int
            k=k+1;
            waitbar(t/t_end);
        end
    end
    close(w);
    exc_time= toc;
    
    fprintf('Simulation total time: %f seconds\n', exc_time);

    %% simulation statistics
    mean_mRNA_simulation = mean(arn_array);
    stdev_mRNA_simulation = std(arn_array);
    noise_mRNA_simulation = stdev_mRNA_simulation/mean_mRNA_simulation;
    fprintf('\n-----------------Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  mean_mRNA_simulation , noise_mRNA_simulation);
    
    %% Plotting result in a fancy way
    cla;
    hold on;
    handle{1}= subplot(2,1,1);
    plot(t_array,arn_array); %mRNA time series
    xlabel('time (s)');
    ylabel('mRNA no.');
    handle{2}= subplot(2,1,2);
    %         length(arn_array), max(arn_array)
    %         for j=0:max(arn_array)
    %             sum(arn_array==j)/ length(arn_array)
    %         end
    %         hist(arn_array, max(arn_array)); %mRNA histogram
    x= min(arn_array):max(arn_array);
    liste.arn= [0:max(arn_array); hist(arn_array, max(arn_array)+1)]';
    bar(x, hist(arn_array, max(arn_array)+1-min(arn_array))./ sum(hist(arn_array)), 'hist');%mRNA histogram
    set(gca,'XTick',min(arn_array):round((max(arn_array)-min(arn_array))/2):max(arn_array));
    legend(sprintf(' mean = %g \n stdev = %g',mean_mRNA_simulation,stdev_mRNA_simulation))
    
    prob_distribution = fitdist(arn_array(:),'Poisson');
    hold on;
    probdf = pdf(prob_distribution,x);
    plot(x,probdf,'LineWidth',2, 'color', 'r');
    xlabel('mRNA no.');
    ylabel('percent');
    
    hold off;
    fprintf('\n-----------------Fitting poisson distribution-----------------\n');
    [H P STATS] = chi2gof(arn_array,'cdf', @(z)poisscdf(z,mean(arn_array)), 'nparams', 1, 'nbins', max(arn_array)+1-min(arn_array))
    save('data.mat', 't_array', 'prom_array', 'arn_array');
end
end


%%% Function burst
function [liste, handle]= divburst(time, trans, darn ,arnnumber,  ag , trad, dprot, protnumber, ig , celldiv)

t=0;
t_end= time;  % End Time
t_int= 1;  %Time interval to save a point, You can change it
k=1; %wait bar update counter
tot_wait=10^2; %waitbar percent total

if(trad==0 && dprot==0)
    
    %% Initialization
    %%%  reaction prom. on -> prom. off
    %%%  reaction prom. off -> prom. on
    %%% reaction prom -> ARN + prom (k1)
    %%% reaction ARN -> 0 (k2)
    prom.number=1; % Promoter number
    prom.state= 1;
    cell=1;
    %%% array to store results
    i=1; %indice of arrayi=1; %indice of array
    t_array(i)=t;
    prom_array(i)=prom.state;
    arn_array(i)=arnnumber;
    tic;
    w= waitbar(0, 'updating molecule number...');
    
    while t< t_end
        
        %% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
        % am = hm*cm for the current molecular population numbers.
        if(prom.state==1)
            a= [prom.number*ig trans*prom.number arnnumber*darn cell*celldiv];
        elseif(prom.state==0)
            a= [prom.number*ag arnnumber*darn cell*celldiv];
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
                elseif mu == 3 %celldiv
                    arnnumber=round(arnnumber*rand);
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
                elseif mu == 4 %celldiv
                    arnnumber=round(arnnumber*rand);
                end
        end
        
        % update  molecule number, species and time for later output
        if t >= i*t_int
            i=i+1;
            t_array(i)=i;
            prom_array(i)=prom.state;
            arn_array(i)=arnnumber;
        end
        
        % update waitbar
        if t>=k*tot_wait*t_int
            k=k+1;
            waitbar(t/t_end);
        end
    end
    close(w);
    exc_time= toc;
    
  
    %% simulation statistics
    mean_mRNA_simulation = mean(arn_array);
    stdev_mRNA_simulation = std(arn_array);
    noise_mRNA_simulation = stdev_mRNA_simulation/mean_mRNA_simulation;
    fprintf('\n-----------------Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\n',  mean_mRNA_simulation , noise_mRNA_simulation);
    
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
    legend(sprintf(' mean = %g \n stdev = %g',mean_mRNA_simulation,stdev_mRNA_simulation))
    
    xlabel('mRNA no.');
    ylabel('count');
    
    handle{3}= subplot(2,2,[3 4]);
    t=1;
    for i=1:numel(prom_array)-1
        timearray(t)=t_array(i);
        timearray(t+1)=t_array(i);
        promoarray(t)=prom_array(i);
        if(prom_array(i)==prom_array(i+1))
            promoarray(t+1)=prom_array(i);
        else
            promoarray(t+1)=prom_array(i+1);
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
    
%     fprintf('\n-----------------Fitting poisson distribution-----------------\n');
%     [H P STATS] = chi2gof(arn_array,'cdf', @(z)poisscdf(z,mean(arn_array)), 'nparams', 1, 'nbins', max(arn_array)+1-min(arn_array))
    save('data.mat', 't_array', 'prom_array', 'arn_array');
    
else %%%---> with protein model
    %% Initialization
    %%% reaction prom. on -> prom. off
    %%% reaction prom. off -> prom. on
    %%% reaction prom -> ARN + prom (k1)
    %%% reaction ARN -> 0 (k2)
    %%% reaction ARN-> protein
    %%% reaction protein->0;
    
    prom.number=1; % Promoter number
    prom.state= 1;
    cell=1;
    %%% array to store results
    i=1; %indice of arrayi=1; %indice of array
    t_array(i)=t;
    prom_array(i)=prom.state;
    arn_array(i)=arnnumber;
    prot_array(i)=protnumber; %protein array and initial value
    tic;
    w= waitbar(0, 'updating molecule number...');
    
    while t< t_end
        
        %% Step 1 : Calculate and store the M quantities a1 = h1*c1 ...
        % am = hm*cm for the current molecular population numbers.
        if(prom.state==1)
            a= [prom.number*ig trans*prom.number arnnumber*darn arnnumber*trad protnumber*dprot cell*celldiv];
        elseif(prom.state==0)
            a= [prom.number*ag arnnumber*darn arnnumber*trad protnumber*dprot cell*celldiv];
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
                elseif mu==3 %traduction
                    protnumber=protnumber+1;
                elseif mu == 4 % prot decay
                    protnumber = protnumber - 1;
                    if(protnumber<=0)
                        protnumber=0;
                    end
                
                elseif mu == 5 % celldiv
                    arnnumber=round(arnnumber*rand);
                    protnumber=round(protnumber*rand);
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
               elseif mu == 6 % celldiv
                    arnnumber=round(arnnumber*rand);
                    protnumber=round(protnumber*rand);
                end
        end
        
        % update  molecule number, species and time for later output
        if t >= i*t_int
            i=i+1;
            t_array(i)=i;
            prom_array(i)=prom.state;
            arn_array(i)=arnnumber;
            prot_array(i)=protnumber;
        end
        
        % update waitbar
        if t>=k*tot_wait*t_int
            k=k+1;
            waitbar(t/t_end);
        end
    end
    close(w);
    exc_time= toc;
    
    fprintf('Simulation total time: %f seconds\n', exc_time);

    %% simulation statistics
    mean_mRNA_simulation = mean(arn_array);
    stdev_mRNA_simulation = std(arn_array);
    noise_mRNA_simulation = stdev_mRNA_simulation/mean_mRNA_simulation;
    mean_protein_simulation = mean(prot_array);
    stdev_protein_simulation = std(prot_array);
    noise_protein_simulation = stdev_protein_simulation/mean_protein_simulation;
    
    fprintf('\n-----------------Simulation statistics-----------------\n\nmean(mRNA)= %.2f\nnoise(mRNA)= %.2f\nmean(Prot)= %.2f\nnoise(Prot)= %.2f\n',  mean_mRNA_simulation , noise_mRNA_simulation,mean_protein_simulation, noise_protein_simulation);
    
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
    legend(sprintf(' mean = %g \n stdev = %g',mean_mRNA_simulation,stdev_mRNA_simulation))
    
    xlabel('mRNA no.');
    ylabel('count');
    
    handle{3}= subplot(3,2,[5 6]);
    t=1;
    for i=1:numel(prom_array)-1
        timearray(t)=t_array(i);
        timearray(t+1)=t_array(i);
        promoarray(t)=prom_array(i);
        if(prom_array(i)==prom_array(i+1))
            promoarray(t+1)=prom_array(i);
        else
            promoarray(t+1)=prom_array(i+1);
        end
        t=t+2;
    end
    plot(timearray, promoarray, '-r*');
    
    ylabel('Promoter state.');
    xlabel('time (s)');
    set(gca, 'YLim', [-1 2]);
    set(gca, 'YTick', [0 1]);
    set(gca,'YTickLabel',{'off', 'on'});
%     
%     fprintf('\n-----------------Fitting poisson distribution-----------------\n');
%     [H P STATS] = chi2gof(arn_array,'cdf', @(z)poisscdf(z,mean(arn_array)), 'nparams', 1, 'nbins', max(arn_array)+1-min(arn_array))
    save('data.mat', 't_array', 'prom_array', 'arn_array');
    
    %%% Protein output
    handle{3}= subplot(3,2,3);
    plot(t_array,prot_array);
    xlabel('time (s)');
    ylabel('protein no.');
    handle{4}=subplot(3,2,4);
    hist(prot_array, max(prot_array)-min(prot_array)+1); %protein histogram
    xlabel('protein no.');
    ylabel('occurrence');
    legend(sprintf(' mean = %g \n stdev = %g',mean_protein_simulation,stdev_protein_simulation))
    liste.prot= [0:max(prot_array); hist(prot_array, max(prot_array)+1)]';
    save('data.mat', 't_array', 'prom_array', 'arn_array', 'prot_array');
    
end
end