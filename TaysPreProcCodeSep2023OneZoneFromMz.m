 %% From Onezone_frommz_preprocessingcode_july2023 - edited only keeping R2
 
    numZones=6; %num zones recorded
    
    I = transpose(data(1:numZones+1:end) * ampsPerVolt);
    
    % 2. Retrieve Raw Voltage Signals -----------------------------------------
    
    % intialize raw voltage matrix
    V = zeros(length(data)/(numZones+1), numZones);
     
    
    % populate matrix with raw voltage signal for pore and each zone
    for i = 1:numZones
        V(:,i) = transpose(data(i+1:numZones+1:end));
    end
     

    % create total voltage vector
    Vtotal = sum(V(:,1:numZones),2);
    
    % 3. Compute Time Vector At Original Sample Rate --------------------------

    

    t = (0:size(I,1)-1).' ./ sampleRate; 
     
    %4. Calculate Resistance Signals -----------------------------------------
    
    f=figure ;
    set(f,'Position', [2000, 500, 2000, 400])
    hold on %plotting raw voltage of all zones 
    for i = 1:numZones 
        plot(V(:,i))

    end
    
    R = V ./ I; %ohms law
    
    R=-R(:,1); %only using first zone

    %if first few samples reach 0:
    %R=R(300161:end);
    
    %% Optional tiled plotting
    f=figure;
    set(f,'Position', [2000, 500, 2000, 400])
    plot(R)
    
    %% 5. Filter Resistance Signals --------------------------------------------
    fpass = 50; % passband frequency of the lowpass filter (in Hz); remove high frequency noise and keep signal under 50 Hz
    window = 71; % averaging window width (in units of samples) originally 31
    window2= 71;
    polynomial=3;
    polynomial2=1;
    polynomial3=2;
    
    numsgo=10;
%     Rfalt1= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial,window);  
     Rfalt2= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial,window);  
%      Rfalt3= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial2,window2);  
%     Rfalt4= sgolayfilt(lowpass(R, fpass, sampleRate),polynomial3,101);  
    for i=numsgo

%         Rfalt1= sgolayfilt(Rfalt1,polynomial,window); 
         Rfalt2= sgolayfilt(Rfalt2,polynomial,window); 
%         Rfalt3= sgolayfilt(Rfalt3,polynomial2,window2); 
%         Rfalt4= sgolayfilt(Rfalt4,polynomial3,101);
    end
%     Rlow=lowpass(R,50,sampleRate);
%     Rlow=Rlow(window2:end-window2,:);
    
%     Rfalt1=Rfalt1(window2:end-window2,:);
     Rfalt2=Rfalt2(window:end-window,:);
%     Rfalt3=Rfalt3(window2:end-window2,:);
%     Rfalt4=Rfalt4(window2:end-window2,:);
    % plot(Rlow(:,1),'LineWidth',2.0)
    % plot(Rfalt1(:,1),'LineWidth',2.0)
    %plot(Rfalt2(:,1),'LineWidth',2.0)
    
%     %% for high noise - try running this 
%     Rfalt2=Rfalt4;
    %% plotting
    
    
   
    f1=figure;
   set(f1,'Position', [2000, 500, 2000, 400]) 

     %plot(Rlow(:,1),'LineWidth',2.0)
    hold on 
    %plot(Rfalt1(:,1),'LineWidth',2.0)
%     plot(Rfalt2(:,1),'LineWidth',1.0)
%     plot(Rfalt3(:,1),'LineWidth',1.0)
%     plot(Rfalt4(:,1),'LineWidth',2.0)
    
    
     Rf=Rfalt2(1:end-200,1); %getting rid of start and end drops
     %Rf = Rf(1:462856);
     hold on
     plot(Rf)
    
    
    %% get stderror
    
     %pick region for stdev calculation --> un
     % 
     % comment when you are doing the
      
     
     
     stderror=zeros(2,1);
      
    
            d=datacursormode(f1);
              d.Enable='on';
             d.DisplayStyle='window';
             startstdi=input('click where you want region to start, then press enter');
             if isempty(startstdi)==1 
                 vals=getCursorInfo(d);
                 startstd=vals.Position(1,1);
            end
             d.Enable='off';
             d.Enable='on';
             endstdi=input('click where you want region to end, then press enter');
             if isempty(endstdi)==1
                 vals=getCursorInfo(d);
                 endstd=vals.Position(1,1);
             end
             
                stderror(1,1)=std(Rf(startstd:endstd,1));
            
    
    %% fit baseline
    %for asls if not downsampled 
    lambda1 = 1e12; %larger means smoother background 
    % lambda2=1e5; 
    %for asls if downsampled
    % lambda2=5e3; %50
    % lambda3=1e5; %1000
    p=0; % less than 0.5 means negative is strongly supressed %0.02
    maxiter=20; % maximum iterations 
    noisez1=stderror(1,1);
    
    
    aslsparam=struct('lambda', lambda1, 'p', p, 'max_iter', maxiter, 'noise_margin', noisez1);
    
    
    
    
    yasls=zeros(length(Rf(:,1)),1);
%     yasls2=yasls;
%     yasls3=yasls;
%     ydetrend=yasls;
     yas2det=yasls;
%     yas3det=yasls;
    
           yasls=(ASLS2(Rf,aslsparam));
           %yasls2(:,i)=ASLS2(Rf(:,i),asls2param(i,1));
           %yasls3(:,i)=ASLS2(Rf(:,i),asls3param(i,1));
    
    %% plot yasls results 
    f=figure;
    set(f,'Position', [2000, 500, 2000, 400])
    plot(yasls)
    
    
    
    
        hold on 
       
        plot(Rf)
        %plot(yasls2(:,i))
        %plot(yasls3(:,i))
%          plot(Rfalt2)
    
    %Rfsmooth=smoothdata(Rf,'movmean',91);
    % plot(Rfsmooth(:,1))
     %yas3det will follow the sizing pulse, but MAY NOT go into the sq pulse depending on transit time, rely on yas2det for that  
    
    
        
     %% subtract baseline_replot
    f1=figure;
    set(f1,'Position', [2000, 500, 2000, 400]) 
    hold on
    
        ydetrend=Rf-yasls;
         yas2det=Rf-yasls;
%         yas3det=Rfalt3-yasls;
%         plot(ydetrend)
%         plot(yas2det)
        plot(yas2det,'LineWidth',1)
      
    %stop here and save your pre-processed data 

    %t_filtered = (0:length(yas2det)-1).' / sampleRate ; %time vector for
    %plotting against time
       
    %% calculate noise 
    %pick region for stdev calculation --> uncomment when

    
           d=datacursormode(f1);
              d.Enable='on';
             d.DisplayStyle='window';
             startstdi=input('click where you want region to start, then press enter');
             if isempty(startstdi)==1
                 vals=getCursorInfo(d);
                 startstd=vals.Position(1,1);
             end
             d.Enable='off';
            
    
             d.Enable='on';
             endstdi=input('click where you want region to end, then press enter');
             if isempty(endstdi)==1
                 vals=getCursorInfo(d);
                 endstd=vals.Position(1,1);
             end
              std4errorplot=zeros(endstd-startstd+1,3);
           
                stderror(2,1)=std(yas2det(startstd:endstd,1));
                std4errorplot(:,1)=yas2det(startstd:endstd,1);
        
    
    
    
    %% reset stderror
    
    stderror(1,1)=stderror(2,1);
    
    
       
    
    % %% pick pulse to analyze 
    % wp=16000;
    %  d=datacursormode(f1);
    %  d.Enable='on';
    %  d.DisplayStyle='window';
    %  startxchosen=input('click where you want region to start, then press enter');
    % if isempty(startxchosen)==1
    %     vals=getCursorInfo(d);
    %     startv=vals.Position(1,1);
    % end
    % d.Enable='off';
    % 
    % 
    % 
    % f2=figure;
    % sz=10;
    % tl=tiledlayout('flow');
    % 
    % ax1=nexttile;
    % hold(ax1,'on')
    % scatter(t(startv-100:startv+wp,1),ydetrend(startv-100:startv+wp,1),sz,'m')
    % plot(t(startv-100:startv+wp,1),ydetrend(startv-100:startv+wp,1),'m')
    % scatter(t(startv-100:startv+wp,1),ydetrend(startv-100:startv+wp,2),sz,'k')
    % plot(t(startv-100:startv+wp,1),ydetrend(startv-100:startv+wp,2),'k')
    % scatter(t(startv-100:startv+wp,1),ydetrend(startv-100:startv+wp,3),sz,'b')
    % plot(t(startv-100:startv+wp,1),ydetrend(startv-100:startv+wp,3),'b')
    % 
    % pulses=zeros(1,8);
    % s=1;
    % e=0;
    % for i=1:8
    %  d=datacursormode(f2);
    %  d.Enable='on';
    %  d.DisplayStyle='window';
    %  startxchosen=input('click to record a point and then press enter');
    % if isempty(startxchosen)==1
    %     vals=getCursorInfo(d);
    %     pulses(1,i)=vals.Position(1,1);
    % end
    % d.Enable='off';
    %        
    %     
    % end
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    % 
    % 
