function [R_inds,Q_inds,S_inds,QRS_On,QRS_Off,ecg_hat,Peak_activities]=F2AMM3(ecg,FS,DisplayProgress)

ecg=ecg-mean(ecg);


%             [N,Wn] = buttord(0.15,0.20,0.5,20) ;
%             [b,a] = butter(N,Wn) ;
%             ecg = filtfilt(b,a,ecg);
%             
%         %%Filtering Section design and implementation
%         
%         Cuttoff=40;
%         Wp=(Cuttoff/FS)*2;    % Cuttoff frequency
%         %FIR
%         b=fir1(100,Wp);
%         a=1;
% %         %IIR
% %         Transition=Cuttoff+8;
% %         Ws=(Transition/FS)*2; % Transtions phase
% %         Rp=0.5; %Ripple
% %         RS=20; % attenuation in the stopband
% %         [N,Wn] = buttord(Wp,Ws,Rp,RS);
% %         [b,a] = butter(N,Wn) ;
%         %%Applying the filter
%         %freqz(b,a,2000,FS);
%         %pause();
%         ecg = filtfilt(b,a,ecg);

            
% ecg=ecg/max(ecg);
ecglength=length(ecg);
sigmean=sign(mean(ecg))*mean(ecg);
ecg=ecg-sigmean;



%% for onset and offset optimization
% Opt_lim=FS;
tolorance=30; %in milliseconds
Opt_lim=max(1,round(tolorance/(1000/FS)));
%% Creating the Artificial QRS complex structure


G=[];
G=[G 0:-0.01:-.2;];
G=[G -0.2:0.01:1.2;];
G=[G 1.2:-0.01:-.3;];
G=[G -0.3:0.01:0;];

%Estimating the number of samples in a normal QRS complex based on Sampling
% Normal QRS length 0.06-0.12 sec (0.1s considered here)
Samples=floor(0.08*FS); if(mod(Samples,2)==0), Samples=Samples+1; end

% Creating a high resolution QRS complex
% G_Vals=[0 -.2 1.2 -.3 0] ;
G_Vals=[0 -.1 8 -.1 0] ;
G_Inds=[0 1 2 3 4]*(Samples/4); % times Samples/4 to make a 10 milliseconds QRS complex
g = interp1(G_Inds,G_Vals,0:G_Inds(end)/(Samples-1):G_Inds(end));
% g=[zeros(1,10) g zeros(1,10)];
struc_len=Samples;%+20; %Length of the Created Structure
gs=g; % gs will keep track of the changes in the adaptiv g structure.

% gs=g; % gs will keep track of the changes in the adaptiv g structure.
% finding fiducial points of the created g structure
% [g_R_val,g_R_ind]=max(g);
% [g_S_val,g_S_ind]=min(g(g_R_ind:end));
% [g_Q_val,g_Q_ind]=min(g(1:g_R_ind));
% g_Offset_ind=find(g(g_R_ind+g_S_ind-1:end)==0);
% G_Vals=[0 g_Q_val g_R_val g_S_val 0] ;
% G_Inds=[g_Q_ind g_R_ind g_R_ind+g_S_ind-1 g_R_ind+g_S_ind-1+g_Offset_ind(1)-1];

%% Adaptive QRS complex detection Initialization
% The running window through ECG
t=0.5; %window length in seconds
window=round(t*FS);

% Fiducial points of the QRS Complex
Peak_Ind=1;
QRS_On=zeros(length(ecg),1); QRS_Off=zeros(length(ecg),1); Q_inds=zeros(length(ecg),1); R_inds=zeros(length(ecg),1); S_inds=zeros(length(ecg),1); 

RR_Intervals=zeros(length(ecg),1);
RR_Ind=1;

% Detected Complexes values and indecies
QRS=zeros(length(ecg),5); QRS_Inds=zeros(length(ecg),4); 

% Adaptive threshold parameters for detecting and removing false QRS complexes
Peak_activities=zeros(length(ecg),1);
ecg_hat=zeros(length(ecg),1);
ecg_hat1=zeros(length(ecg),1);
ecg_diff=zeros(length(ecg),1);
% ecg_Activity=zeros(length(ecg),1);
% Inactivity=85;%round(.15*FS); %floor(Samples*.3);
Inactivity=10;%round(.1*FS); %floor(Samples*.3);

activity_th=.8;
QRS_len=0.5*Samples;
QRS_NO=1; % number of QRS complexes participating in creating the new QRS structure for MM operations.
% Learning rate for the QRS Structure in MM operations.
Etha=0.9;%Etha=1 then QRS_NO= number of QRS complexes determining the new QRS shape


global g_amplitude;
g_amplitude=0.7;
global g_length;

%% QRS detector parameters
R_Reversed=0; %for negative R peaks
not_changed=0; Onset=0; checked=1; wrong_R=0;
Q_pos_changed=0; S_pos_changed=0; Qon_second_chance=1; Qoff_second_chance=1;
Qon_temp=0; Qoff_temp=0;

%% How much of the singal is parsed
ECG_Parsed_Ind=struc_len+1;
%% Detecting QRS fiducial points and Complexes 
N=length(ecg);

%% Displaying the progress in Terminal
PERCENT=1;

for i=1:window:N-window-2*struc_len-1
    if(DisplayProgress==1 && mod(PERCENT,300)==0)
        disp([num2str((i/(N-window-2*struc_len))*100), '% of the tape DONE.']);
        PERCENT=1;
    else
        PERCENT=PERCENT+1;
    end
    %% Performing Morphological Opertions
    Window=ecg(i:i+window+2*struc_len);
%     if(i==1), g=g/max(g); g=g*max(Window)*g_amplitude; end
try
%     if(i==1),
        Window_now=ecg(i:i+2*FS);
         detrended_window=detrend(Window_now);
        %detrended_window=Window_now;
        g=g/max(g); g=g*max(detrended_window)*g_amplitude; G_Vals=G_Vals*max(detrended_window)*g_amplitude; gs=g;
%     end
catch
    
end

%     Window=Window/max(Window);
%     Window=dilation(Window,g)+erosion(Window,g);
%     Window=Window(:);
    opened=open_ecg(Window,g);
    closed=close_ecg(Window,g);
    That=Window-opened';
    Bhat=Window-closed';
    Hat=(Bhat+That)/2;

    % Suppressing small activities in HAT and computing diff and Activity
    Hat(abs(Hat)<=.1*std(Hat))=0;
    Hat1=Hat;
    diff1=Hat(2:end)-Hat(1:end-1);
    diff1(abs(diff1)<=0.995*std(diff1))=0;
    rules=find(diff1<0.01*std(Hat));
    Hat1(rules+1)=0;    
%     Activity=abs(diff);

    Hat_Len=length(Hat(struc_len+2:end-struc_len));
    ecg_hat(ECG_Parsed_Ind:ECG_Parsed_Ind+Hat_Len-1)=Hat(struc_len+2:end-struc_len);
    ecg_hat1(ECG_Parsed_Ind:ECG_Parsed_Ind+Hat_Len-1)=Hat1(struc_len+2:end-struc_len);
    ecg_diff(ECG_Parsed_Ind:ECG_Parsed_Ind+Hat_Len-1)=diff1(struc_len+1:end-struc_len);
%     ecg_Activity(ECG_Parsed_Ind:ECG_Parsed_Ind+Hat_Len-1)=Activity(struc_len+1:end-struc_len);
    ECG_Parsed_Ind=ECG_Parsed_Ind+Hat_Len;
        
      
    %% Detecting Fiducial points
    j=struc_len+2;
    while j<=length(Hat)-struc_len
        if (Hat(j)==0), not_changed=not_changed+1;
        else
            if(Onset==0), not_changed=0; Onset_ind=j; QRS_On(Peak_Ind)= Onset_ind+i-1; Onset=1; checked=0; end % -1 because of starting index=1 not 0
             not_changed=0; 
            Offset_ind=j; 
            Offset_Window_ind=i-1;
        end
%             if(checked==0 && (Offset_ind + Offset_Window_ind - QRS_On(Peak_Ind))>60)
%             Inactivity=Inactivity-5;
%             end
        if(not_changed > Inactivity && checked==0)
            checked=1;
           %% Inactivity=EthaD;
            not_changed=0;
            Onset=0;
            QRS_Off(Peak_Ind)= Offset_ind + Offset_Window_ind;    
            Peak_activities(Peak_Ind)= sum(abs(ecg_hat(QRS_On(Peak_Ind):QRS_Off(Peak_Ind))));
            if(sum(abs(ecg_hat(QRS_On(Peak_Ind):QRS_Off(Peak_Ind))))>0.01)
              if(~isempty(find(ecg_hat(QRS_On(Peak_Ind):QRS_Off(Peak_Ind))>0, 1)) || ~isempty(find(ecg_hat1(QRS_On(Peak_Ind):QRS_Off(Peak_Ind))~=0, 1)))
                  if(~isempty(find(ecg_hat(QRS_On(Peak_Ind):QRS_Off(Peak_Ind))<0, 1)) || ~isempty(find(ecg_hat1(QRS_On(Peak_Ind):QRS_Off(Peak_Ind))~=0, 1)))                  
                        QRS_temp=zeros(1,5);
                        QRS_Ind_temp=zeros(1,4);

                        %% R peak Extraction                             
                        [R_max_val,R_max_ind]=max(ecg_hat(QRS_On(Peak_Ind):QRS_Off(Peak_Ind)));
                        [R_min_val,R_min_ind]=min(ecg_hat(QRS_On(Peak_Ind):QRS_Off(Peak_Ind)));
                                                                      
%                         [S_max_val,S_max_ind]=min(ecg(QRS_On(Peak_Ind):QRS_On(Peak_Ind)+R_max_ind-1));
%                         [Q_max_val,Q_max_ind]=min(ecg(QRS_On(Peak_Ind)+R_max_ind-1:QRS_Off(Peak_Ind)));
%                         [S_min_val,S_min_ind]=max(ecg(QRS_On(Peak_Ind):QRS_On(Peak_Ind)+R_min_ind-1));
%                         [Q_min_val,Q_min_ind]=max(ecg(QRS_On(Peak_Ind)+R_min_ind-1:QRS_Off(Peak_Ind)));
%                         temp1=R_max_val-S_max_val-Q_max_val+abs(S_max_val)+abs(Q_max_val);
%                         temp2=abs(R_min_val-S_min_val-Q_min_val)+abs(S_min_val)+abs(Q_min_val);
                        %if(abs(ecg(QRS_On(Peak_Ind)+R_max_ind-1)-ecg(QRS_On(Peak_Ind))) >= .5*abs(ecg(QRS_On(Peak_Ind))-ecg(QRS_On(Peak_Ind)+R_min_ind-1)))
                        if(abs(ecg(QRS_On(Peak_Ind)+R_max_ind-1)-ecg(QRS_On(Peak_Ind))) >= .4*abs(ecg(QRS_On(Peak_Ind))-ecg(QRS_On(Peak_Ind)+R_min_ind-1)))
%                         if(R_max_val>0.15)
                            R_Reversed=0;
%                             [R_max_val,R_max_ind]=max(ecg(QRS_On(Peak_Ind):QRS_Off(Peak_Ind)));
                            R_ind= QRS_On(Peak_Ind) + R_max_ind(1)-1 ;
                        else
                            R_Reversed=1;
%                             [R_min_val,R_min_ind]=min(ecg(QRS_On(Peak_Ind):QRS_Off(Peak_Ind)));
                            R_ind= QRS_On(Peak_Ind) + R_min_ind(1)-1 ;
                        end
                        if(R_Reversed==1)
                            if(ecg(R_ind+1)<=ecg(R_ind))
                                while ecg(R_ind+1)<=ecg(R_ind)
                                   R_ind=R_ind+1;
                                end
                                R_ind=min(R_ind,length(ecg_hat));
                                QRS_Off(Peak_Ind)= max(R_ind,QRS_Off(Peak_Ind));
                            elseif(ecg(R_ind-1)<=ecg(R_ind))
                                while ecg(R_ind-1)<=ecg(R_ind)
                                   R_ind=R_ind-1;
                                end
                                QRS_On(Peak_Ind)= R_ind;
                            end
                        else
                            if(ecg(R_ind+1)>=ecg(R_ind))
                                while ecg(R_ind+1)>=ecg(R_ind)
                                   R_ind=R_ind+1;
                                end
                                R_ind=min(R_ind,length(ecg_hat));
                                QRS_Off(Peak_Ind)= max(R_ind,QRS_Off(Peak_Ind));
                            elseif(ecg(R_ind-1)>=ecg(R_ind))
                                while ecg(R_ind-1)>=ecg(R_ind)
                                   R_ind=R_ind-1;
                                end
                                QRS_On(Peak_Ind)= R_ind;
                            end                    
                        end

                        %Checking for Wrong Rs
                        if (Peak_Ind > 1 && ((R_ind - R_inds(Peak_Ind-1))/FS)<0.25), %Humanly impossible
                            if(Peak_activities(Peak_Ind-1)>Peak_activities(Peak_Ind))                                
                                wrong_R=1;
                            else
                                QRS_On(Peak_Ind-1)= QRS_On(Peak_Ind);
                                QRS_Off(Peak_Ind-1)=QRS_Off(Peak_Ind);
                                Peak_activities(Peak_Ind-1)=Peak_activities(Peak_Ind);
                                R_inds(Peak_Ind-1)=R_ind;
                                QRS_temp(3)=ecg(R_ind);
                                QRS_Ind_temp(2)=R_ind;
                                Peak_Ind=Peak_Ind-1;
                            end
                        else         
                            R_inds(Peak_Ind)=R_ind;
                            QRS_temp(3)=ecg(R_ind);
                            QRS_Ind_temp(2)=R_ind;
                        end
%%heuristically removing FPs!
%                         if (Peak_Ind > 1)
%                             No_Beats=min(Peak_Ind-1,6);
%                             Beat_Interval=(R_inds(Peak_Ind)-R_inds(Peak_Ind-1))/FS;
% %                             Beat_Interval=(QRS_On(Peak_Ind)-QRS_Off(Peak_Ind-1))/FS;
%                             if(Beat_Interval<2 && Peak_activities(Peak_Ind)*Beat_Interval<.2*mean(Peak_activities(Peak_Ind-No_Beats:Peak_Ind-1)))
%                                 wrong_R=1;
%                             end                            
%                         end
                        if(wrong_R==0)
                            %% Q peak detection
                            if(R_Reversed==0)
                                [val,ind]=min(ecg_hat(QRS_On(Peak_Ind):R_ind)); %-1
                                Q_inds(Peak_Ind)= QRS_On(Peak_Ind)+ind(1)-1; % -1 because of starting index=1 not 0
                                        %% Fine Tune Q and onset
                                        while (Q_inds(Peak_Ind)-1>0) && ecg(Q_inds(Peak_Ind)-1)<=ecg(Q_inds(Peak_Ind))
                                            Q_inds(Peak_Ind)=Q_inds(Peak_Ind)-1;
                                            Q_pos_changed=1;
                                        end
%                                         Q_pos_changed=1;
                                        lim=Opt_lim;
                                        if(Q_pos_changed==1)
                                            Q_pos_changed=0;
                                            QRS_On(Peak_Ind)=Q_inds(Peak_Ind);
                                            while lim>0 && (QRS_On(Peak_Ind)-1>0) && (ecg(QRS_On(Peak_Ind)-1)>=ecg(QRS_On(Peak_Ind)) || Qon_second_chance==1)
                                                if(ecg(QRS_On(Peak_Ind)-1)>=ecg(QRS_On(Peak_Ind)))
                                                    QRS_On(Peak_Ind)=QRS_On(Peak_Ind)-1;
                                                    lim=lim-1;
                                                    %Qon_second_chance=1;
                                                else
                                                    Qon_second_chance=0;
                                                    Qon_temp=QRS_On(Peak_Ind);
                                                    QRS_On(Peak_Ind)=QRS_On(Peak_Ind)-1;
                                                end
                                            end
                                            if(Qon_temp~=0)
                                                if(ecg(QRS_On(Peak_Ind))-ecg(Qon_temp)<0.06)
                                                    QRS_On(Peak_Ind)=Qon_temp;
                                                end
                                                Qon_temp=0;
                                                Qon_second_chance=1;
                                            end
                                        end
                            else
                                [val,ind]=max(ecg_hat(QRS_On(Peak_Ind):R_ind)); %-1
                                Q_inds(Peak_Ind)=QRS_On(Peak_Ind)+ind(1)-1; % -1 because of starting index=1 not 0
                                        %% Fine Tune Q and onset
                                        while (Q_inds(Peak_Ind)-1>0) && ecg(Q_inds(Peak_Ind)+1)>=ecg(Q_inds(Peak_Ind))
                                            Q_inds(Peak_Ind)=Q_inds(Peak_Ind)+1;
                                            Q_pos_changed=1;
                                        end
%                                         Q_pos_changed=1;
                                        lim=Opt_lim;
                                        if(Q_pos_changed==1)
                                            Q_pos_changed=0;
                                            onset=QRS_On(Peak_Ind);
                                            QRS_On(Peak_Ind)=Q_inds(Peak_Ind);
                                            while lim>0 && (QRS_On(Peak_Ind)-1>onset) && (QRS_On(Peak_Ind)-1>0) && (ecg(QRS_On(Peak_Ind)-1)<=ecg(QRS_On(Peak_Ind)) || Qon_second_chance==1)
                                                if(ecg(QRS_On(Peak_Ind)-1)<=ecg(QRS_On(Peak_Ind)))
                                                    QRS_On(Peak_Ind)=QRS_On(Peak_Ind)-1;
                                                    lim=lim-1;
                                                    %Qon_second_chance=1;
                                                else
                                                    Qon_second_chance=0;
                                                    Qon_temp=QRS_On(Peak_Ind);
                                                    QRS_On(Peak_Ind)=QRS_On(Peak_Ind)-1;
                                                end
                                            end
                                            if(Qon_temp~=0)
                                                if(ecg(QRS_On(Peak_Ind))-ecg(Qon_temp)< -0.06)
                                                    QRS_On(Peak_Ind)=Qon_temp;
                                                end
                                                Qon_temp=0;
                                                Qon_second_chance=1;
                                            end
                                        end
                            end
                            QRS_temp(2)= ecg(Q_inds(Peak_Ind));
                            QRS_Ind_temp(1)= Q_inds(Peak_Ind);
                            %% S peak detection
                            if(R_Reversed==0)
                                [val,ind]=min(ecg_hat(R_ind:QRS_Off(Peak_Ind)));% -1
                                S_inds(Peak_Ind)=R_ind+ind(1)-1;
                                        %% Fine Tune S and offset
                                        while (S_inds(Peak_Ind)+1 <= ecglength) && ecg(S_inds(Peak_Ind)+1)<=ecg(S_inds(Peak_Ind))
                                            S_inds(Peak_Ind)=S_inds(Peak_Ind)+1;
                                             S_pos_changed=1;
                                        end
%                                         S_pos_changed=1;
                                         lim=Opt_lim;
                                         if(S_pos_changed==1)
                                            S_pos_changed=0;
                                            offset=QRS_Off(Peak_Ind);
                                            QRS_Off(Peak_Ind)=S_inds(Peak_Ind);
                                            while lim>0 && (QRS_Off(Peak_Ind)+1 <= offset) && ((QRS_Off(Peak_Ind)+1 <= ecglength) && ecg(QRS_Off(Peak_Ind)+1)>=ecg(QRS_Off(Peak_Ind)) || Qoff_second_chance==1)
                                                if(ecg(QRS_Off(Peak_Ind)+1)>=ecg(QRS_Off(Peak_Ind)))
                                                    QRS_Off(Peak_Ind)=QRS_Off(Peak_Ind)+1;
                                                    lim=lim-1;
                                                    %Qoff_second_chance=1;
                                                else
                                                    Qoff_second_chance=0;
                                                    Qoff_temp=QRS_Off(Peak_Ind);
                                                    QRS_Off(Peak_Ind)=QRS_Off(Peak_Ind)+1;
                                                end
                                            end
                                            if(Qoff_temp~=0)
                                                if(ecg(QRS_Off(Peak_Ind))-ecg(Qoff_temp)<0.06)
                                                    QRS_Off(Peak_Ind)=Qoff_temp;
                                                end
                                                Qoff_temp=0;
                                                Qoff_second_chance=1;
                                            end
                                         end
                            else
                                [val,ind]=max(ecg_hat(R_ind:QRS_Off(Peak_Ind)));% -1
                                S_inds(Peak_Ind)=R_ind+ind(1)-1;
                                        %% Fine Tune S and offset
                                        while (S_inds(Peak_Ind)+1 <= ecglength) &&  ecg(S_inds(Peak_Ind)+1)>=ecg(S_inds(Peak_Ind))
                                            S_inds(Peak_Ind)=S_inds(Peak_Ind)+1;
                                             S_pos_changed=1;
                                        end
%                                         S_pos_changed=1;
                                         lim=Opt_lim;
                                         if(S_pos_changed==1)
                                            S_pos_changed=0;
                                            QRS_Off(Peak_Ind)=S_inds(Peak_Ind);
                                            while lim>0 && (ecg(QRS_Off(Peak_Ind)+1)<=ecg(QRS_Off(Peak_Ind)) || Qoff_second_chance==1)
                                                if(ecg(QRS_Off(Peak_Ind)+1)<=ecg(QRS_Off(Peak_Ind)))
                                                    QRS_Off(Peak_Ind)=QRS_Off(Peak_Ind)+1;
                                                    lim=lim-1;
                                                    %Qoff_second_chance=1;
                                                else
                                                    Qoff_second_chance=0;
                                                    Qoff_temp=QRS_Off(Peak_Ind);
                                                    QRS_Off(Peak_Ind)=QRS_Off(Peak_Ind)+1;
                                                end
                                            end
                                            if(Qoff_temp~=0)
                                                if(ecg(QRS_Off(Peak_Ind))-ecg(Qoff_temp)< -0.06)
                                                    QRS_Off(Peak_Ind)=Qoff_temp;
                                                end
                                                Qoff_temp=0;
                                                Qoff_second_chance=1;
                                            end
                                         end
                            end
                                    QRS_temp(4)=ecg(S_inds(Peak_Ind));
                                    QRS_Ind_temp(3)=S_inds(Peak_Ind);
                                % QRS Onset and Offset    
                                QRS_temp(1)=ecg(QRS_On(Peak_Ind));
                                QRS_temp(5)= ecg(QRS_Off(Peak_Ind));
                                QRS_Ind_temp(4)= QRS_Off(Peak_Ind);                                               
                            %% Checking if the detected Complex is false positive
                            update=1;
                            if((QRS_Off(Peak_Ind)-QRS_On(Peak_Ind))<=.2*Samples)
                                update=0;
                            end
                            if((QRS_Off(Peak_Ind)-QRS_On(Peak_Ind))<=(71/(1000/FS)))
                                wrong_R=1;
                            end
                            if(wrong_R==0)  
%                                 if(R_Reversed==1)
%                                     update=0;
%                                 end
                                if(update==1)
                                    if(R_Reversed==1)
                                        QRS_temp=QRS_temp.*-1;
                                    end
                                    QRS(Peak_Ind,:)=QRS_temp;
                                    QRS_Ind_temp=QRS_Ind_temp-QRS_On(Peak_Ind);
                                    QRS_Inds(Peak_Ind,:)=QRS_Ind_temp;
                                    %% Reconstructing the QRS structure
                                     if(Peak_Ind>QRS_NO+1)
%                                         Activity_Ratio=mean(Peak_activities(Peak_Ind-1-QRS_NO:Peak_Ind-1))/Peak_activities(Peak_Ind);
%                                     end
%                                     if(Peak_Ind>QRS_NO+1 && (Activity_Ratio>0.67 && Activity_Ratio<15))
                                         % Adaptive Learning Coefficient, Adapting the
                                         % Cooef
%                                          if(mean(Peak_activities(Peak_Ind-1-QRS_NO:Peak_Ind-1))*.8>Peak_activities(Peak_Ind))
%                                              Etha=Etha+0.05;
%                                              if(mean(Peak_activities(Peak_Ind-1-QRS_NO:Peak_Ind-1))*.8>Peak_activities(Peak_Ind))
%                                                  Etha=Etha+0.5;
%                                              end
%                                              Etha=min(Etha,1);
%                                          elseif(mean(Peak_activities(Peak_Ind-1-QRS_NO:Peak_Ind-1))<Peak_activities(Peak_Ind)*.8)
%                                              Etha=Etha-0.05;
%                                              if(mean(Peak_activities(Peak_Ind-1-QRS_NO:Peak_Ind-1))<Peak_activities(Peak_Ind)*.8)
%                                                   Etha=Etha-0.05;
%                                              end
%                                              Etha=max(Etha,0.1);
%                                          else
%                                              Etha=0.25;
%                                          end
 Etha=0.25;
                                         QRS_NOW=QRS(Peak_Ind-QRS_NO-1:Peak_Ind-1,:);
                                         QRS_Inds_NOW=QRS_Inds(Peak_Ind-QRS_NO-1:Peak_Ind-1,:);
    %                                         QRS_Now=ecg(QRS_On(Peak_Ind):QRS_Off(Peak_Ind));
    %                                         QRS_Now=resample(QRS_Now,Samples,length(QRS_Now))';
    %                                         QRS_Now=[zeros(1,10) QRS_Now*0.5 zeros(1,10)];                                        
    %                                         g=(1-Etha)*g + Etha*QRS_Now;
                                           try 
%                                                     [G_Vals,G_Inds,g]=Rebuild_Structure(g,G_Vals,G_Inds,QRS_NOW(end,:)-sigmean,QRS_Inds_NOW(end,:),Samples,Etha);
                                           catch
                                           end
                                          gs=[gs g ones(1,20)*g(end)];
                                     end
                                end
                                Peak_Ind=Peak_Ind+1;
                            end
                        end
                        wrong_R=0;
                  end
              end
            end
        elseif(not_changed > Inactivity && checked==0)
            Onset=0;
            checked=1;
        end
        j=j+1;
    end
    
end
%% trimming the unused allocation
% for Peak_Ind Datatypes
QRS_On=QRS_On(1:Peak_Ind-1); QRS_Off=QRS_Off(1:Peak_Ind-1); 
Q_inds=Q_inds(1:Peak_Ind-1); R_inds=R_inds(1:Peak_Ind-1); S_inds=S_inds(1:Peak_Ind-1); 
QRS=QRS(1:Peak_Ind-1,:); QRS_Inds=QRS_Inds(1:Peak_Ind-1,:); 
Peak_activities=Peak_activities(1:Peak_Ind-1,:);

%  figure(1);
%  plot(gs);

function Out=dilation(X,s)
    len=floor(length(s)/2);
    S=s';    
    for i=len+1:length(X)-len
        out(i)=max(X(i-len:i+len)+S);           
    end
    Out=[out zeros(1,len)];
return

function Out=erosion(X,s)
    len=floor(length(s)/2);
    S=s';
     S=flip(S);
    for i=len+1:length(X)-len
        out(i)=min(X(i-len:i+len)-S);
    end
    Out=[out zeros(1,len)];
return

function opened=open_ecg(x,s)
    opened=dilation(erosion(x,s)',s);
return
    
function closed=close_ecg(x,s)
    closed=erosion(dilation(x,s)',s);
return

function [VALS,POS,New_G]=Rebuild_Structure(G,G_Vals,G_Inds,QRS,Indecies,Samples,Etha)
    % Creating the rough estimate of the Structres positions and values.
    %Values
    global g_amplitude;
    global g_length;
     G_Vals=G_Vals*(1/g_amplitude);
    if(G_Vals==zeros(1,5)), Etha=1; end        
     On_Val=(1-Etha)*G_Vals(1) + Etha*mean(QRS(:,1)); 
     Q_Val=(1-Etha)*G_Vals(2) + Etha*mean(QRS(:,2));
     R_Val=(1-Etha)*G_Vals(3) + Etha*mean(QRS(:,3));
     S_Val=(1-Etha)*G_Vals(4) + Etha*mean(QRS(:,4));
     Off_Val=(1-Etha)*G_Vals(5) + Etha*mean(QRS(:,5));    
    %Pos
    Q_Pos=floor((1-Etha)*G_Inds(2)+ Etha*mean(Indecies(:,1)));
    R_Pos=floor((1-Etha)*G_Inds(3)+ Etha*mean(Indecies(:,2)));
    S_Pos=floor((1-Etha)*G_Inds(4)+ Etha*mean(Indecies(:,3)));
    Off_Pos=floor((1-Etha)*G_Inds(5)+ Etha*mean(Indecies(:,4)));
    
      VALS=[On_Val Q_Val R_Val S_Val Off_Val];
      VALS=VALS*g_amplitude;
     POS=round([0 Q_Pos R_Pos S_Pos Off_Pos]*.8); 
%            POS=G_Inds;
     if((Q_Val-On_Val)~=0 && Q_Pos~=0), New_G = interp1(POS,VALS,0:POS(end)/(Samples-1):POS(end));
     else New_G = interp1(POS(2:end),VALS(2:end),0:POS(end)/(Samples-1):POS(end));
     end
     New_G=New_G-sign(New_G(1))*New_G(1);
%        New_G=[zeros(1,10) New_G zeros(1,10)];

%       glen=min(5,round(length(New_G)*(g_length*G_Inds(4))));
%       New_G=resample(New_G,round(length(New_G)*(g_length*G_Inds(4))),length(New_G));


%       New_G=New_G*g_amplitude;
      
      
      
%      plot(New_G);
%      pause;


%     if((S_Pos-Q_Pos)<=Samples)
%         %Ploting the new QRS structure's schematic
%         VALS=[On_Val Q_Val R_Val S_Val Off_Val];
%         POS=[0 Q_Pos R_Pos S_Pos Off_Pos];
% 
%         %%Building the new Structure based on the Sampling rate of the ECG singal
%         New_G=[];
%         if((Q_Val-On_Val)~=0 && Q_Pos~=0), New_G=[New_G On_Val:sign(Q_Val-On_Val)/(Q_Pos):Q_Val;]; end
%         New_G=[New_G Q_Val:sign(R_Val-Q_Val)/(R_Pos-Q_Pos):R_Val;];
%         New_G=[New_G R_Val:sign(S_Val-R_Val)/(S_Pos-R_Pos):S_Val;];
%         New_G=[New_G S_Val:sign(Off_Val-S_Val)/(Off_Pos-S_Pos):Off_Val;];  
% 
%         New_G=resample(New_G,length(G),length(New_G));
%         New_G=New_G-sign(New_G(1))*New_G(1);
%     else
%         VALS=G_Vals;
%         POS=G_Inds;
%         New_G=G;
%     end

return