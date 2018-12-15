clc
clear all
close all

number = 5;
for trial=1:5
counter_MUSE = zeros(5,65);
counter_SMA = zeros(5,65);
counter_HW = zeros(5,65);
for num_index=1:number
T = 10000;
% num = [4:1:8];
num = [5:15:65];
num_users = num(num_index);
pkts = round(1*rand(1,num_users));
pkts_MUSE = pkts;
pkts_SMA = pkts;
% pkts_WMA = pkts;
pkts_HW = pkts;
colFlag_MUSE = zeros(1,num_users); %collision flag for MUSE
colFlag_SMA = zeros(1,num_users); %collision flag for SMA
% colFlag_WMA = zeros(1,num_users); %collision flag for WMA
colFlag_HW = zeros(1,num_users); %collision flag for HW
cw = 32;
backoff = round((cw-1)*rand(1,num_users));
backoff_MUSE = backoff;
backoff_SMA = backoff;
% backoff_WMA = backoff;
backoff_HW = backoff;
time_index_MUSE = zeros(2,T); %0=bf, 1=collision, 2=transmission occurs
time_index_SMA = zeros(2,T);
% time_index_WMA = zeros(2,T);
time_index_HW = zeros(2,T);
ExpTraffic = 1*rand(1,num_users);
ExpTraffic_SMA = ExpTraffic;
% ExpTraffic_WMA = ExpTraffic;
ExpTraffic_HW = ExpTraffic;
traffic_MUSE = zeros(num_users,T);
traffic_SMA = zeros(num_users,T);
% traffic_WMA = zeros(num_users,T);
traffic_HW = zeros(num_users,T);
succ = 0;
% Create "rates" at which packets will arrive, these stay constant for
% whole test? This seems like it would make more sense since we are doing 
% predictive modeling 
packet_generation_rates = 0.9 + (0.25)*rand(num_users,1); % U[0.75,1]

for time=1:T
    [trial, time]
    winner_MUSE = [];
    winner_SMA = [];
%     winner_WMA = [];
    winner_HW = [];
    groupMembers_MUSE = [];
    groupMembers_SMA = [];
%     groupMembers_WMA = [];
    groupMembers_HW = [];
    traffic_MUSE(:,time) = pkts_MUSE;
    traffic_SMA(:,time) = pkts_SMA;
%     traffic_WMA(:,time) = pkts_WMA;
    traffic_HW(:,time) = pkts_HW;
    % Add packets to each user as a geometric distribution
    for i = 1:length(packet_generation_rates)
        pkts_gen = geornd(packet_generation_rates(i));
        pkts_MUSE(i) = pkts_MUSE(i) + pkts_gen;
        pkts_SMA(i) = pkts_SMA(i) + pkts_gen;
%         pkts_WMA(i) = pkts_WMA(i) + pkts_gen;
        pkts_HW(i) = pkts_HW(i) + pkts_gen;
    end    
    %% MUSE
    if(ismember(0,backoff_MUSE)) % backoff to zero
        winner_MUSE = find(backoff_MUSE==0);
        if(length(winner_MUSE)>1) %collision
            colFlag_MUSE(winner_MUSE) = colFlag_MUSE(winner_MUSE) + 1; % collision times for clients
            time_index_MUSE(1,time) = 1; % index of collision time
        else
            % transmission starts
            time_index_MUSE(1,time) = 2; % index of transmission time
            [packets_MUSE,hasData_MUSE,groupMembers_MUSE] = baseline(num_users,winner_MUSE,pkts_MUSE);
            pkts_MUSE = packets_MUSE;
            time_index_MUSE(2,time) = sum(hasData_MUSE);
            
            for l=1:length(groupMembers_MUSE)
                counter_MUSE(num_index,groupMembers_MUSE(l)) = counter_MUSE(num_index,groupMembers_MUSE(l)) + 1;
            end
            
            colFlag_MUSE(winner_MUSE) = 0;
            colFlag_MUSE(groupMembers_MUSE) = 0;
        end
    end
%     winner
%     groupMembers_MUSE
%     backoff_MUSE
%     colFlag_MUSE
    backoff_MUSE = backoff_MUSE - 1;
    if(length(winner_MUSE)==1) % reselect bf after transmission
        backoff_MUSE(winner_MUSE) = round((cw-1)*rand(1,1));
        if(sum(groupMembers_MUSE>0)>0)
            backoff_MUSE(groupMembers_MUSE) = round((cw-1)*rand(1,sum(groupMembers_MUSE>0)));
        end
    elseif(length(winner_MUSE)>1 && sum(colFlag_MUSE(winner_MUSE))>0) % exponential bf due to collision
        for i=1:length(winner_MUSE)
            backoff_MUSE(winner_MUSE(i)) = round((cw*pow2(colFlag_MUSE(winner_MUSE(i)))-1)*rand(1,1));
        end
    else
        backoff_MUSE = backoff_MUSE;
    end
    
    %% SMA
    mode = 1;
    if(ismember(0,backoff_SMA)) % backoff to zero
        winner_SMA = find(backoff_SMA==0);
        if(length(winner_SMA)>1) %collision
            colFlag_SMA(winner_SMA) = colFlag_SMA(winner_SMA) + 1; % collision times for clients
            time_index_SMA(1,time) = 1; % index of collision time
        else
            % transmission starts
            time_index_SMA(1,time) = 2; % index of transmission time
            succ = sum(time_index_SMA(1,:)==2);
            [packets_SMA,hasData_SMA,groupMembers_SMA,TrafficBuffer_SMA] = user_selection(num_users,winner_SMA,pkts_SMA,mode,ExpTraffic_SMA,succ);
            pkts_SMA = packets_SMA;
            ExpTraffic_SMA = TrafficBuffer_SMA;
            time_index_SMA(2,time) = sum(hasData_SMA);
            
            for l=1:length(groupMembers_SMA)
                counter_SMA(num_index,groupMembers_SMA(l)) = counter_SMA(num_index,groupMembers_SMA(l)) + 1;
            end
            
            colFlag_SMA(winner_SMA) = 0;
            colFlag_SMA(groupMembers_SMA) = 0;
        end
    end
%     winner_SMA
%     groupMembers_SMA
%     backoff_SMA
%     colFlag_SMA
    backoff_SMA = backoff_SMA - 1;
    if(length(winner_SMA)==1) % reselect bf after transmission
        backoff_SMA(winner_SMA) = round((cw-1)*rand(1,1));
        if(sum(groupMembers_SMA>0)>0)
            backoff_SMA(groupMembers_SMA) = round((cw-1)*rand(1,sum(groupMembers_SMA>0)));
        end
    elseif(length(winner_SMA)>1 && sum(colFlag_SMA(winner_SMA))>0) % exponential bf due to collision
        for i=1:length(winner_SMA)
            backoff_SMA(winner_SMA(i)) = round((cw*pow2(colFlag_SMA(winner_SMA(i)))-1)*rand(1,1));
        end
    else
        backoff_SMA = backoff_SMA;
    end  
    
    %% WMA
%     mode = 2;
%     if(ismember(0,backoff_WMA)) % backoff to zero
%         winner_WMA = find(backoff_WMA==0);
%         if(length(winner_WMA)>1) %collision
%             colFlag_WMA(winner_WMA) = colFlag_WMA(winner_WMA) + 1; % collision times for clients
%             time_index_WMA(1,time) = 1; % index of collision time
%         else
%             % transmission starts
%             time_index_WMA(1,time) = 2; % index of transmission time
%             [packets_WMA,hasData_WMA,groupMembers_WMA,TrafficBuffer_WMA] = user_selection(num_users,winner_WMA,pkts_WMA,mode,ExpTraffic_WMA,succ);
%             pkts_WMA = packets_WMA;
%             ExpTraffic_WMA = TrafficBuffer_WMA;
%             time_index_WMA(2,time) = sum(hasData_WMA);
%             
%             colFlag_WMA(winner_WMA) = 0;
%             colFlag_WMA(groupMembers_WMA) = 0;
%         end
%     end
% %     winner_WMA
% %     groupMembers_WMA
% %     backoff_WMA
% %     colFlag_WMA
%     backoff_WMA = backoff_WMA - 1;
%     if(length(winner_WMA)==1) % reselect bf after transmission
%         backoff_WMA(winner_WMA) = round((cw-1)*rand(1,1));
%         if(sum(groupMembers_WMA>0)>0)
%             backoff_WMA(groupMembers_WMA) = round((cw-1)*rand(1,sum(groupMembers_WMA>0)));
%         end
%     elseif(length(winner_WMA)>1 && sum(colFlag_WMA(winner_WMA))>0) % exponential bf due to collision
%         for i=1:length(winner_WMA)
%             backoff_WMA(winner_WMA(i)) = round((cw*pow2(colFlag_WMA(winner_WMA(i)))-1)*rand(1,1));
%         end
%     else
%         backoff_WMA = backoff_WMA;
%     end  
    
    %% HW
    mode = 3;
    if(ismember(0,backoff_HW)) % backoff to zero
        winner_HW = find(backoff_HW==0);
        if(length(winner_HW)>1) %collision
            colFlag_HW(winner_HW) = colFlag_HW(winner_HW) + 1; % collision times for clients
            time_index_HW(1,time) = 1; % index of collision time
        else
            % transmission starts
            time_index_HW(1,time) = 2; % index of transmission time
            [packets_HW,hasData_HW,groupMembers_HW,TrafficBuffer_HW] = user_selection(num_users,winner_HW,pkts_HW,mode,ExpTraffic_HW,succ);
            pkts_HW = packets_HW;
            ExpTraffic_HW = TrafficBuffer_HW;
            time_index_HW(2,time) = sum(hasData_HW);
            
            for l=1:length(groupMembers_HW)
                counter_HW(num_index,groupMembers_HW(l)) = counter_HW(num_index,groupMembers_HW(l)) + 1;
            end
            
            colFlag_HW(winner_HW) = 0;
            colFlag_HW(groupMembers_HW) = 0;
        end
    end
%     winner_HW
%     groupMembers_HW
%     backoff_HW
%     colFlag_HW
    backoff_HW = backoff_HW - 1;
    if(length(winner_HW)==1) % reselect bf after transmission
        backoff_HW(winner_HW) = round((cw-1)*rand(1,1));
        if(sum(groupMembers_HW>0)>0)
            backoff_HW(groupMembers_HW) = round((cw-1)*rand(1,sum(groupMembers_HW>0)));
        end
    elseif(length(winner_HW)>1 && sum(colFlag_HW(winner_HW))>0) % exponential bf due to collision
        for i=1:length(winner_HW)
            backoff_HW(winner_HW(i)) = round((cw*pow2(colFlag_HW(winner_HW(i)))-1)*rand(1,1));
        end
    else
        backoff_HW = backoff_HW;
    end
end

t_DIFS = 34*1e-06;
t_SIFS = 16*1e-06;
t_BF = 9*1e-06;
t_ACK = 24*1e-06;
t_TRIG = 24*1e-06;
t_DATA = 248*1e-06;

%% MUSE
throughput_MUSE(trial,num_index) = 0;
success_MUSE = find(time_index_MUSE(1,:)==2);
for muse=1:length(success_MUSE)-1
   Data_MUSE = time_index_MUSE(2,success_MUSE(muse))*t_DATA;
   OH_MUSE = time_index_MUSE(2,success_MUSE(muse))*t_DATA+t_TRIG+5*t_SIFS+4*t_ACK+t_DIFS;
   if(success_MUSE(muse+1) == success_MUSE(muse)+1)
       Col_MUSE = 0;
       BF_MUSE = 0;
   else
       Col_MUSE = sum(time_index_MUSE(1,success_MUSE(muse)+1:success_MUSE(muse+1)-1)==1)*t_DIFS;
       BF_MUSE = sum(time_index_MUSE(1,success_MUSE(muse)+1:success_MUSE(muse+1)-1)==0)*t_BF;
   end
   throughput_MUSE(trial,num_index) = throughput_MUSE(trial,num_index) + Data_MUSE/(OH_MUSE+Col_MUSE+BF_MUSE);
end
throughput_MUSE(trial,num_index) = throughput_MUSE(trial,num_index)/(length(success_MUSE)-1);

%% SMA
throughput_SMA(trial,num_index) = 0;
success_SMA = find(time_index_SMA(1,:)==2);
for sma=1:length(success_SMA)-1
   Data_SMA = time_index_SMA(2,success_SMA(sma))*t_DATA;
   OH_SMA = time_index_SMA(2,success_SMA(sma))*t_DATA+2*t_TRIG+6*t_SIFS+4*t_ACK+t_DIFS;
   if(success_SMA(sma+1) == success_SMA(sma)+1)
       Col_SMA = 0;
       BF_SMA = 0;
   else
       Col_SMA = sum(time_index_SMA(1,success_SMA(sma)+1:success_SMA(sma+1)-1)==1)*t_DIFS;
       BF_SMA = sum(time_index_SMA(1,success_SMA(sma)+1:success_SMA(sma+1)-1)==0)*t_BF;
   end
   throughput_SMA(trial,num_index) = throughput_SMA(trial,num_index) + Data_SMA/(OH_SMA+Col_SMA+BF_SMA);
end
throughput_SMA(trial,num_index) = throughput_SMA(trial,num_index)/(length(success_SMA)-1);

%% HW
throughput_HW(trial,num_index) = 0;
success_HW = find(time_index_HW(1,:)==2);
for hw=1:length(success_HW)-1
   Data_HW = time_index_HW(2,success_HW(hw))*t_DATA;
   OH_HW = time_index_HW(2,success_HW(hw))*t_DATA+2*t_TRIG+6*t_SIFS+4*t_ACK+t_DIFS;
   if(success_HW(hw+1) == success_HW(hw)+1)
       Col_HW = 0;
       BF_HW = 0;
   else
       Col_HW = sum(time_index_HW(1,success_HW(hw)+1:success_HW(hw+1)-1)==1)*t_DIFS;
       BF_HW = sum(time_index_HW(1,success_HW(hw)+1:success_HW(hw+1)-1)==0)*t_BF;
   end
   throughput_HW(trial,num_index) = throughput_HW(trial,num_index) + Data_HW/(OH_HW+Col_HW+BF_HW);
end
throughput_HW(trial,num_index) = throughput_HW(trial,num_index)/(length(success_HW)-1);


% % total time of data
% Data_MUSE = sum(time_index_MUSE(2,time_index_MUSE(2,:)>0))*t_DATA;
% Data_SMA = sum(time_index_SMA(2,time_index_SMA(2,:)>0))*t_DATA;
% % Data_WMA = sum(time_index_WMA(2,time_index_WMA(2,:)>0))*t_DATA;
% Data_HW = sum(time_index_HW(2,time_index_HW(2,:)>0))*t_DATA;
% % total time of data including overhead
% OH_MUSE = sum(time_index_MUSE(2,time_index_MUSE(2,:)>0))*t_DATA+t_TRIG+5*t_SIFS+4*t_ACK+t_DIFS;
% OH_SMA = sum(time_index_SMA(2,time_index_SMA(2,:)>0))*t_DATA+2*t_TRIG+6*t_SIFS+4*t_ACK+t_DIFS;
% % OH_WMA = sum(time_index_WMA(2,time_index_WMA(2,:)>0))*t_DATA+2*t_TRIG+6*t_SIFS+4*t_ACK+t_DIFS;
% OH_HW = sum(time_index_HW(2,time_index_HW(2,:)>0))*t_DATA+2*t_TRIG+6*t_SIFS+4*t_ACK+t_DIFS;
% % collision time
% Col_MUSE = sum(time_index_MUSE(1,time_index_MUSE(1,:)==1))*t_DIFS;
% Col_SMA = sum(time_index_SMA(1,time_index_SMA(1,:)==1))*t_DIFS;
% % Col_WMA = sum(time_index_WMA(1,time_index_WMA(1,:)==1))*t_DIFS;
% Col_HW = sum(time_index_HW(1,time_index_HW(1,:)==1))*t_DIFS;
% % backoff time
% BF_MUSE = sum(time_index_MUSE(1,time_index_MUSE(1,:)==0))*t_BF;
% BF_SMA = sum(time_index_SMA(1,time_index_SMA(1,:)==0))*t_BF;
% % BF_WMA = sum(time_index_WMA(1,time_index_WMA(1,:)==0))*t_BF;
% BF_HW = sum(time_index_HW(1,time_index_HW(1,:)==0))*t_BF;

% TH_MUSE = Data_MUSE/(OH_MUSE+Col_MUSE+BF_MUSE)
% TH_SMA = Data_SMA/(OH_SMA+Col_SMA+BF_SMA)
% TH_WMA = Data_WMA/(OH_WMA+Col_WMA+BF_WMA)
% TH_HW = Data_HW/(OH_HW+Col_HW+BF_HW)
% TH_MUSE(trial,num_index) = Data_MUSE/(OH_MUSE+Col_MUSE+BF_MUSE)
% TH_SMA(trial,num_index) = Data_SMA/(OH_SMA+Col_SMA+BF_SMA)
% TH_WMA(trial,num_index) = Data_WMA/(OH_WMA+Col_WMA+BF_WMA)
% TH_HW(trial,num_index) = Data_HW/(OH_HW+Col_HW+BF_HW)
end
end


%Jain's fairness index
for i=1:length(num)
    Jain_MUSE(i) = sum(counter_MUSE(i,1:num(i)))^2/(num(i)*sum(counter_MUSE(i,1:num(i)).^2));
    Jain_SMA(i) = sum(counter_SMA(i,1:num(i)))^2/(num(i)*sum(counter_SMA(i,1:num(i)).^2));
    Jain_HW(i) = sum(counter_HW(i,1:num(i)))^2/(num(i)*sum(counter_HW(i,1:num(i)).^2));
end

figure
plot(mean(throughput_MUSE),'-o');hold on;
plot(mean(throughput_SMA),'-x');hold on;
plot(mean(throughput_HW),'-p');
legend('MUSE','UserSelect-MA','UserSelect-HW');
ylabel('Normalized Throughput(%)');
xlabel('Number of Users');
title('Throughput P=0.9');
xx = num;
set(gca,'XTick',1:length(xx),'fontsize',18);
set(gca,'XTickLabel',xx);

figure
plot(Jain_MUSE,'-o');hold on;
plot(Jain_SMA,'-x');hold on;
plot(Jain_HW,'-p');
legend('MUSE','UserSelect-MA','UserSelect-HW');
% ylabel('Normalized Phy Throughput(%)');
xlabel('Number of Users');
title('0.9 Packet Generation Rates - Jains Fairness');
xx = num;
set(gca,'XTick',1:length(xx),'fontsize',18);
set(gca,'XTickLabel',xx);


